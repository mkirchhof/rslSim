# Bayesian Network based probabilistic Rule Stacking Learner
# Author: michael.kirchhof@udo.edu
# Created: 07.01.2021
# Version: 0.5.1 "Hide It"

# Dependencies: (not loaded into namespace due to style guide)
# library(bnlearn) # for constructing bayesian networks
# library(gRain) # for exact inference on bayesian networks
# library(MASS) # for ginv()
source("norn.R")


# createRSL - creates an empty Rule Stacking Learner
createRSL <- function(){
  # labels: list of lists. Each list represents one XOR-related group of labels,
  #         each entry contains the label names, the internal ID in the BN 
  #         and the prior
  # classifiers: list of lists. Each list contains a classifier's name, its
  #              internal ID in the BN, its confusion matrix and its prior
  # rules: dataframe containing a rules's name, probability, and the IDs of 
  #        the nodes representing the rule and its auxiliary node in the BN
  # bayesNet: the bnlearn object representing the network (always updated when
  #           the rsl is updated)
  # compiledNet: the grain object representing the network (only compiled when
  #              needed)
  # needsCompilation: logical, flags if the compiled net is in sync with rsl
  bn <- bnlearn::empty.graph(nodes = "placeholder")
  cpt <- list(placeholder = array(c(1, 0), dim = 2, dimnames = list(placeholder = c("true", "false"))))
  bn <- bnlearn::custom.fit(bn, cpt)
  rsl <- list(labels = list(),
              classifiers = list(),
              rules = data.frame("name" = character(0), 
                                 "prob" = numeric(0), 
                                 "ruleID" = character(0),
                                 "auxID" = character(0),
                                 stringsAsFactors = FALSE),
              bayesNet = bn,
              norn = create.norn(),
              compiledNet = NULL,
              needsCompilation = TRUE)
  class(rsl) <- "rsl"
  
  return(rsl)
}


# .normalizeRule - removes whitespaces from a rule
.normalizeRule <- function(rule){
  return(gsub("[[:blank:]]", "", rule))
}


# addRule - adds a probabilistic rule to a rsl
# Input:
#  rsl - an rsl object
#  rule - a character string describing the rule in the form "A, B <- !C, D"
#         which corresponds to "((not C) and D) -> (A or B)", where A, B, C and 
#         D are names of classes in the network
#  prob - numeric in [0, 1], the marginal probability of the rule, or NA to
#         learn optimal prob given data
# Output:
#  the updated rsl object
addRule <- function(rsl, rule, prob = 0.9){
  # TODO: Add argument data to automatically learn the prob
  
  if(.ruleAlreadyExists(rsl, rule)){
    warning(paste0("Rule ", rule, " already exists in the rsl. ",
                   "Please adjust its prob instead of adding it multiple times. Skipping this rule."))
    return(rsl)
  }
  # TODO: Add more type checks
  
  # Find label nodes associated to the rule
  rule <- .normalizeRule(rule)
  labels <- .decomposeRule(rule)
  rID <- .getNewRuleID(rsl)
  aID <- .getNewAuxID(rsl)
  
  # Add to rsl$rules
  rsl$rules <- rbind(rsl$rules, 
                     data.frame(name = rule, prob = prob, ruleID = rID, auxID = aID,
                                stringsAsFactors = FALSE))
  
  # save for each label node the labels that fulfill the rule
  # So we get a list like list(L1 = c("cat", "mouse"), L2 = c("hungry"), L3 = c())
  allowedStates <- list()
  for(i in seq(along = labels)){
    labelID <- .labelToID(rsl, names(labels)[i])
    if(!labelID %in% names(allowedStates)){
      # if we do not have this label node yet, add it
      allowedStates[[labelID]] <- rsl$labels[[labelID]]$names
    }
    if(labels[i] == TRUE){ # The rare occasion where "== TRUE" adds to understanding 
      allowedStates[[labelID]] <- intersect(allowedStates[[labelID]], names(labels)[i])
    } else {
      allowedStates[[labelID]] <- setdiff(allowedStates[[labelID]], names(labels)[i])
    }
  }
  
  # build rule CPT
  # for every combination of variables, check if they are %in% their allowedStates
  dimlist <- lapply(rsl$labels[names(allowedStates)], function(x) x$names)
  ruleStates <- list(c("fulfilled", "not_fulfilled"))
  names(ruleStates) <- rID
  allCombs <- expand.grid(dimlist)
  isFulfilled <- matrix(FALSE, ncol = ncol(allCombs), nrow = nrow(allCombs))
  for(labelIndex in seq(ncol(allCombs))){
    isFulfilled[, labelIndex] <- allCombs[, labelIndex] %in% allowedStates[[labelIndex]]
  }
  # We have a disjunctive form, so check row-wise if any label is ok
  ruleFulfilled <- rowSums(isFulfilled) > 0
  rProbTable <- array(as.numeric(c(t(matrix(c(ruleFulfilled, !ruleFulfilled), ncol = 2)))), 
                      dim = sapply(c(ruleStates, dimlist), length), 
                      dimnames = c(ruleStates, dimlist))
  
  # Build aux CPT
  auxStates <- list(c("fulfilled", "not_fulfilled"))
  names(auxStates) <- aID
  aProbTable <- array(c(prob, 1 - prob, 1 - prob, prob),
                      dim = c(2, 2), dimnames = c(auxStates, ruleStates))
  
  # add rule and aux to rsl$bayesNet
  tables <- lapply(rsl$bayesNet, "[[", "prob")
  tables[[rID]] <- rProbTable
  tables[[aID]] <- aProbTable
  nodes <- c(bnlearn::nodes(rsl$bayesNet), rID, aID)
  arcs <- rbind(bnlearn::arcs(rsl$bayesNet), 
                cbind(names(allowedStates), rID), 
                c(rID, aID))
  rsl$bayesNet <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::arcs(rsl$bayesNet) <- arcs
  rsl$bayesNet <- bnlearn::custom.fit(rsl$bayesNet, tables)
  rsl$needsCompilation <- TRUE
  
  # add rule and aux to rsl$norn
  # cast allowedStates to a list of inhProbs
  inhProbs <- getLabels(rsl)
  inhProbs <- inhProbs[names(allowedStates)]
  for(i in seq(along = inhProbs)){
    inhProbs[[i]] <- ifelse(inhProbs[[i]] %in% allowedStates[[i]], 0, 1)
  }
  rsl$norn <- .addRule.norn(rsl$norn, inhProbs, rID, aID, prob)
  
  return(rsl)
}


# .ruleAlreadyExists - checks if a rule already exists in the rsl
.ruleAlreadyExists <- function(rsl, rule){
  rule <- .normalizeRule(rule)
  return(any(rsl$rules$name == rule))
}


# .labelToID - returns the label node a label or subset of labels is associated to
.labelToID <- function(rsl, labels){
  isLabel <- sapply(rsl$labels, function(x) all(labels %in% x$names))
  if(length(isLabel) > 0 && sum(isLabel) == 1){
    return(rsl$labels[[which(isLabel)]]$id)
  } else {
    return(NULL)
  }
}


# .isValidRule - returns if a character is in the form "A, B <- !C, D" etc.
.isValidRule <- function(rule){
  rule <- .normalizeRule(rule)
  # rule has correct scheme
  isValid <- grepl("^(!?[^!,<>]+,)*(!?[^!,<>]+)?<-(!?[^!,<>]+,)*(!?[^!,<>]+)?$", rule)
  # rule has at least one variable somewhere
  isValid <- isValid && grepl("^(!?[^!,<>]+){1}(,!?[^!,<>]+)*<-|<-(!?[^!,<>]+){1}(,!?[^!,<>]+)*$", rule)
  
  return(isValid)
}


# .decomposeRule - reads a rule and returns a named boolean vector representing
#                  that rule as a disjunction. The names are the variables and
#                  the boolean value is the true/false state of the variable.
.decomposeRule <- function(rule){
  if(!.isValidRule(rule)){
    stop("Rule ", rule, " is not in a valid format.")
  }
  
  rule <- .normalizeRule(rule)
  head <- gsub("^([^<]*)<-.*$", "\\1", rule)
  body <- gsub("^([^<]*)<-(.*)$", "\\2", rule)
  head <- .textToVariables(head)
  body <- .textToVariables(body)
  variables <- c(head, !body) # A <- B is eq. to "A or not B"
  
  return(variables)
}


# .textToVariables - turns a character in the form "A,!B,C" into a named boolean
#                    vector like c(A = TRUE, B = FALSE, C = TRUE)
.textToVariables <- function(text){
  entries <- strsplit(text, ",")[[1]]
  variables <- !grepl("^!", entries)
  names(variables) <- gsub("^!", "", entries)
  
  return(variables)
}


# .optimizeRuleProb - uses a line search to find a rule probability that
#                     optimizes the a-posteriori likelihood of the given data
.optimizeRuleProb <- function(rsl, rule, data){
  # TODO: Implement (medium priority)
}


# .addJointRule - adds a joint rule (learned in the learnRules intermediate step)
#                 CAUTION: Only for internal testing purposes!
# Input: 
#  rsl - an rsl object
#  weights - the weights of the joint rule in standard order 
#           (see .generateStandardOrder())
# Output:
#  the rsl object with the added joint rule
.addJointRule <- function(rsl, weights){
  # get IDs for the joint rule
  rID <- .getNewRuleID(rsl)
  aID <- .getNewAuxID(rsl)
  
  # Add to rsl$rules
  rsl$rules <- rbind(rsl$rules, 
                     data.frame(name = "jointRule", prob = 1, ruleID = rID, auxID = aID,
                                stringsAsFactors = FALSE))
  
  # Build weights into a CPT
  dimlist <- lapply(rsl$labels, function(x) x$names)
  ruleStates <- list(c("fulfilled", "not_fulfilled"))
  names(ruleStates) <- rID
  rProbTable <- array(as.numeric(c(t(matrix(c(weights, 1 - weights), ncol = 2)))), 
                      dim = sapply(c(ruleStates, dimlist), length), 
                      dimnames = c(ruleStates, dimlist))
  
  # Build aux CPT
  auxStates <- list(c("fulfilled", "not_fulfilled"))
  names(auxStates) <- aID
  aProbTable <- array(c(1, 0, 0, 1),
                      dim = c(2, 2), dimnames = c(auxStates, ruleStates))
  
  # add rule and aux to rsl$bayesNet
  tables <- lapply(rsl$bayesNet, "[[", "prob")
  tables[[rID]] <- rProbTable
  tables[[aID]] <- aProbTable
  nodes <- c(bnlearn::nodes(rsl$bayesNet), rID, aID)
  arcs <- rbind(bnlearn::arcs(rsl$bayesNet), 
                cbind(names(getLabels(rsl)), rID), 
                c(rID, aID))
  rsl$bayesNet <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::arcs(rsl$bayesNet) <- arcs
  rsl$bayesNet <- bnlearn::custom.fit(rsl$bayesNet, tables)
  rsl$needsCompilation <- TRUE
  
  return(rsl)
}


# .preprocessInhProbs - takes a named vector of inhibition probs and turns it into
#                       a list with all inhibition probs per label node
.preprocessInhProbs <- function(rsl, probs){
  if(length(probs) == 0){
    return(list())
  }
  
  # match probs to labels (https://stackoverflow.com/a/51298361)
  labelIDs <- sapply(names(probs), .labelToID, rsl = rsl)
  probList <- split.default(probs, labelIDs)
  relLabels <- .getAllLabelNodes(rsl)
  relLabels <- relLabels[relLabels %in% labelIDs]
  probList <- probList[match(relLabels, names(probList))]
  
  # check if the nodes contain all of their labels
  # Add labels and re-order if necessary
  for(i in seq(along = probList)){
    allLabels <- .IDtoLabels(rsl, names(probList)[i])
    order <- match(allLabels, names(probList[[i]]))
    if(any(is.na(order))){
      # add missing labels
      missing <- setdiff(allLabels, names(probList[[i]]))
      missingProbs <- rep(1, length(missing))
      names(missingProbs) <- missing
      probList[[i]] <- c(probList[[i]], missingProbs)
      order <- match(allLabels, names(probList[[i]]))
    }
    probList[[i]] <- probList[[i]][order]
  }
  
  return(probList)
}


# .buildNoisyORCPT - builds a noisy-OR cpt given a list of vectors that contain
#                    the inhibition probs of each label
.buildNoisyORCPT <- function(rID, probs){
  nParents <- length(probs)
  
  # Create the dimlist
  dimlist <- list()
  dimlist[[1]] <- c("fulfilled", "not_fulfilled")
  names(dimlist)[1] <- rID
  for(i in seq(along = probs)){
    dimlist[[i + 1]] <- names(probs[[i]])
    names(dimlist)[i + 1] <- names(probs)[i]
  }
  
  if(nParents > 0){
    # create CPT via noisy-or: P(R = 0 | L1, ... Ln) = prod((inhProb_{l_i})^{l_i})
    cpt <- try(array(c(0, 1), dim = c(2, sapply(probs, length)), dimnames = dimlist))
    combs <- expand.grid(probs)
    probFalse <- exp(rowSums(log(combs)))
    cpt <- array(c(rbind(1 - probFalse, probFalse)), dim = c(2, sapply(probs, length)), dimnames = dimlist)
  } else if(nParents == 0){
    # If a noisy-or node has no parents, it is never activated. 
    # Thus add an almost-zero Prior for P(R) (to avoid numeric problems)
    cpt <- try(array(c(0.0001, 0.9999), dim = 2, dimnames = dimlist))
  }
  
  return(cpt)
}


# .addNoisyOR - adds a noisy-OR node to the rsl
# Input:
#  rsl - an rsl object
#  inhProbs - a named numeric vector, where each name gives a label and the 
#             number gives it inhibition probability (has to be in [0, 1])
#  prob - the rule probability
# Output:
#  the updated rsl
.addNoisyOR <- function(rsl, inhProbs, prob = 1){
  # Add to rsl$rules
  name <- paste0(names(inhProbs), " (", 1 - inhProbs, ")", collapse = " | ")
  rID <- .getNewRuleID(rsl)
  aID <- .getNewAuxID(rsl)
  rsl$rules <- rbind(rsl$rules, 
                     data.frame(name = name, prob = 1, ruleID = rID, auxID = aID,
                                stringsAsFactors = FALSE))
  
  # build noisy or prob Table
  probList <- .preprocessInhProbs(rsl, inhProbs)
  rProbTable <- .buildNoisyORCPT(rID, probList)
  
  # Build aux CPT
  auxStates <- list(c("fulfilled", "not_fulfilled"))
  names(auxStates) <- aID
  ruleStates <- list(c("fulfilled", "not_fulfilled"))
  names(ruleStates) <- rID
  aProbTable <- array(c(prob, 1 - prob, 1 - prob, prob),
                      dim = c(2, 2), dimnames = c(auxStates, ruleStates))
  
  # add rule and aux to rsl$bayesNet
  tables <- lapply(rsl$bayesNet, "[[", "prob")
  tables[[rID]] <- rProbTable
  tables[[aID]] <- aProbTable
  nodes <- c(bnlearn::nodes(rsl$bayesNet), rID, aID)
  arcs <- rbind(bnlearn::arcs(rsl$bayesNet), 
                cbind(names(probList), rep(rID, length(probList))), 
                c(rID, aID))
  rsl$bayesNet <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::arcs(rsl$bayesNet) <- arcs
  rsl$bayesNet <- bnlearn::custom.fit(rsl$bayesNet, tables)
  rsl$needsCompilation <- TRUE
  
  # add rule and aux to rsl$norn
  rsl$norn <- .addRule.norn(rsl$norn, probList, rID, aID, prob)
  
  return(rsl)
}


# .addAllNoisyOR - takes a matrix of inhibition probabilities for several 
#                  noisy-or rules, preprocesses them (throws out unnecessary ones)
#                  and adds them all to the rsl
# Input:
#  rsl - an rsl object
#  inhProbs - a matrix of inhibition probabilities (in [0, 1]) 
#             where each row is a rule and each column is a label. Columns need
#             to be named
#  ruleProbs - a numeric containing the rule prob for each rule. If only one prob
#              is given, it is used for all rules.
# Output:
#  the updated rsl object
.addAllNoisyOR <- function(rsl, inhProbs, ruleProbs = 1){
  if(nrow(inhProbs) == 0){
    return(rsl)
  }
  
  if(length(ruleProbs) == 1){
    ruleProbs <- rep(ruleProbs, nrow(inhProbs))
  }
  
  for(i in seq(nrow(inhProbs))){
    probs <- inhProbs[i, ]
    probs <- probs[probs != 1]
    rsl <- .addNoisyOR(rsl, probs, prob = ruleProbs[i])
  }
  
  return(rsl)
}


# removeRule - removes a rule with a given rule ID (something like R1, R2, ...)
#              from an rsl
removeRule <- function(rsl, rID){
  if(!rID %in% .getAllRules(rsl)){
    warning(paste("Rule", rID, "is not in the RSL. Skipping its removal."))
    return(rsl)
  }
  
  aID <- .ruleToAux(rsl, rID)
  
  # remove rule from rsl$rules
  rsl$rules <- rsl$rules[rsl$rules$ruleID != rID, ]
  
  # remove rule and aux from rsl$bayesNet
  tables <- lapply(rsl$bayesNet, "[[", "prob")
  tables[[rID]] <- NULL
  tables[[aID]] <- NULL
  nodes <- setdiff(bnlearn::nodes(rsl$bayesNet), c(rID, aID))
  arcs <- bnlearn::arcs(rsl$bayesNet)
  arcs <- arcs[!(arcs[, 1] %in% c(rID, aID) | arcs[, 2] %in% c(rID, aID)), ]
  rsl$bayesNet <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::arcs(rsl$bayesNet) <- arcs
  rsl$bayesNet <- bnlearn::custom.fit(rsl$bayesNet, tables)
  rsl$needsCompilation <- TRUE
  
  # remove rule from rsl$norn
  rsl$norn <- .removeRule.norn(rsl$norn, rID)
  
  return(rsl)
}


# .removeAllRules - removes all rules from an rsl
.removeAllRules <- function(rsl){
  for(rule in .getAllRules(rsl)){
    rsl <- removeRule(rsl, rule)
  }
  
  return(rsl)
}


# getRules - returns a dataframe with information about all rules
getRules <- function(rsl){
  return(rsl$rules[, c("name", "prob")])
}


.coerceToMultilabel <- function(labels){
  if(length(labels) == 1){
    warning(paste0("label seems to be binary. Coerced it to multilabel by adding ",
                   "not_", labels, " as second label."))
    # Treat as a binary label and add a "not" case
    labels <- c(labels, paste0("not_", labels))
  }
  
  return(labels)
}


# addLabels - adds an XOR-related group of labels (or a binary label)
addLabels <- function(rsl, labels, prior = NA){
  if(!is.character(labels) || length(labels) < 1){
    stop("Labels vector is invalid.")
  }
  if(!is.na(prior) && (!is.numeric(prior) || length(prior) != length(labels) || sum(prior) != 1)){
    stop("Prior vector is invaid.")
  }
  if(.labelsAlreadyExist(rsl, labels)){
    stop("Some of the labels already exist in the rsl. Please choose unique names.")
    # Note: If you want to change the code to support non-unique labels, 
    # code for learnRules() has to be changed aswell.
  }
  if(any(duplicated(labels))){
    stop("Label names are not unique.")
  }
  if(any(grepl("[!<>,]", labels))){
    stop("Label names must not contain the following characters: <>!,")
  }
  
  labels <- .coerceToMultilabel(labels)
  
  if(length(prior) == 1 && !is.na(prior)){
    prior <- c(prior, 1 - prior)
  }
  if(length(prior) == 1 && is.na(prior)){
    prior <- rep(1/length(labels), length(labels))
  }
  
  # Add to rsl$labels
  id <- .getNewLabelID(rsl)
  rsl$labels[[id]] <- list(names = labels, 
                           id = id,
                           prior = prior)
  
  # build the probTable for the label node
  dimlist <- list(labels)
  names(dimlist)[1] <- id
  probTable <- array(prior, dim = length(labels), dimnames = dimlist)
  
  # rebuild the rsl$bayesNet
  tables <- lapply(rsl$bayesNet, "[[", "prob")
  tables[[id]] <- probTable
  nodes <- c(bnlearn::nodes(rsl$bayesNet), id)
  arcs <- bnlearn::arcs(rsl$bayesNet)
  rsl$bayesNet <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::arcs(rsl$bayesNet) <- arcs
  rsl$bayesNet <- bnlearn::custom.fit(rsl$bayesNet, tables)
  rsl$needsCompilation <- TRUE
  
  # add to rsl$norn
  rsl$norn <- .addLabel.norn(rsl$norn, id, prior)
  
  return(rsl)
}


# .labelsAlreadyExist - checks if any of the labels in a given vector already
#                       exists in the rsl
.labelsAlreadyExist <- function(rsl, labels){
  isInLabelSet <- sapply(rsl$labels, function(x) any(labels %in% x$names))
  
  return(any(isInLabelSet))
}


# .classifierAlreadyExists - checks if there is an existing classifier called name
.classifierAlreadyExists <- function(rsl, name){
  hasEqualName <- sapply(rsl$classifiers, function(x) x$name == name)
  
  return(any(hasEqualName))
}


# removeLabels - removes a group of labels from the rsl
removeLabels <- function(rsl, labels){
  # TODO: Implement (low priority)
}


# getLabels - returns a list of all XOR-related label groups
getLabels <- function(rsl){
  return(lapply(rsl$labels, function(x) x$names))
}


# addClassifier - adds a classifier that connects to one XOR related group of 
#                 labels
# Input:
#  rsl - an rsl object
#  name - character, name of the classifier
#  labels - character, labels the classifier gives an output for
#  confusionMatrix - numeric matrix of size labels X labels, if NULL, will be 
#                    constructed using the accuracy (each row is an actual
#                    label, each column a prediction). Can either give the raw
#                    confusion matrix or a normalized version where columns
#                    sum to 1.
#  accuracy - numeric in [0, 1], used to calculate confusionMatrix if it is not
#             given. If NULL, standard value of 0.9 is used
#  prior - numeric, the priors of the actual labels
# Output:
#  updated rsl object
addClassifier <- function(rsl, name, labels, confusionMatrix = NULL, 
                          accuracy = NULL, prior = NA){
  if(.classifierAlreadyExists(rsl, name)){
    stop("name is already given to a classifier. Please select a new one.")
  }
  
  labels <- .coerceToMultilabel(labels)
  
  if(length(prior) == 1 && is.na(prior)){
    prior <- rep(1/length(labels), length(labels))
  }
  if(length(prior) == 1 && !is.na(prior)){
    # TODO: This could be dangerous if the labels are given by the user in a 
    # non-intuitive order. Maybe we should force the user to give all priors.
    prior <- c(prior, 1 - prior)
  }
  
  # Verify confusion matrix
  if(is.null(confusionMatrix)){
    # Construct confusion matrix from accuracy
    if(is.null(accuracy)){
      accuracy <- 0.95
    }
    confusionMatrix <- matrix((1 - accuracy) / (length(labels) - 1), 
                              ncol = length(labels), nrow = length(labels))
    diag(confusionMatrix) <- accuracy
  } else {
    if(!is.matrix(confusionMatrix)){
      stop("confusionMatrix is not a matrix.")
    }
    if(!is.numeric(confusionMatrix)){
      stop("confusionMatrix is not numeric.")
    }
    if(!all(dim(confusionMatrix) == c(length(labels), length(labels)))){
      stop("confusionMatrix has the wrong dimensions.")
    }
    if(any(colSums(confusionMatrix) == 0)){
      # We don't know the confusion structure of some label the classifier outputs.
      # Make them non-informative
      warning(paste0("confusionMatrix does not contain information on the case that ",
                     "the classifier outputs the label(s) ", 
                     colnames(confusionMatrix)[colSums(confusionMatrix) == 0],
                     ". Adding uniform distribution on those cases.", collapse = ", "))
      confusionMatrix[, colSums(confusionMatrix) == 0] <- 1 / nrow(confusionMatrix)
    }
    if(any(confusionMatrix > 1)){
      # Assume the confusionMatrix is raw
      confusionMatrix <- confusionMatrix / rep(colSums(confusionMatrix), each = nrow(confusionMatrix))
    }
    if(!all(colSums(confusionMatrix) == 1)){
      stop("Some cols of confusionMatrix do not sum to 1.")
    }
    if(!all(0 <= confusionMatrix & confusionMatrix <= 1)){
      stop("confusionMatrix contains values not in [0, 1]")
    }
    if(!is.null(rownames(confusionMatrix)) && !all(rownames(confusionMatrix) == labels)){
      stop("confusionMatrix has different rownames than labels (or different order)")
    }
    if(!is.null(colnames(confusionMatrix)) && !all(colnames(confusionMatrix) == labels)){
      stop("confusionMatrix has different colnames than labels (or different order)")
    }
  }
  colnames(confusionMatrix) <- labels
  rownames(confusionMatrix) <- labels
  
  # Find out if the labels already exist and have no classifier connected to
  # them yet:
  labelNode <- .labelSetToID(rsl, labels)
  if(!is.na(labelNode)){
    if(length(bnlearn::parents(rsl$bayesNet, labelNode)) != 0){
      stop("Labels already exist and have a classifier connected to them.")
    } else {
      # I think this is expected behaviour, so removed the warning
      # warning(paste0("labels ", labels, " already found in the rsl.", 
      #                "Connecting the classifier to those labels. Re-using old prior.", 
      #                collapse = ", "))
      prior <- rsl$labels[[labelNode]]$prior
    }
  } else {
    if(.labelsAlreadyExist(rsl, labels)){
      stop("Some of the labels already exist. Please choose new names.")
    } else {
      rsl <- addLabels(rsl, labels, prior)
      labelNode <- .labelSetToID(rsl, labels)
    }
  }
  
  # add to rsl$bayesNet
  # build cpt for the classificator node
  cID <- .getNewClassifierID(rsl)
  dimlist <- list(labels)
  names(dimlist)[1] <- cID
  # If the label prior is uniform, we can choose the classifier prior uniform, too
  # (which can avoid problems by the numerical inversion of ginv() )
  if(all(prior == 1 / length(prior))){
    cPrior <- prior
  } else {
    # Use ginv instead of solve to handle confusion matrices without full rank
    cPrior <- c(MASS::ginv(confusionMatrix) %*% prior)
    # TODO: When using ginv, make sure the result is projected into the [0,1]
    # space when the rank is not full
  }
  cProbTable <- array(cPrior, dim = length(labels), dimnames = dimlist)
  if(any(cPrior < 0 | cPrior > 1) | !all.equal(sum(cPrior), 1)){
    stop("The given confusion matrix and prior do not work together.")
  }
  
  # build cpt for the labels node
  dimlist <- list(labels, labels)
  names(dimlist) <- c(labelNode, cID)
  lProbTable <- array(confusionMatrix, dim = c(length(labels), length(labels)), dimnames = dimlist)
  
  # rebuild the rsl$bayesNet
  tables <- lapply(rsl$bayesNet, "[[", "prob")
  tables[[cID]] <- cProbTable
  tables[[labelNode]] <- lProbTable
  nodes <- c(bnlearn::nodes(rsl$bayesNet), cID)
  arcs <- rbind(bnlearn::arcs(rsl$bayesNet), cbind(cID, labelNode))
  rsl$bayesNet <- bnlearn::empty.graph(nodes = nodes)
  bnlearn::arcs(rsl$bayesNet) <- arcs
  rsl$bayesNet <- bnlearn::custom.fit(rsl$bayesNet, tables)
  rsl$needsCompilation <- TRUE
  
  # add to rsl$classifiers
  rsl$classifiers[[cID]] <- list(name = name,
                                 id = cID,
                                 confusionMatrix = confusionMatrix,
                                 prior = cPrior)
  
  # add to rsl$norn
  rsl$norn <- .addClassifier.norn(rsl$norn, cID, labelNode, confusionMatrix, cPrior)
  
  return(rsl)
}


# removeClassifier - removes a classifier from an rsl
removeClassifier <- function(rsl, name){
  # TODO: Implement (low priority)
}


# getClassifiers - returns a list with information about all classifiers
getClassifiers <- function(rsl){
  res <- lapply(rsl$classifiers, function(x){
    return(list(confusionMatrix = x$confusionMatrix,
                labels = rownames(x$confusionMatrix)))
  })
  names(res) <- sapply(rsl$classifiers, "[[", "name")
  
  return(res)
}


print.rsl <- function(rsl){
  classifiers <- getClassifiers(rsl)
  labels <- getLabels(rsl)
  rules <- getRules(rsl)
  
  cat("Rule Stacking Learner comprises ", 
      length(classifiers), " classifiers, ",
      sum(sapply(labels, length)), " labels and ",
      nrow(rules), " rules.\n\n",
      "Classifiers:\n",
      paste0(names(classifiers), rep(": ", length(classifiers)), 
             sapply(classifiers, function(x) paste0(x$labels, collapse = ", ")), 
             collapse = "\n"), "\n\n",
      "Labels:\n", 
      paste(sapply(labels, paste, collapse = " | "), collapse = "\n"), "\n\n",
      "Rules:\n",
      paste(rules$name, " (prob = ", round(rules$prob, 4), ")", sep = "", collapse = "\n"), "\n",
      sep = "")
}


plot.rsl <- function(rsl){
  # TODO: Make the plot nicer (include node names instead of IDs)
  # -> maybe with igraph?
  rsl <- .compile(rsl)
  plot(rsl$compiledNet)
}


# compile - creates a grain network for inference in the RSL
.compile <- function(rsl){
  if(rsl$needsCompilation){
    rsl$compiledNet <- bnlearn::as.grain(rsl$bayesNet)
    rsl$compiledNet <- gRain:::compile.grain(rsl$compiledNet)
    rsl$needsCompilation <- FALSE
  }
  
  return(rsl)
}


# .getAllAuxNodes - returns the IDs of all auxiliary rule nodes
.getAllAuxNodes <- function(rsl){
  return(rsl$rules$auxID)
}


# .getNewAuxID - returns an unused ID for an Aux node (e.g. "A123")
.getNewAuxID <- function(rsl){
  auxs <- rsl$rules$auxID
  auxIDs <- gsub("^A([[:digit:]]+)$", "\\1", auxs)
  maxID <- max(c(0, as.integer(auxIDs)))
  return(paste0("A", maxID + 1))
}


# .getAllRules - returns a character vector containing all rules by their ID
.getAllRules <- function(rsl){
  return(rsl$rules$ruleID)
}


# .getNewRuleID - returns an unused ID for a rule node (e.g. "R123")
.getNewRuleID <- function(rsl){
  rules <- .getAllRules(rsl)
  ruleIDs <- gsub("^R([[:digit:]]+)$", "\\1", rules)
  maxID <- max(c(0, as.integer(ruleIDs)))
  return(paste0("R", maxID + 1))
}


# .getAllLabelNodes - returns a character vector including all label group nodes
.getAllLabelNodes <- function(rsl){
  return(unname(sapply(rsl$labels, "[[", "id")))
}


# .getNewLabelID - returns an unused ID for a label node (e.g. "L123")
.getNewLabelID <- function(rsl){
  labels <- .getAllLabelNodes(rsl)
  labelIDs <- gsub("^L([[:digit:]]+)$", "\\1", labels)
  maxID <- max(c(0, as.integer(labelIDs)))
  return(paste0("L", maxID + 1))
}


# .getAllClassifiers - returns a character vector including all classifier names
.getAllClassifiers <- function(rsl){
  return(unname(sapply(rsl$classifiers, "[[", "id")))
}


# .getNewClassifierID - returns an unused ID for a classifier node (e.g. "C123")
.getNewClassifierID <- function(rsl){
  classifiers <- .getAllClassifiers(rsl)
  classifierIDs <- gsub("^C([[:digit:]]+)$", "\\1", classifiers)
  maxID <- max(c(0, as.integer(classifierIDs)))
  return(paste0("C", maxID + 1))
}


.classifierToID <- function(rsl, name){
  isClassifier <- sapply(rsl$classifiers, function(x) x$name == name)
  if(length(isClassifier) > 0 && sum(isClassifier) == 1){
    return(rsl$classifiers[[which(isClassifier)]]$id)
  } else {
    return(NA)
  }
}


# .labelsToClassiferID - for a (subset of) label(s), returns the ID of a 
#                        classifier that includes all of them (and possibly more)
.labelsToClassifierID <- function(rsl, labels){
  isClassif <- sapply(rsl$classifiers, function(x){
    classifLabels <- colnames(x$confusionMatrix)
    return(all(labels %in% classifLabels))
  })
  if(length(isClassif) > 0 && sum(isClassif) == 1){
    return(names(rsl$classifiers)[isClassif])
  } else {
    return(NA)
  }
}

# .classifierIDtoLabelID - for a given classifier ID, returns the ID of the 
#                          label set it classifies
.classifierIDtoLabelID <- function(rsl, id){
  return(bnlearn::children(rsl$bayesNet, id))
}


#.labelIDtoClassifierID - for a given label ID, returns the ID of the connected
#                         classifier (or NULL if none is connected)
.labelIDtoClassifierID <- function(rsl, id){
  return(bnlearn::parents(rsl$bayesNet, id))
}


.ruleToID <- function(rsl, name){
  # TODO: Implement (on demand)
}


# .labelSetToID - for a given vector of labels, gives back the label node id
.labelSetToID <- function(rsl, labels){
  isLabel <- sapply(rsl$labels, function(x) 
    length(x$names) == length(labels) && all(x$names == labels))
  if(length(isLabel) > 0 && sum(isLabel) == 1){
    return(rsl$labels[[which(isLabel)]]$id)
  } else {
    return(NA)
  }
}


.IDtoClassifier <- function(rsl, id){
  return(rsl$classifiers[[id]]$name)
}


.IDtoRule <- function(rsl, id){
  return(rsl$rules$name[rsl$rules$ruleID == id])
}


.IDtoLabelNode <- function(rsl, id){
  # TODO: Make the rsl actually save the label node names
  return(rsl$labels[[id]]$id)
}


.IDtoLabels <- function(rsl, id){
  return(rsl$labels[[id]]$names)
}


.IDtoClassifierLabels <- function(rsl, id){
  return(colnames(rsl$classifiers[[id]]$confusionMatrix))
}


.ruleToAux <- function(rsl, id){
  return(rsl$rules$auxID[rsl$rules$ruleID == id])
}

.auxIDtoRuleID <- function(rsl, id){
  return(rsl$rules$ruleID[rsl$rules$auxID == id])
}


# .setAuxEvidence - sets evidence of all auxiliary rule nodes to 1
.setAuxEvidence <- function(rsl, exclude = c(), propagate = FALSE){
  auxs <- .getAllAuxNodes(rsl)
  auxs <- setdiff(auxs, exclude)
  rsl$compiledNet <- gRain::setEvidence(rsl$compiledNet, nodes = auxs, 
                                        states = rep("fulfilled", length(auxs)),
                                        propagate = propagate)
  
  return(rsl)
}


# .retractAuxEvidence - removes evidence from all auxiliary nodes
.retractAuxEvidence <- function(rsl, propagate = FALSE){
  auxs <- .getAllAuxNodes(rsl)
  rsl$compiledNet <- gRain::retractEvidence(rsl$compiledNet, nodes = auxs,
                                            propagate = propagate)
  return(rsl)
}


# .makeAuxEvidence - returns a list in the format for gRain::querygrain
#                    to set all aux nodes to "fulfilled"
.makeAuxEvidence <- function(rsl, exclude = c()){
  auxs <- .getAllAuxNodes(rsl)
  auxs <- setdiff(auxs, exclude)
  auxList <- lapply(auxs, function(x) "fulfilled")
  names(auxList) <- auxs
  
  return(auxList)
}


# .setEvidence - sets (soft) evidence to all classifier inputs by editing the
#                internal compiled gRain network
# Input:
#  rsl - an rsl object
#  evidence - a list where each entry corresponds to a classifier node and has a 
#             1-row dataframe with probabilities of all labels of that node
.setEvidence <- function(rsl, evidence){
  if(length(evidence) == 0){
    return(rsl)
  }
  
  # edit classifier cpts
  rsl$compiledNet$cptlist[names(evidence)] <- lapply(evidence, function(x) as.array(as.matrix(x)))
  
  # edit clique tree
  cliques <- sapply(rsl$compiledNet$potential$pot_orig, function(x) names(dimnames(x))[1])
  for(i in seq(along = evidence)){
    cl <- which(cliques == names(evidence)[i])
    dimnames <- dimnames(rsl$compiledNet$potential$pot_orig[[cl]])
    newCPT <- t(rsl$classifiers[[names(evidence)[i]]]$confusionMatrix) * unlist(evidence[[i]])
    dimnames(newCPT) <- dimnames
    rsl$compiledNet$potential$pot_orig[[cl]] <- rsl$compiledNet$potential$pot_temp[[cl]] <-
      newCPT
  }
  
  return(rsl)
}


# .removeEvidence - resets all evidence in the grain classifier
.removeEvidence <- function(rsl){
  # Reset all prior distributions of classifiers to their original state
  priorList <- lapply(rsl$classifiers, "[[", "prior")
  rsl <- .setEvidence(rsl, priorList)
  
  return(rsl)
}


# .preprocessData - adds missing labels and missing probabilities and reorders
#                   the columns of a dataframe into the standard scheme
# Input:
#  rsl - an rsl object
#  data - a dataframe where each column is a label and contains probabilties
#  imputeMissings - if an observation contains NAs for a whole label group, should
#                   they be replaced with their prior? (note: when input is given
#                   for a couple of labels in the group, the rest will always be
#                   imputed regardless of this argument)
# Output:
#  a list where each entry corresponds to a node and contains a dataframe with
#  the imputed prior weights of the labels in that node for each observation
.preprocessData <- function(rsl, data, imputeMissings = TRUE){
  # TODO: Make sure the output is in the correct order
  # TODO: Currently removes all labels that are not connected to a classifier (low priority)
  
  # match columns of data to classifiers (https://stackoverflow.com/a/51298361)
  labelIDs <- sapply(colnames(data), .labelsToClassifierID, rsl = rsl)
  dataList <- split.default(data, labelIDs)
  relLabels <- .getAllClassifiers(rsl)
  relLabels <- relLabels[relLabels %in% labelIDs]
  dataList <- dataList[match(relLabels, names(dataList))]
  
  # check if the nodes contain all of their labels
  # Add labels and re-order if necessary
  for(i in seq(along = dataList)){
    allLabels <- .IDtoClassifierLabels(rsl, names(dataList)[i])
    order <- match(allLabels, colnames(dataList[[i]]))
    if(any(is.na(order))){
      # add missing labels
      missing <- setdiff(allLabels, colnames(dataList[[i]]))
      missingData <- matrix(NA, ncol = length(missing))
      colnames(missingData) <- missing
      dataList[[i]] <- cbind(dataList[[i]], as.data.frame(missingData))
      order <- match(allLabels, colnames(dataList[[i]]))
    }
    dataList[[i]] <- dataList[[i]][, order]
  }
  
  # Impute probabilities to missing labels and NAs
  # per observation: (1-sum(known probabilities)) / (number of missing probabilities)
  for(i in seq(along = dataList)){
    existingProb <- rowSums(dataList[[i]], na.rm = TRUE)
    nNA <- rowSums(is.na(dataList[[i]]))
    
    # If only some are missing, impute them with equal probability
    missingProb <- (1 - existingProb) / nNA
    missingData <- matrix(missingProb, nrow = nrow(dataList[[i]]), ncol = ncol(dataList[[i]]))
    
    # If all are missing, impute them with the classifier's prior
    isAllNA <- nNA == ncol(dataList[[i]])
    prior <- unlist(rsl$classifiers[[names(dataList)[i]]]$prior)
    if(imputeMissings){
      missingData[isAllNA, ] <- rep(prior, each = sum(isAllNA))
    } else {
      missingData[isAllNA, ] <- NA
    }
    
    dataList[[i]][is.na(dataList[[i]])] <- missingData[is.na(dataList[[i]])]
  }
  
  return(dataList)
}


# .predictExact.rsl - computes a-posteriori estimates of all labels
# Input:
#  rsl - an rsl object
#  data - a dataframe where each column corresponds to a label and gives the 
#         probability  of that label (not all labels have to be given, NA are allowed)
#  type - "marginal" or "joint", whether the a-posteriori estimates should be
#         marginal a-posteriori probabilities or joint MAP estimates
#         Note: joint estimates will be 0/1 encoded and may have a long runtime
#  cluster - a parallel cluster
#  showProgress - logical indicating whether progress should be printed to console
# Output:
#   a dataframe where each column gives the estimates of each label
.predictExact.rsl <- function(rsl, data, type = "marginal", cluster = NULL, showProgress = FALSE){
  # TODO: Add type checks
  # TODO: Optionally also output the rule a-posteriori probabilities
  
  if(showProgress) cat("Compiling rsl...")
  rsl <- .compile(rsl)
  
  if(showProgress) cat("Preprocessing data...")
  dataList <- .preprocessData(rsl, data)
  
  # prepare for later:
  rsl <- .setAuxEvidence(rsl)
  relevantNodes <- names(rsl$labels)
  
  # compute a-posteriori probabilities
  labelNodes <- getLabels(rsl)
  if(type == "joint"){
    labelStart <- cumsum(sapply(labelNodes, length))
    labelStart <- c(0, labelStart[-length(labelStart)])
  }
  labels <- unlist(labelNodes)
  post <- matrix(NA_real_, ncol = length(labels), nrow = nrow(data))
  colnames(post) <- labels
  rownames(post) <- rownames(data)
  post <- as.data.frame(post)
  
  clusterExport(cluster, c("rsl", "dataList", "relevantNodes", "type", "labels", "showProgress"), 
                envir = environment())
  if(showProgress) cat("Predicting...\n")
  post[] <- t(parSapply(cluster, seq(nrow(data)), function(i){
    if(showProgress) cat(i)
    observation <- lapply(dataList, "[", i, , drop = FALSE) # argument left blank on purpose
    
    # No need to remove the evidence, because in every iteration the same 
    # classificator evidences will be overridden
    #rsl <- .removeEvidence(rsl)
    rsl <- .setEvidence(rsl, observation)
    
    est <- gRain::querygrain(rsl$compiledNet, nodes = relevantNodes, type = type)
    if(type == "marginal"){
      # bring est to the order of labels
      names(est) <- NULL
      est <- unlist(est)
      # TODO: Do this matching only once in the end to save runtime
      est <- est[match(labels, names(est))]
    } else if(type == "joint"){
      maxInd <- which(est == max(est), arr.ind = TRUE)[1, ]
      # TODO: Do this only once in the end to save runtime
      maxInd <- maxInd[match(names(labelNodes), names(maxInd))]
      est <- rep(0, length(labels))
      est[labelStart + maxInd] <- 1
    }
    return(est)
  }))
  
  return(post)
}


# .predict.rsl - computes a-posteriori estimates of all labels
# Input:
#  rsl - an rsl object
#  data - a dataframe where each column corresponds to a label and gives the 
#         probability  of that label (not all labels have to be given, NA are allowed)
#  method - "exact" or "approximate" to use exact calculation or pearls linear
#           approximate algorithm (useful in big networks). Default ("auto")
#           chooses exact for small and approximate for big networks
#  cluster - a parallel cluster
#  type - "marginal" or "joint", whether the a-posteriori estimates should be
#         marginal a-posteriori probabilities or joint MAP estimates
#         Note: joint estimates will be 0/1 encoded and may have a long runtime
#  showProgress - logical indicating whether progress should be printed to console
# Output:
#   a dataframe where each column gives the estimates of each label
predict.rsl <- function(rsl, data, method = "auto", type = "marginal", cluster = NULL,
                        showProgress = FALSE){
  if(method == "auto"){
    labs <- getLabels(rsl)
    size <- sum(log(sapply(labs, length)))
    method <- ifelse(size <= 25, "exact", "approximate")
  } else if(!method %in% c("exact", "approximate")){
    stop('method must be one of "auto", "exact" or "approximate".')
  }
  
  if(method == "exact"){
    pred <- .predictExact.rsl(rsl, data, type, cluster, showProgress)
  } else if(method == "approximate"){
    pred <- predict.norn(rsl$norn, rsl, data, type, cluster, showProgress)
  }
  
  return(pred)
}


# .generateStandardOrder - creates a matrix of all possible label combinations.
#                          Each column is a label node and iterates through its
#                          possible labels using expand.grid order
.generateStandardOrder <- function(rsl){
  return(expand.grid(getLabels(rsl), stringsAsFactors = FALSE))
}


# .computeGradientHamming - computes the gradient for hamming loss of a 
#                               given combination of priors and actuals
# Input:
#  weights - the current weights of the rules in standard order
#  jointPrior - the joint prior distribution of the observation in standard order
#  actual - character vector containing the actual values for each label node
#  standardOrder - the dataframe with the standard order of labels (this could
#                  also just be generated by .generateStandardOrder(), but re-
#                  creating it each time we compute a gradient is too costy)
# Output:
#  the gradient, a numeric vector of the same length as weights
# NOTE:
#  For a better intuition: standardOrder is a matrix. jointPrior, weights and
#  the gradient are in the same order as the COLUMNS of that matrix. actual is
#  in the order of the ROWS of that matrix.
.computeJointRuleGradient <- function(weights, jointPrior, actual, standardOrder){
  # TODO: re-check if this is implemented correctly (with Lena)
  grad <- rep(0, length(weights))
  logLik <- 0
  for(l in seq(ncol(standardOrder))){
    # For each label, add its gradient to the overall log loss gradient
    labelIsCorrect <- standardOrder[, l] == actual[l]
    probLabelCorrect <- sum((jointPrior * weights)[labelIsCorrect]) / sum(jointPrior * weights)
    prefactor <- probLabelCorrect^2 / sum((weights * jointPrior)[labelIsCorrect]) * jointPrior
    gradCorrect <- prefactor * (1 - probLabelCorrect) / probLabelCorrect
    gradIncorrect <- - prefactor
    labelGrad <- ifelse(labelIsCorrect, gradCorrect, gradIncorrect)
    
    # TODO: Check if this way of avoiding NaNs in the grad is ok (high priority)
    # TODO: Change logLik etc to other names corresponding to the new method
    gradAdd <- labelGrad
    if(!any(is.nan(gradAdd))){
      # Because d/dx log(f(x)) = 1 / x * d/dx f(x)
      grad <- grad + 1 / probLabelCorrect * gradAdd
    }
    
    logLik <- logLik + log(probLabelCorrect)
  }
  
  # Note that we have to use a "* (-1)" in order to have the gradient showing
  # towards the steepest ascent (not descent), because we want to maximize 
  # likelihood, not minimize it
  grad <- -grad
  
  return(list(grad = grad, logLik = logLik))
}


# .cosSimilarity - computes the cosine similarity between two numeric vectors
.cosSimilarity <- function(a, b){
  return(sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))))
}


#.findOptJointRule - uses adam optimizer to find the weights of an optimal joint
#                    rule
# Input:
#  rsl - an rsl object
#  prior - dataset with priors of the labels (assumed to be preprocessed already)
#  actual - dataset of correct labels (assumed to be preprocessed already)
#  nRules - the desired number of rules to be learned
#  standardOrder - result of .generateStandardOrder(rsl) (handed over as 
#                  argument to avoid duplicate generation)
#  batchsize - batchsize for adam optimizer
#  alpha - hyperparameter for adam optimizer
#  beta1 - hyperparameter for adam optimizer
#  beta2 - hyperparameter for adam optimizer
#  maxIter - maximum iterations before adam optimizer is forcefully stopped
#  eps - term to avoid dividing by zero for adam optimizer
# Output:
#  a numeric vector giving the 
.findOptJointRule <- function(rsl, prior, actual, nRules, standardOrder,
                              batchsize, alpha, beta1, beta2, eps, maxIter){
  # Adam optimizer
  # The goal is to find the best "joint rule" that produces the best 
  # a-posteriori probabilities for the given inputs.
  cat("Finding best joint rule...")
  # TODO: Start with random weights or with 0.5^L?
  weights <- rep(0.5^nRules, nrow(standardOrder))
  # Adam parameters as in https://towardsdatascience.com/10-gradient-descent-optimisation-algorithms-86989510b5e9
  t <- 1 # iteration
  m <- 0 # momentum
  v <- 0 # exponential moving average of squared gradients
  repeat{
    # Compute the minibatch gradient:
    selectedObs <- sample(nrow(prior), batchsize)
    grad <- rep(0, length(weights))
    logLik <- 0
    for(obs in selectedObs){
      # It might make sense to compute the joint prior once for all observations
      # in the beginning, but it might run out of memory
      # The joint Prior gives the joint prior weight for each possible combination
      # of labels (= each row) in standardOrder
      jointPrior <- apply(standardOrder, 1, function(x)
        prod(prior[obs, colnames(prior) %in% x]))
      # This uses that labels are unique. If that is changed, we have to change
      # the computation of jointPrior here too
      ham <- .computeJointRuleGradient(weights = weights, 
                                       jointPrior = jointPrior, 
                                       actual = unlist(actual[obs, ]),
                                       standardOrder = standardOrder)
      logLik <- logLik + ham$logLik
      grad <- grad + ham$grad
    }
    
    m <- beta1 * m + (1 - beta1) * grad
    v <- beta2 * v + (1 - beta2) * grad^2
    mHat <- m / (1 - beta1^t)
    vHat <- v / (1 - beta2^t)
    # TODO: Is truncating to [0, 1] a good idea? We could also rescale to [0, 1].
    newWeights <- pmin(1, pmax(0, weights - alpha / (sqrt(vHat) + eps) * mHat))
    # TODO: Should we use another loss here because of the invariance of the BN
    #       to constant factors in the weight vector?
    diff <- sum((weights -  newWeights)^2)
    relDiff <- sqrt(diff) / sqrt(sum(weights^2))
    weights <- newWeights
    
    cat("Diff:", diff, ", relDiff:", relDiff, "logLik:", logLik, "\n")
    # temp <- list(weights = weights, diff = diff, relDiff = relDiff, logLik = logLik)
    # save(temp, file = paste0("gradDesc", t, ".RData"))
    
    t <- t + 1
    # if(diff < delta1){
    #   cat("Converged.\n")
    #   break
    # }
    if(t > maxIter){
      cat("Reached maxIter without converging.\n")
      break
    } 
  }
  
  return(weights)
}

# .splitJointRule - Splits the joint rule into local rules.
#                   Find rules such that the weights generated by that rule set 
#                   is as similar as possible to given joint rule weights
# Input:
#  rsl - an rsl object
#  weights - the weights vector of the joint rule
#  nRules - integer giving the desired number of local rules
#  standardOrder - as generated by .generateStandardOrder() (given as argument
#                  to avoid duplicated generation)
# onlyPositiveRules - logical whether all local rules should be forced to have a 
#                     probability of > 0.5
# Output:
#  a list where each entry is a list consisting of
#    ruleHead - the labels that are in the rule head
#    p - the probability of the rule
#    weights - the vector of rule weights in standard order
.splitJointRule <- function(rsl, weights, nRules, standardOrder, onlyPositiveRules){
  cat("Splitting joint rule into", nRules, "local rules...")
  # start with non-informative rules. 
  rules <- list()
  allLabels <- c(unlist(getLabels(rsl)), paste0("!", unlist(getLabels(rsl))))
  for(rule in seq(nRules)){
    # -> cos distance only uses direction, not length. So, it does not matter how
    #    big the values in the vector are as long as they are all the same. So we
    #    will just use p = 0.5
    rules[[rule]] <- list(ruleHead = c(),
                          p = 0.5,
                          weights = rep(0.5, nrow(standardOrder)))
  }
  
  # Greedy local rule extraction:
  pLower <- ifelse(onlyPositiveRules, 0.5, 0)
  rulesWeights <- rep(0.5^nRules, nrow(standardOrder))
  bestSimilarity <- .cosSimilarity(rulesWeights, weights)
  # -> why cos distance? Because it ignores the vector's length and only direction counts
  #    (and the BN always normalizes the weights - so to say - too because of 
  #    the normalization over the sum of all combinations)
  # For each rule (greedy)
  for(rule in seq(nRules)){
    cat("Rule", rule, "...")
    # Forward search through all labels for a rule head 
    best <- rules[[rule]]
    for(i in seq(length(allLabels))){
      bestChanged <- FALSE
      for(label in setdiff(allLabels, best$ruleHead)){
        proposal <- rules[[rule]]
        proposal$ruleHead <- c(proposal$ruleHead, label)
        startsWithNot <- grepl("^!", proposal$ruleHead)
        cleanedHead <- gsub("^!", "", proposal$ruleHead)
        trueLabels <- cleanedHead[!startsWithNot]
        falseLabels <- cleanedHead[startsWithNot]
        isFulfilled <- apply(standardOrder, 1, function(x){
          return(any(x %in% trueLabels) | any(!falseLabels %in% x))
        })
        # Select the p that optimizes the similarity
        p <- optimize(function(p){
          curWeights <- p * isFulfilled + (1 - p) * (!isFulfilled)
          # Replace current best proposal with our new one
          proposedWeights <- rulesWeights / rules[[rule]]$weights * curWeights
          return(.cosSimilarity(proposedWeights, weights))
        }, lower = pLower, upper = 1, maximum = TRUE)
        proposal$p <- p$maximum
        proposal$weights <- proposal$p * isFulfilled + (1 - proposal$p) * (!isFulfilled)
        
        # Check if the rule has increased the similarity
        if(p$objective > bestSimilarity){
          bestChanged <- TRUE
          best <- proposal
          bestSimilarity <- p$objective
        }
      }
      
      # Stop if similarity has not increased by adding a label
      if(bestChanged){
        rulesWeights <- rulesWeights / rules[[rule]]$weights * best$weights
        rules[[rule]] <- best
      } else {
        break
      }
    }
  }
  
  return(rules)
}


# .learnRulesJointRule - learns rules for an rsl object in order to maximize an 
#                        a-posteriori loss given data
# Input:
#  rsl - an rsl object without any rules yet (existing rules will be deleted)
#  prior - a dataframe where each column corresponds to a label and gives the 
#          probability  of that label (not all labels have to be given, NA are allowed)
#  actual - a dataframe where each column corresponds to a label node (!) and 
#           gives its correct label (all nodes have to be given, NAs are not allowed)
#           So, it has to has as many columns as there are label groups.
#  nRules - the desired number of rules to be learned
#  method - "hamming" for optimizing the a-posteriori hamming loss
#  onlyPositiveRules - logical indicating whether only rules with p in [0.5, 1] 
#                      should be searched (TRUE) or with p in [0, 1] (FALSE).
#                      The former is easier to interpret, the latter will have
#                      better prediction accuracy.
#  batchsize - batchsize for adam optimizer
#  alpha - hyperparameter for adam optimizer
#  beta1 - hyperparameter for adam optimizer
#  beta2 - hyperparameter for adam optimizer
#  maxIter1 - maximum iterations before adam optimizer is forcefully stopped
#  eps - term to avoid dividing by zero for adam optimizer
#  delta1 - convergence threshold for adam optimizer
# Output:
#  rsl object, but with added rules
.learnRulesJointRule <- function(rsl, prior, actual, nRules = 20, method = "hamming",
                                 onlyPositiveRules = FALSE, 
                                 batchsize = 20, alpha = 0.002, beta1 = 0.9, beta2 = 0.999,
                                 maxIter1 = 10000, delta1 = 1e-6, eps = 1e-8){
  # TODO: Allow rsl to have existing rules and use them as starting point
  # TODO: Allow to auto-tune the nRules by setting it to NA
  # TODO: Implement method "joint loss"
  # TODO: Make less restrictions on the actual input and impute if possible
  # TODO: check that actual is in a correct format
  # TODO: Make the gradient descent work with incomplete data
  
  # Convert the classifier priors into the actual-label priors
  # (which corresponds to predicting without any rules)
  cat("Preparing data...")
  prior <- predict(rsl, prior)
  
  # Apply Adam optimizer to find best weights for joint rule
  standardOrder <- .generateStandardOrder(rsl)
  weights <- .findOptJointRule(rsl, prior, actual, nRules, standardOrder, 
                               batchsize, alpha, beta1, beta2, eps, maxIter1)
  
  # Split the joint rule into local rules:
  rules <- .splitJointRule(rsl, weights, nRules, standardOrder, onlyPositiveRules)
  
  # Build found rules into the rsl
  cat("\nAdding local rules to rsl...\n")
  for(rule in seq(nRules)){
    # TODO: Refactor most of this into .splitJointRule()
    if(length(rules[[rule]]$ruleHead) == 0 | isTRUE(all.equal(rules[[rule]]$p, 0.5))){
      warning(paste0("Local rule ", rule, " is noninformative. Not adding this rule to rsl."))
    } else {
      startsWithNot <- grepl("^!", rules[[rule]]$ruleHead)
      cleanedHead <- gsub("^!", "", rules[[rule]]$ruleHead)
      trueLabels <- cleanedHead[!startsWithNot]
      falseLabels <- cleanedHead[startsWithNot]
      rsl <- addRule(rsl, 
                     rule = paste(paste(trueLabels, collapse = ", "), 
                                  "<-", 
                                  paste(falseLabels, collapse = ", ")), 
                     prob = rules[[rule]]$p)
    }
  }
  
  return(rsl)
}


# .logActivation - returns the logistic activation function f(x) = 1 / (1 + e^-x)
.logActivation <- function(x){
  return(1 / (1 + exp(-x)))
}


# .smoothMin - returns smoothed minimum of a number of values using logSumExp
# Input:
#  x - the values
#  alpha - a factor indicating how sharp the smooth min should approximate the
#          actual minimum (higher = sharper, closer to 0 = smoother)
.smoothMin <- function(x, alpha = 20){
  x <- -x
  m <- max(x)
  return(-(m + log(sum(exp((x - m) * alpha))) / alpha))
}


# .smoothArgMin - returns the smooth argmin (or derivation of .smoothMin)
#                 (http://erikerlandson.github.io/blog/2018/05/28/computing-smooth-max-and-its-gradients-without-over-and-underflow/)
# Input:
#  x - the values
#  alpha - a factor indicating how sharp the smooth min should approximate the
#          actual minimum (higher = sharper, closer to 0 = smoother)
.smoothArgMin <- function(x, alpha = 20){
  x <- -x
  m <- max(x)
  exped <- exp((x - m) * alpha)
  return((exped / sum(exped)))
}


# .splitProbs - splits a matrix of (inhibition) probs by the label the colnames
#               belong to
.splitProbs <- function(rsl, probs){
  # match probs to labels (https://stackoverflow.com/a/51298361)
  probs <- as.data.frame(probs) # so that it is treated as list for split.default
  
  labelIDs <- sapply(colnames(probs), .labelToID, rsl = rsl)
  probList <- split.default(probs, labelIDs)
  relLabels <- .getAllLabelNodes(rsl)
  relLabels <- relLabels[relLabels %in% labelIDs]
  probList <- probList[match(relLabels, names(probList))]
  
  return(probList)
}


# .gradSmoothMin - calculates the gradient of .smoothMin for a matrix of 
#                  inhibition probs
# Input: 
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
# Output:
#  a matrix of the same size as probs containing the gradients
.gradSmoothMin <- function(rsl, probs){
  probList <- .splitProbs(rsl, probs)
  probList <- lapply(probList, function(x) t(apply(x, 1, .smoothArgMin)))
  grad <- do.call(cbind, probList)
  
  return(grad)
}


# .gradLinRegularizer - returns the gradient of the regularizer that punishes
#                       inhibition probabilities linearly
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
# Output:
#  a matrix of the same size as probs containing the gradients
.gradLinRegularizer <- function(rsl, probs){
  return(-.gradSmoothMin(rsl, probs))
}


# .linRegularizer - returns the value of the regularizer that punishes 
#                   inhibition probs linearly
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
# Output:
#  Numeric giving the regularizer penalty
.linRegularizer <- function(rsl, probs){
  probList <- .splitProbs(rsl, probs)
  cost <- sum(1 - sapply(probList, function(x) apply(x, 1, .smoothMin)))
  
  return(cost)
}


# .gradNonLinRegularizer - returns the gradient of the regularizer that punishes
#                          inhibition probabilities non-linearly
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
# Output:
#  a matrix of the same size as probs containing the gradients
.gradNonLinRegularizer <- function(rsl, probs){
  gradF <- function(x){2 * x}
  g <- function(x){1 - x^4}
  gradG <- function(x){-4 * x^3}
  
  # prepare argmin and min
  probList <- .splitProbs(rsl, probs)
  argMin <- lapply(probList, function(x) t(apply(x, 1, .smoothArgMin)))
  min <- lapply(probList, function(x) apply(x, 1, .smoothMin))
  
  # compute g'(...) and transform into a matrix that we can multiply with argMin
  gDeriv <- list()
  for(i in seq(along = min)){
    gDeriv[[i]] <- matrix(rep(gradG(min[[i]]), ncol(argMin[[i]])), ncol = ncol(argMin[[i]]))
  }
  gDeriv <- do.call(cbind, gDeriv)
  
  # compute f'(...)
  fDeriv <- gradF(rowSums(do.call(cbind, lapply(min, g))))
  
  # grad = f'(...) * g'(...) * argMin
  grad <- fDeriv * gDeriv * do.call(cbind, argMin)
  
  return(grad)
}


# .nonLinRegularizer - returns the value of the regularizer that punishes 
#                      inhibition probs non-linearly
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
# Output:
#  Numeric giving the regularizer penalty
.nonLinRegularizer <- function(rsl, probs){
  g <- function(x){1 - x^4}
  f <- function(x){x^2}
  
  probList <- .splitProbs(rsl, probs)
  cost <- do.call(cbind, lapply(probList, function(x) g(apply(x, 1, .smoothMin))))
  cost <- sum(f(rowSums(cost)))
  
  return(cost)
}


# .gradBetaRegularizer - returns the gradient of the regularizer that punishes
#                        inhibition probabilities via a beta(alpha, beta) distribution
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
#  eps - constant to prevent dividing by 0
#  alpha - parameter for underlying beta distribution
#  beta - parameter for underlying beta distribution
# Output:
#  a matrix of the same size as probs containing the gradients
.gradBetaRegularizer <- function(rsl, probs, eps = 1e-4, alpha = 0.15, beta = 0.6){
  # prepare argmin and min
  probList <- .splitProbs(rsl, probs)
  argMin <- - do.call(cbind, lapply(probList, function(x) t(apply(x, 1, .smoothArgMin))))
  min <- lapply(probList, function(x) apply(x, 1, .smoothMin))
  min <- rep(min, times = sapply(probList, length))
  min <- 1 - do.call(cbind, min)
  
  grad <- (alpha - 1) / (min + eps) * argMin - (beta - 1) / (1 - min + eps) * argMin
  # we want to maximize the beta-likelihood, but ADAM minimizes, so take grad * (-1)
  grad <- -grad
  
  return(grad)
}


# .betaRegularizer - returns the value of the regularizer that punishes 
#                    inhibition probs via a beta(alpha, beta) distribution
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
#  eps - constant to prevent getting log(0)
#  alpha - parameter for underlying beta distribution
#  beta - parameter for underlying beta distribution
# Output:
#  Numeric giving the regularizer penalty
.betaRegularizer <- function(rsl, probs, eps = 1e-4, alpha = 0.15, beta = 0.6){
  probList <- .splitProbs(rsl, probs)
  min <- 1 - sapply(probList, function(x) apply(x, 1, .smoothMin))
  cost <- (alpha - 1) * sum(log(min + eps)) + (beta - 1) * sum(log(1 - min + eps))
  # cost is a linear transformation of the likelihood, so higher = better.
  # to make it smaller = better, take it * (-1)
  cost <- -cost
  
  return(cost)
}


# .computeNoisyORGradient - computes the gradient of hamming loss for an
#                           rsl with noisy-or rules
.computeNoisyORGradient <- function(rsl, inhProbs, actual, obs, exactness){
  # Exact method requires setting the evidence once in the beginning
  if(exactness == "exact"){
    rsl <- .setEvidence(rsl, obs)
  }
  
  # Precalculate some things we need repeatedly later
  allAuxs <- .getAllAuxNodes(rsl)
  allLabels <- .getAllLabelNodes(rsl)
  isUnknownActual <- is.na(actual)
  isAnyActualUnknown <- any(isUnknownActual)
  knownActuals <- actual[!isUnknownActual]
  if(isAnyActualUnknown){
    obsWithActuals <- obs
    for(a in seq(along = knownActuals)){
      cID <- .labelIDtoClassifierID(rsl, names(knownActuals)[a])
      obsWithActuals[[cID]][] <- as.numeric(colnames(obsWithActuals[[cID]]) == knownActuals[a])
    }
  }
  isLabelCorrect <- sapply(colnames(inhProbs), function(x) actual[.labelToID(rsl, x)] == x)
  
  # Compute gradient
  grad <- matrix(0, nrow = nrow(inhProbs), ncol = ncol(inhProbs))
  for(rule in seq(nrow(inhProbs))){
    # Condition the network on all other rules
    aID <- allAuxs[rule]
    rID <- .auxIDtoRuleID(rsl, aID)
    if(exactness == "exact"){
      rsl <- .retractAuxEvidence(rsl)
      rsl <- .setAuxEvidence(rsl, exclude = aID, propagate = TRUE)
    } else if(exactness == "approximate"){
      auxs <- allAuxs[allAuxs != aID]
      auxObs <- lapply(auxs, function(x) c(0, 1))
      names(auxObs) <- auxs
    }
    
    # compute P(rule = 1 | other rules = 1, x)
    if(exactness == "exact"){
      pRule <- gRain::querygrain(rsl$compiledNet, nodes = rID, type = "marginal")[[1]][1]
    } else if(exactness == "approximate"){
      pRule <- .beliefPropagation(rsl$norn, c(obs, auxObs), outNodes = rID)[[1]][2]
    }
    
    # Compute P(rule = 1 | all known actuals, other rules = 1, x)
    if(isAnyActualUnknown){
      if(exactness == "exact"){
        pRuleCondOnActual <- gRain::querygrain(
          gRain::setEvidence(rsl$compiledNet, nodes = names(knownActuals), states = knownActuals, propagate = TRUE), 
          nodes = rID, type = "marginal")[[1]][1]
      } else if(exactness == "approximate"){
        pRuleCondOnActual <- .beliefPropagation(rsl$norn, c(obsWithActuals, auxObs), outNodes = rID)[[1]][2]
      }
    }
    
    # Make sure we do not divide by zero or so
    if(is.nan(pRule) || pRule == 0){
      warning("Possible numeric instability while computing gradients.")
    } else {
      # compute P(each label | current rule = 0, all other rules = 1, x)
      if(exactness == "exact"){
        rsl$compiledNet <- gRain::setEvidence(rsl$compiledNet, nodes = aID, states = "not_fulfilled", propagate = TRUE)
        pLabels <- gRain::querygrain(rsl$compiledNet, nodes = allLabels, type = "marginal")
      } else if(exactness == "approximate"){
        pLabels <- .beliefPropagation(rsl$norn, c(obs, auxObs, list(aID = c(1, 0))), outNodes = allLabels)
      }
      # TODO: Make sure this is always the same order as actual
      pLabels <- pLabels[match(allLabels, names(pLabels))]
      pLabels <- unlist(pLabels)
      
      # compute P(each label | known actuals, current rule = 0, all other rules = 1)
      if(isAnyActualUnknown){
        if(exactness == "exact"){
          pLabelsCondOnActual <- gRain::querygrain(
            gRain::setEvidence(rsl$compiledNet, nodes = names(knownActuals), states = knownActuals, propagate = TRUE), 
            nodes = allLabels, type = "marginal")
          pLabelsCondOnActual[names(knownActuals)] <- obsWithActuals[!isUnknownActual]
        } else if(exactness == "approximate"){
          pLabelsCondOnActual <- .beliefPropagation(rsl$norn, c(obsWithActuals, auxObs, list(aID = c(1, 0))), outNodes = allLabels)
        }
        # TODO: Make sure this is always the same order as actual
        pLabelsCondOnActual <- pLabelsCondOnActual[match(allLabels, names(pLabelsCondOnActual))]
        pLabelsCondOnActual <- unlist(pLabelsCondOnActual)
      }
      
      for(label in seq(ncol(inhProbs))){
        pLabel <- pLabels[label]
        # compute gradient
        gr <- pLabel * (1 - pRule) / (inhProbs[rule, label] * pRule)
        if(is.na(isLabelCorrect[label])){
          pLabelCondOnActual <- pLabelsCondOnActual[label]
          gr <- gr - pLabelCondOnActual * (1 - pRuleCondOnActual) /
            (inhProbs[rule, label] * pRuleCondOnActual)
        } else if(isLabelCorrect[label]){
          if(isAnyActualUnknown){
            gr <- gr - (1 - pRuleCondOnActual) / 
              (pRuleCondOnActual * inhProbs[rule, label])
          } else {
            gr <- gr - prod(inhProbs[rule, isLabelCorrect]) /
              ((1 - prod(inhProbs[rule, isLabelCorrect])) * inhProbs[rule, label]) 
          }
        }
        grad[rule, label] <- gr
      }
    }
  }
  
  # Note that we have to use a "* (-1)" in order to have the gradient showing
  # towards the steepest ascent (not descent), because we want to maximize 
  # likelihood, not minimize it
  grad <- -grad
  
  return(grad)
}


# .truncProbs - limits the allowed number of label nodes per rule to a maximum 
#               number and sets all other label nodes' inhibition probs to 1
# Input:
#  rsl - an rsl object
#  probs - a matrix containing inhibition probs per rule (row) and label (col).
#          columns have to be named with their label.
#  maxLabels - how many label nodes are allowed per Rule
#  trunc - logical, should probs close to 1 be truncated to 1?
# Output:
#  the modified probs matrix
.truncProbs <- function(rsl, probs, maxLabels = 5, trunc = TRUE){
  # See which label nodes have the smallest inhProbs (per Rule)
  probList <- .splitProbs(rsl, probs)
  minProb <- lapply(probList, function(x) apply(x, 1, min))
  minProb <- do.call(cbind, minProb)
  ranks <- t(apply(minProb, 1, rank, ties.method = "random"))
  isSmall <- ranks <= maxLabels
  
  # Set the inhibition probs of all that are not small to 1
  for(i in seq(along = probList)){
    probList[[i]][!isSmall[, i], ] <- 1
  }
  
  newProbs <- as.matrix(do.call(cbind, probList))
  colnames(newProbs) <- colnames(probs)
  
  # Set all inhibition probs that are close to 1 to 1 to avoid numerical instability
  if(trunc){
    newProbs[newProbs > 0.99] <- 1
  }
  
  return(newProbs)
}


# . findOptNoisyOR - applies ADAM optimization to find the set inhibition 
#                    probabilities for noisy-or node that produce the best 
#                    a-posteriori probabilities for the given inputs
# output:
#  a matrix with nRules rows and each column gives a inhibition prob per label
.findOptNoisyOR <- function(rsl, prior, actual, nRules, exactness, maxIter, maxHours, 
                            batchsize, alpha, beta1, beta2, eps, initValues, reg, 
                            lambda, maxLabelsPerRule, alphaReg, betaReg, regDecay, cluster){
  # TODO: This might not work if classifiers and label nodes have different labels
  startTime <- Sys.time()
  
  # Generate start values
  labels <- getLabels(rsl)
  nLabelNodes <- length(labels)
  labels <- unlist(labels)
  if(is.na(maxLabelsPerRule) || maxLabelsPerRule > nLabelNodes){
    maxLabelsPerRule <- nLabelNodes
  }
  if(!is.null(initValues) && nrow(initValues) == nRules && ncol(initValues) == length(labels)){
    inhProbs <- initValues
  } else {
    inhProbs <- matrix(1, nrow = nRules, ncol = length(labels))
    colnames(inhProbs) <- labels
    inhProbs <- .splitProbs(rsl, inhProbs)
    # For each rule, select maxLabelsPerRule labels to have inhibition probs < 1
    for(i in seq(nRules)){
      if(maxLabelsPerRule >= length(inhProbs)){
        selectedLabels <- seq(length(inhProbs))
      } else {
        selectedLabels <- sample(length(inhProbs), maxLabelsPerRule)
      }
      for(j in selectedLabels){
        inhProbs[[j]][i, ] <- runif(ncol(inhProbs[[j]][i, ]), 0.9, 1)
      }
    }
    inhProbs <- as.matrix(do.call(cbind, inhProbs))
  }
  colnames(inhProbs) <- labels
  
  # prepare cluster
  clusterExport(cl, c("prior", "actual"), envir = environment())
  # Adam parameters as in https://towardsdatascience.com/10-gradient-descent-optimisation-algorithms-86989510b5e9
  t <- 0 # iteration
  m <- 0 # momentum
  v <- 0 # exponential moving average of squared gradients
  repeat{
    t <- t + 1
    if(t > maxIter){
      cat("Reached maxIter without converging.\n")
      break
    }
    if(difftime(Sys.time(), startTime, units = "hours") > maxHours){
      cat("Reached maxHours without converging.\n")
      break
    }
    
    # Add the new noisy-or rules to the rsl
    rsl <- .removeAllRules(rsl)
    rsl <- .addAllNoisyOR(rsl, inhProbs, 0.999)
    if(exactness == "exact"){
      rsl <- .compile(rsl)
    }
    
    
    # Compute the minibatch gradient:
    selectedObs <- sample(nrow(prior[[1]]), batchsize)
    grad <- rep(0, length(inhProbs))
    # prepare cluster
    clusterExport(cl, c("rsl", "inhProbs", "exactness"), envir = environment())
    regs <- parLapply(cluster, selectedObs, function(obs){
      # This uses that labels are unique. If that is changed, we have to change
      # the computation here too
      observation <- lapply(prior, "[", obs, , drop = FALSE) # argument left blank on purpose
      curGrad <- .computeNoisyORGradient(rsl,
                                         inhProbs = inhProbs, 
                                         actual = unlist(actual[obs, ]),
                                         obs = observation,
                                         exactness = exactness)
      
      # In case gradient explodes (due to "not possible" observations for the 
      # rules), replace with the average gradient
      isGradOk <- curGrad != Inf & curGrad != -Inf & !is.nan(curGrad)
      curGrad[!isGradOk] <- 0
      # We don't have to weight the current gradient by anything; that already
      # happens inside .computeNoisyORGradient
      
      return(curGrad)
    })
    
    for(i in seq(along = regs)){
      grad <- grad + regs[[i]]
    }
    
    # Compute the regularizer gradient
    regCost <- switch(reg,
                      "none" = 0,
                      "linear" = .linRegularizer(rsl, inhProbs),
                      "nonlinear" = .nonLinRegularizer(rsl, inhProbs),
                      "beta" = .betaRegularizer(rsl, inhProbs, alphaReg, betaReg))
    regGrad <- switch(reg,
                      "none" = 0,
                      "linear" = batchsize * nLabelNodes * lambda * .gradLinRegularizer(rsl, inhProbs),
                      "nonlinear" = batchsize * nLabelNodes * lambda * .gradNonLinRegularizer(rsl, inhProbs),
                      "beta" = batchsize * lambda * .gradBetaRegularizer(rsl, inhProbs, alphaReg, betaReg))
    
    # Put the gradient into the ADAM formula
    grad <- grad + regDecay^t * regGrad
    m <- beta1 * m + (1 - beta1) * grad
    v <- beta2 * v + (1 - beta2) * grad^2
    mHat <- m / (1 - beta1^t)
    vHat <- v / (1 - beta2^t)
    # 0.00001 to prevent dividing by 0 when computing the gradient
    newProbs <- pmin(1, pmax(0.00001, inhProbs - alpha / (sqrt(vHat) + eps) * mHat))
    inhProbs[] <- newProbs
    
    # set "almost 1" inhProbs to 1 to reduce computational workload
    inhProbs <- .truncProbs(rsl, inhProbs, maxLabels = maxLabelsPerRule,
                            trunc = FALSE)
    
    cat(as.character(Sys.time()), "\n")#, "regCost:", regCost, "\n")
    #print(inhProbs)
  }
  
  return(inhProbs)
}


# .learnRulesNoisyOR - learns rules via the noisy-or algorithm
.learnRulesNoisyOR <- function(rsl, prior, actual, nRules, exactness = "approx", 
                               maxIter = 1000, maxHours = 48,
                               batchsize = 20, alpha = 0.001, beta1 = 0.9, 
                               beta2 = 0.999, eps = 1e-8, initValues = NULL,
                               reg = "none", lambda = 1e-3, maxLabelsPerRule = Inf, 
                               alphaReg = 0.15, betaReg = 0.6, regDecay = 1,
                               cluster = cluster){
  # TODO: Allow rsl to have existing rules and use them as starting point
  # TODO: Allow to auto-tune the nRules by setting it to NA
  # TODO: Make the gradient descent work with incomplete data
  
  # Apply Adam optimizer to find best noisy OR rules
  cat("Preprocessing data...\n")
  prior <- .preprocessData(rsl, prior)
  cat("Learning...\n")
  probs <- .findOptNoisyOR(rsl, prior, actual, nRules, exactness, maxIter, maxHours,
                           batchsize, alpha, beta1, beta2, eps, initValues, reg, 
                           lambda, maxLabelsPerRule, alphaReg, betaReg, regDecay, cluster)
  
  # Add the learned rules to the rsl
  # TODO: Prune the learned rule into a "normal" rule
  rsl <- .addAllNoisyOR(rsl, probs, 0.999)
  
  return(rsl)
}


# learnRules - learns rules for an rsl object in order to maximize an 
#              a-posteriori loss given data
# Input:
#  rsl - an rsl object without any rules yet (existing rules will be deleted)
#  prior - a dataframe where each column corresponds to a label and gives the 
#          probability  of that label (not all labels have to be given, NA are allowed)
#  actual - a dataframe where each column corresponds to a label node (!) and 
#           gives its correct label (all nodes have to be given, NAs are not allowed)
#           So, it has to has as many columns as there are label groups.
#  nRules - the desired number of rules to be learned
#  method - "noisyor" to learn rules via noisy or, 
#           "jointRule" to learn via a joint rule
#  exactness - "exact" to calculate exact gradients, "approx" to approximate via
#              Pearls linear algorithm (same as in predict.rsl). "auto" to 
#              automatically decide based on size of model
#  maxIter - integer, maximum number of training iterations
#  maxHours - limit for the training time in hours
#  batchsize - batchsize for adam optimizer
#  alpha - hyperparameter for adam optimizer
#  beta1 - hyperparameter for adam optimizer
#  beta2 - hyperparameter for adam optimizer
#  eps - term to avoid dividing by zero for adam optimizer
#  initValues - matrix or vector of initial weights for the rule learners
#  reg - character giving whether to use no regularizer ("none"), the linear
#        regularizer ("linear"),the nonlinear regularizer ("nonlinear") or
#        the beta-distribution based regularizer ("beta")
#  lambda - strength of regularizer. Increase to make rules smaller. 1 / lambda
#           roughly gives the number of label nodes per rule.
#  maxLabelsPerRule - integer giving the hard limit of how many label nodes may be
#                     used per rule (as the computation cost rises exponentially).
#                     Not capped if set to Inf.
#  alphaReg - hyperparameter for beta regularizer
#  betaReg - hyperparameter for beta regularizer
#  regDecay - Factor by which regularization strength should be multiplied with 
#             after each iteration (to achieve a decay of e.g. 0.98^t). Use 1
#             to not apply a decay.
#  cluster - a cluster of the package parallel
# Output:
#  rsl object, but with added rules
learnRules <- function(rsl, prior, actual, nRules = 10, method = "noisyor", 
                       exactness = "auto", maxIter = 1000, maxHours = 48,
                       batchsize = 20, alpha = 0.001, beta1 = 0.9, 
                       beta2 = 0.999, eps = 1e-8, initValues = NULL, 
                       reg = "beta", lambda = 5, maxLabelsPerRule = 5, 
                       alphaReg = 0.15, betaReg = 0.6, regDecay = 0.98,
                       cluster = NULL){
  # TODO: Add more type checks
  if(exactness == "auto"){
    labs <- getLabels(rsl)
    size <- sum(log(sapply(labs, length)))
    exactness <- ifelse(size <= 25, "exact", "approximate")
  } else if(!exactness %in% c("exact", "approximate")){
    stop('exactness must be one of "auto", "exact" or "approximate".')
  }
  
  if(nrow(getRules(rsl)) > 0){
    stop("The rsl object must not have any rules in it already.")
  }
  if(nrow(prior) != nrow(actual)){
    stop("Unequal number of observations in prior and actual.")
  }
  # Check that actual is in a correct format
  allLabels <- .getAllLabelNodes(rsl)
  if(!is.data.frame(actual) ||
     !setequal(allLabels, colnames(actual)) ||
     !all(colnames(actual) == allLabels)){
    stop("actual has a wrong format.")
  }
  
  batchsize <- min(nrow(actual), batchsize)
  
  if(method == "jointRule"){
    return(.learnRulesJointRule(rsl, prior, actual, nRules, maxIter1 = maxIter,
                                batchsize = batchsize, alpha = alpha, beta1 = beta1,
                                beta2 = beta2, eps = eps))
  } else if(method == "noisyor"){
    return(.learnRulesNoisyOR(rsl, prior, actual, nRules, exactness = exactness,
                              maxIter = maxIter, maxHours = maxHours, batchsize = batchsize, 
                              alpha = alpha, beta1 = beta1, beta2 = beta2, eps = eps, 
                              initValues = initValues, 
                              reg = reg, lambda = lambda, 
                              maxLabelsPerRule = maxLabelsPerRule,
                              alphaReg = alphaReg, betaReg = betaReg, regDecay = regDecay,
                              cluster = cluster))
  } else {
    stop(paste0("Method ", method, " is no known learning method."))
  }
}


# .crispToProbabilisticData - turns a dataframe that columns are label nodes and
#                             is filled with the (crisp) labels of those nodes
#                             into a dataframe where each column is a label and
#                             is 1 or 0 depending on whether it is active or not
.crispToProabilisticData <- function(rsl, data){
  # TODO: Check for incomplete data and make this a public function
  # TODO: Implement Unit Tests
  datalist <- lapply(seq(ncol(data)), function(i){
    labels <- .IDtoLabels(rsl, colnames(data)[i])
    isLabel <- matrix(as.numeric(rep(data[, i], each = length(labels)) == labels), 
                      ncol = length(labels), byrow = TRUE)
    isLabel <- as.data.frame(isLabel)
    colnames(isLabel) <- labels
    return(isLabel)
  })
  
  return(do.call(cbind, datalist))
}


# .probabilisticToCrispData - turns a dataframe where each column is a label and
#                             gives its label probability into a dataframe where
#                             each column is a label node and has the label with
#                             the highest probability
# Input:
#  tieBreak - if two labels have the same probability, which should be chosen as
#             label. "random" to randomize it, "first" to always use the first
#             and "NA" to report an NA in that case.
.probabilisticToCrispData <- function(rsl, data, tieBreak = "NA"){
  # TODO: Check for incomplete data and make this a public function
  # TODO: Implement Unit Tests
  
  .whichMax <- function(probs, tieBreak){
    isMax <- which(probs == max(probs))
    if(length(isMax) == 0){
      return(NA)
    } else if(length(isMax) == 1){
      return(isMax)
    } else if(length(isMax) > 1){
      if(tieBreak == "random"){
        return(sample(isMax, 1))
      } else if(tieBreak == "first"){
        return(isMax[1])
      } else if(is.na(tieBreak) || tieBreak == "NA"){
        return(NA)
      }
    }
  }
  
  dataList <- .preprocessData(rsl, data)
  labelList <- lapply(dataList, function(x){
    colnames(x)[apply(x, 1, .whichMax, tieBreak)]
  })
  labelDataframe <- do.call(cbind, labelList)
  labelDataframe <- as.data.frame(labelDataframe, stringsAsFactors = FALSE)
  colnames(labelDataframe) <- sapply(colnames(labelDataframe), function(x) 
    .classifierIDtoLabelID(rsl, x))
  
  return(labelDataframe)
}


# hammingLoss - computes the relative amount of labels that are wrong
# Input:
#  pred - dataframe where each column is a labelset and includes its label
#         prediction per observation
#  actual - dataframe where each column is a labelset and includes its true label
#           per observation
#  na.rm - logical indicating whether NAs in pred or in actual should be ignored
hammingLoss <- function(pred, actual, na.rm = TRUE){
  if(any(sapply(pred, class) != "character") | any(sapply(actual, class) != "character")){
    stop("Wrong types. Please provide dataframes with the label names per observation.")
  }
  if(any(dim(pred) != dim(actual))){
    stop("pred and actual have different dimensions.")
  }
  if(!all(colnames(pred) %in% colnames(actual)) | !all(colnames(actual) %in% colnames(pred))){
    stop("pred and actual have different colnames.")
  }
  pred <- pred[, match(colnames(actual), colnames(pred))]
  
  return(mean(pred != actual, na.rm = na.rm))
}


# accuracy - computes the relative amount of observations in which all labels are
#            correct (keep in mind that here more is better!)
# Input:
#  pred - dataframe where each column is a labelset and includes its label
#         prediction per observation
#  actual - dataframe where each column is a labelset and includes its true label
#           per observation
#  na.rm - logical indicating whether NAs in pred or in actual should be ignored
accuracy <- function(pred, actual, na.rm = TRUE){
  if(any(sapply(pred, class) != "character") | any(sapply(actual, class) != "character")){
    stop("Wrong types. Please provide dataframes with the label names per observation.")
  }
  if(any(dim(pred) != dim(actual))){
    stop("pred and actual have different dimensions.")
  }
  if(!all(colnames(pred) %in% colnames(actual)) | !all(colnames(actual) %in% colnames(pred))){
    stop("pred and actual have different colnames.")
  }
  pred <- pred[, match(colnames(actual), colnames(pred))]
  
  return(mean(rowSums(pred == actual) == ncol(pred), na.rm = na.rm))
}


# .labelwiseLogLikelihood - given some predictions and actual labels, gives the
#                          summed labelwise log-likelihood of the true labels
# Input:
#  pred - dataframe where each column is a label and includes its predicted 
#         likelihood
#  actual - dataframe where each column is a labelset and includes its true label
#           per observation
.labelwiseLogLikelihood <- function(pred, actual){
  # TODO: Add type checks to make this a public function 
  
  logL <- numeric(nrow(pred))
  for(i in seq(along = logL)){
    logL[i] <- sum(log(pred[i, unlist(actual[i, ])]))
  }
  
  logL[logL == -Inf] <- -23
  
  return(logL)
}


# .labelwiseLikelihood - given some predictions and actual labels, gives the
#                        summed labelwise likelihood of the true labels
# Input:
#  pred - dataframe where each column is a label and includes its predicted 
#         likelihood
#  actual - dataframe where each column is a labelset and includes its true label
#           per observation
.labelwiseLikelihood <- function(pred, actual){
  # TODO: Add type checks to make this a public function 
  
  lik <- numeric(nrow(pred))
  for(i in seq(along = lik)){
    lik[i] <- sum(pred[i, unlist(actual[i, ])])
  }
  
  return(lik)
}


# .likelihood - given an rsl, priors and actual labels, calulate the likelihood
#               of the actual labels
.likelihood <- function(rsl, prior, actual, cluster){
  # TODO: implement type checks to make this a public function
  rsl <- .compile(rsl)
  
  # Preprocess prior data
  prior <- .preprocessData(rsl, prior)
  
  clusterExport(cluster, c("rsl", "prior", "actual"), 
                envir = environment())
  lik <- parSapply(cluster, seq(along = lik), function(i){
    rsl <- .retractAuxEvidence(rsl)
    observation <- lapply(prior, "[", i, , drop = FALSE) # argument left blank on purpose
    rsl <- .setEvidence(rsl, observation)
    
    # Compute likelihood of correct labels
    rsl <- .setAuxEvidence(rsl, propagate = TRUE)
    pBefore <- gRain::pEvidence(rsl$compiledNet)
    pAfter <- gRain::pEvidence(gRain::setEvidence(rsl$compiledNet, nodes = .getAllLabelNodes(rsl), 
                                                  states = unlist(actual[i, ]), propagate = TRUE))
    lik[i] <- pAfter / pBefore
  })
  
  return(lik)
}


# .avgLogLikelihood - given an rsl, priors and actual labels, calculates the 
#                     avg log Likelihood of the true labels.
#                     NOTE: give priors and actuals, do not give predictions and actuals
.avgLogLikelihood <- function(rsl, prior, actual, cluster = NULL){
  logLik <- log(.likelihood(rsl, prior, actual, cluster))
  logLik[logLik == -Inf] <- -23
  return(mean(logLik, na.rm = TRUE))
}


# .colnamesIDtoNames - replaces the IDs that might be in the colnames of a 
#                      dataframe to their corresponding names
.colnamesIDtoNames <- function(rsl, data){
  if(ncol(data) < 1){
    return(data)
  }
  
  for(i in seq(ncol(data))){
    cur <- colnames(data)[i]
    name <- c(.IDtoClassifier(rsl, cur), 
              .IDtoLabelNode(rsl, cur),
              .IDtoRule(rsl, cur))
    if(sum(!is.null(name)) == 1){
      cur <- name[!is.null(name)]
    }
    colnames(data)[i] <- cur
  }
  
  return(data)
}



# simulate - samples observations from a given rsl
# Input:
#  rsl - an rsl object
#  n - number of observations to simulate
#  outputClassifiers - boolean indicating whether the classifier inputs should 
#                      be included in the output
#  outputLabels - boolean indicating whether the (true) labels should be included 
#                 in the output
#  outputRules - boolean indicating whether the active rules should be included 
#                in the output
# Output:
#  a dataframe containing the simulated observations
simulate <- function(rsl, n, outputClassifiers = TRUE, outputLabels = TRUE, 
                     outputRules = TRUE){
  # simulate data
  rsl <- .compile(rsl)
  data <- gRain::simulate.grain(rsl$compiledNet, n)
  
  # throw out variables the user did not request
  outputVars <- character(0)
  if(outputClassifiers) outputVars <- c(outputVars, .getAllClassifiers(rsl))
  if(outputLabels) outputVars <- c(outputVars, .getAllLabelNodes(rsl))
  if(outputRules) outputVars <- c(outputVars, .getAllRules(rsl))
  data <- data[, outputVars]
  
  # make the dataframe a bit more user friendly
  data <- .colnamesIDtoNames(rsl, data)
  data[] <- lapply(data, as.character) # factors to character
  
  return(data)
}
