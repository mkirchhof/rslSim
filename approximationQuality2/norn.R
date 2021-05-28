# Implementation of Pearl's noisy-or message passing algorithm in linear time.
# Only proven to work exactly in polytrees, but appears to have some approximate
# power also in general bayesian networks.
# This version can be used in standalone, but is written to expect the inputs
# that rsl.R hands over to it and performs no further type checks, assuming that
# the inputs (given by rsl.R) are already type-checked.
# Author: michael.kirchhof@udo.edu
# Created: 07.01.2021
# version: 0.3.1 "I will always slice you"

# Dependencies:
# library(gRbase) # for tabSlice2


# create.norn - creates an empty noisyornetwork
create.norn <- function(){
  norn <- list()
  norn$classifs <- list()
  norn$labels <- list()
  norn$rules <- list()
  norn$auxs <- list()
  class(norn) <- "norn"
  
  return(norn)
}


# .addLabel.norn - adds a label node to a norn
.addLabel.norn <- function(norn, lID, prior){
  norn$labels[[lID]] <- list(prior = prior, children = character(0))
  
  return(norn)
}


# .addClassifier.norn - adds a classifier node to a norn
.addClassifier.norn <- function(norn, cID, lID, confusionMatrix, prior){
  norn$classifs[[cID]] <- list(conf = confusionMatrix,
                              prior = prior,
                              children = lID)
  norn$labels[[lID]]$parents <- cID
  
  return(norn)
}


# as.noisyornetwork - takes the labels of an rsl and casts them into a 
#                     noisyornetwork (does not copy the rules, because they are
#                     not saved in a proper structure in an rsl)
# Input:
#  rsl - an rsl object
# Output:
#  a norn (noisyornetwork) object
as.norn.rsl <- function(rsl){
  if(class(rsl) != "rsl"){
    stop("Input is not an rsl object.")
  }
  # We don't need type checks, because those all happened in the rsl already
  
  norn <- create.norn()
  
  # Add label nodes
  for(i in seq(along = rsl$labels)){
    norn <- .addLabel.norn(norn, rsl$labels[[i]]$id, rsl$labels[[i]]$prior)
  }
  
  # Add classifier nodes
  for(i in seq(along = rsl$classifiers)){
    cID <- names(rsl$classifiers)[i]
    lID <- .classifierIDtoLabelID(rsl, cID)
    norn <- .addClassifier.norn(norn, cID, lID, rsl$classifiers[[i]]$confusionMatrix,
                                rsl$classifiers[[i]]$prior)
  }
  
  return(norn)
}


# .addRule.norn - adds a noisy or rule to the network
# Input:
#  norn - a noisyornetwork
#  inhProbs - a list of numeric vectors, each containing the inhibition probs of
#             one label node to be included in the noisy-or rule. Names of the
#             list elements must be the names of the label nodes
#  rID - character, a unique rule ID (usually handed over by rsl)
#  aID - character, a unique aux ID (usually handed over by rsl)
#  p - the rule's probability
# Output:
#  the updated norn object
.addRule.norn <- function(norn, inhProbs, rID, aID, p = 1){
  # Add aux node and connect to rule node
  aCPT <- matrix(1 - p, nrow = 2, ncol = 2)
  diag(aCPT) <- p
  norn$auxs[[aID]] <- list(cpt = aCPT,
                           parents = rID)
  
  # Add new rule
  norn$rules[[rID]] <- list()
  norn$rules[[rID]]$probs <- inhProbs
  norn$rules[[rID]]$parents <- names(inhProbs)
  norn$rules[[rID]]$children <- aID
  
  # Connect to the label nodes
  for(p in norn$rules[[rID]]$parents){
    norn$labels[[p]]$children <- c(norn$labels[[p]]$children, rID)
  }
  
  return(norn)
}


# .removeRule.norn - removes a rule and its aux node from a norn
# Input:
#  norn - a norn object
#  rID - the ID of the rule to be removed
# Output:
#  the updated norn object
.removeRule.norn <- function(norn, rID){
  # Remove all pointers from label nodes to this rule
  lIDs <- norn$rules[[rID]]$parents
  for(l in lIDs){
    norn$labels[[l]]$children <- norn$labels[[l]]$children[norn$labels[[l]]$children != rID]
  }
  
  # Remove rule and aux node
  aID <- norn$rules[[rID]]$children
  norn$rules[[rID]] <- NULL
  norn$auxs[[aID]] <- NULL
  
  return(norn)
}


# .beliefPropagation - given soft or crisp input on (some) label and rule nodes, 
#                      returns posteriors of all nodes in the network using
#                      loopy belief propagation
#                 CAUTION: Does not do any error checking and expects complete input
#                 CAUTION: May not converge and produce bad approximations
# Inputs:
#  norn - a noisyornetwork object
#  input - a list where priors on classifiers and auxs can be specified. 
#          entries are vectors giving priors for classifiers or c(0, 1) for aux.
#          NOTE: for auxs, expects the order "not_fulfilled", "fulfilled", for
#                classifiers their individual correct order
#  maxit - maximum number of iterations to run the loopy belief propagation
#  convThresh - if the beliefs of all variables change by at most this much 
#               the algorithm is said to have converged
#  outnodes - a character vector containing the node ids of the nodes whose 
#             marginals shall be outputted (NULL means all nodes)
# Output:
#  a list just like input, but with the a-posteriori estimates
.beliefPropagation <- function(norn, input, maxit = 20, convThresh = 1e-5, 
                               outNodes = NULL){
  classifs <- names(norn$classifs)
  labels <- names(norn$labels)
  rules <- names(norn$rules)
  auxs <- names(norn$auxs)
  if(is.null(outNodes)){
    outNodes <- c(classifs, labels, rules, auxs)
  }
  
  # Initialize messages
  ownPi <- list()
  ownLambda <- list()
  ownBel <- list()
  piToFrom <- list()
  lambdaToFrom <- list()
  # Initialize all pi messages from label to aux
  for(c in classifs){
    if(!is.null(input[[c]])){
      ownPi[[c]] <- unlist(input[[c]])
    } else {
      ownPi[[c]] <- norn$classifs[[c]]$prior
    }
  }
  for(l in labels){
    piToFrom[[l]] <- list()
    c <- norn$labels[[l]]$parents
    piToFrom[[l]][[c]] <- c(norn$classifs[[c]]$conf %*% ownPi[[c]])
  }
  for(l in labels){
    ownPi[[l]] <- rep(1, length(norn$labels[[l]]$prior))
    for(c in norn$labels[[l]]$parents){
      ownPi[[l]] <- ownPi[[l]] * piToFrom[[l]][[c]]
    }
  }
  for(r in rules){
    piToFrom[[r]] <- list()
    for(l in norn$rules[[r]]$parents){
      piToFrom[[r]][[l]] <- ownPi[[l]]
    }
  }
  for(r in c(auxs, rules)){
    ownPi[[r]] <- rep(1, 2)
  }
  for(a in auxs){
    piToFrom[[a]] <- list()
    r <- norn$auxs[[a]]$parents
    piToFrom[[a]][[r]] <- c(norn$auxs[[a]]$cpt %*% ownPi[[r]])
  }
  for(a in auxs){
    r <- norn$auxs[[a]]$parents
    ownPi[[a]] <- piToFrom[[a]][[r]]
  }
  
  # Now, initialize the lambda messages from auxs to labels
  for(a in auxs){
    if(!is.null(input[[a]])){
      ownLambda[[a]] <- unlist(input[[a]])
    } else {
      ownLambda[[a]] <- c(1, 1)
    }
  }
  for(r in rules){
    lambdaToFrom[[r]] <- list()
    a <- norn$rules[[r]]$children
    lambdaToFrom[[r]][[a]] <- c(norn$auxs[[a]]$cpt %*% ownLambda[[a]])
  }
  for(r in rules){
    a <- norn$rules[[r]]$children
    ownLambda[[r]] <- lambdaToFrom[[r]][[a]]
  }
  for(l in c(labels, classifs)){
    ownLambda[[l]] <- rep(1, length(ownPi[[l]]))
  }
  for(l in labels){
    lambdaToFrom[[l]] <- list()
    for(r in norn$labels[[l]]$children){
      lambdaToFrom[[l]][[r]] <- rep(1, length(ownPi[[l]]))
    }
  }
  for(c in classifs){
    lambdaToFrom[[c]] <- list()
    l <- norn$classifs[[c]]$children
    lambdaToFrom[[c]][[l]] <- c(norn$classifs[[c]]$conf %*% ownLambda[[l]])
  }
  
  # Compute initial own beliefs
  for(n in c(labels)){
    ownBel[[n]] <- ownLambda[[n]] * ownPi[[n]]
    ownBel[[n]] <- ownBel[[n]] / sum(ownBel[[n]])
  }
  
  # Do the loopy belief propagation
  converged <- FALSE
  for(i in seq(maxit)){
    # Do not need to send pi messages from classifiers to labels, because remains
    # unchanged
    
    # Do not need to compute ownPi of labels, because remains unchanged
    
    # Send pi messages from label nodes to rule nodes
    for(l in labels){
      for(r in norn$labels[[l]]$children){
        piToFrom[[r]][[l]] <- ownPi[[l]]
        for(r2 in norn$labels[[l]]$children){
          if(r2 != r){
            piToFrom[[r]][[l]] <- piToFrom[[r]][[l]] * lambdaToFrom[[l]][[r2]]
          }
        }
        piToFrom[[r]][[l]] <- piToFrom[[r]][[l]] / sum(piToFrom[[r]][[l]])
      }
    }
    
    # compute helping variables for the rule nodes
    prodParts <- matrix(1, nrow = length(rules), ncol = length(labels))
    colnames(prodParts) <- labels
    rownames(prodParts) <- rules
    for(r in rules){
      for(l in norn$rules[[r]]$parents){
        prodParts[r, l] <- sum(norn$rules[[r]]$probs[[l]] * piToFrom[[r]][[l]])
      }
    }
    
    # compute ownPi beliefs of rule nodes
    for(r in rules){
      prod <- prod(prodParts[r, ])
      ownPi[[r]] <- c(prod, 1 - prod)
    }
    
    # send pi messages from rules to auxs (which are equal to the auxs ownPi, 
    # because it only has one parent)
    for(r in rules){
      a <- norn$rules[[r]]$children
      ownPi[[a]] <- c(norn$auxs[[a]]$cpt %*% ownPi[[r]])
    }
    
    # No need to compute ownlambda of auxs or to send lambda from aux to rules
    # because evidence on rules does not change
    
    # send lambda messages from rule nodes to label nodes
    for(r in rules){
      for(l in norn$rules[[r]]$parents){
        prod <- prod(prodParts[r, norn$rules[[r]]$parents[norn$rules[[r]]$parents != l]])
        lambdaToFrom[[l]][[r]] <- ownLambda[[r]][1] * norn$rules[[r]]$probs[[l]] * prod +
          ownLambda[[r]][2] * (1 - norn$rules[[r]]$probs[[l]] * prod)
        lambdaToFrom[[l]][[r]] <- lambdaToFrom[[l]][[r]] / sum(lambdaToFrom[[l]][[r]])
      }
    }
    
    # compute ownLambda beliefs of label nodes
    for(l in labels){
      ownLambda[[l]] <- rep(1, length(ownLambda[[l]]))
      for(r in norn$labels[[l]]$children){
        ownLambda[[l]] <- ownLambda[[l]] * lambdaToFrom[[l]][[r]]
      }
    }
    
    # send lambda messages from label nodes to classif nodes
    for(c in classifs){
      l <- norn$classifs[[c]]$children
      lambdaToFrom[[c]][[l]] <- c(t(ownLambda[[l]]) %*% norn$classifs[[c]]$conf)
    }
    
    # compute ownLambda of classifier nodes
    for(c in classifs){
      l <- norn$classifs[[c]]$children
      ownLambda[[c]] <- lambdaToFrom[[c]][[l]]
    }
    
    # compute ownBel beliefs of all nodes
    newOwnBel <- ownBel
    for(n in c(labels)){
      newOwnBel[[n]] <- ownPi[[n]] * ownLambda[[n]]
      newOwnBel[[n]] <- newOwnBel[[n]] / sum(newOwnBel[[n]])
    }
    
    # Check if beliefs are similar to previous iteration (then it has converged)
    converged <- TRUE
    for(i in seq(along = ownBel)){
      if(any(abs(ownBel[[i]] - newOwnBel[[i]]) > convThresh)){
        converged <- FALSE
        break
      }
    }
    
    ownBel <- newOwnBel
    if(converged) break
  }
  
  # Compute own beliefs
  for(n in c(classifs, labels, rules, auxs)){
    ownBel[[n]] <- ownLambda[[n]] * ownPi[[n]]
    ownBel[[n]] <- ownBel[[n]] / sum(ownBel[[n]])
  }
  
  if(!converged) warning("Loopy belief propagation has not converged.")
  return(ownBel[outNodes])
}


# .beliefPropagationJoint - given soft or crisp input on (some) label and rule nodes, 
#                           returns MAP interpretation in the network using
#                           loopy belief propagation
#                 CAUTION: Does not do any error checking and expects complete input
#                 CAUTION: May not converge and produce bad approximations
#                 CAUTION: This function can produce NaNs and errors when you feed
#                          it probabilities of 0 and 1
# Inputs:
#  norn - a noisyornetwork object
#  input - a list where priors on classifiers and auxs can be specified. 
#          entries are vectors giving priors for classifiers or c(0, 1) for aux.
#          NOTE: for auxs, expects the order "not_fulfilled", "fulfilled", for
#                classifiers their individual correct order
#  maxit - maximum number of iterations to run the loopy belief propagation
#  convThresh - if the beliefs of all variables change by at most this much 
#               the algorithm is said to have converged
#  outnodes - a character vector containing the node ids of the nodes whose 
#             marginals shall be outputted (NULL means all nodes)
# Output:
#  a list just like input, but with the a-posteriori estimates
.beliefPropagationJoint <- function(norn, input, maxit = 20, convThresh = 1e-5, 
                                    outNodes = NULL){
  classifs <- names(norn$classifs)
  labels <- names(norn$labels)
  rules <- names(norn$rules)
  auxs <- names(norn$auxs)
  if(is.null(outNodes)){
    outNodes <- c(classifs, labels, rules, auxs)
  }
  
  # Initialize messages
  ownPi <- list()
  ownLambda <- list()
  ownBel <- list()
  piToFrom <- list()
  lambdaToFrom <- list()
  # Initialize all pi messages from label to aux
  for(c in classifs){
    if(!is.null(input[[c]])){
      ownPi[[c]] <- unlist(input[[c]])
    } else {
      ownPi[[c]] <- norn$classifs[[c]]$prior
    }
  }
  for(l in labels){
    piToFrom[[l]] <- list()
    c <- norn$labels[[l]]$parents
    piToFrom[[l]][[c]] <- apply(norn$classifs[[c]]$conf, 1, function(x) max(x * ownPi[[c]]))
  }
  for(l in labels){
    ownPi[[l]] <- rep(1, length(norn$labels[[l]]$prior))
    for(c in norn$labels[[l]]$parents){
      ownPi[[l]] <- ownPi[[l]] * piToFrom[[l]][[c]]
    }
  }
  for(r in rules){
    piToFrom[[r]] <- list()
    for(l in norn$rules[[r]]$parents){
      piToFrom[[r]][[l]] <- ownPi[[l]]
    }
  }
  for(r in c(auxs, rules)){
    ownPi[[r]] <- rep(1, 2)
  }
  for(a in auxs){
    piToFrom[[a]] <- list()
    r <- norn$auxs[[a]]$parents
    piToFrom[[a]][[r]] <- apply(norn$auxs[[a]]$cpt, 1, function(x) max(x * ownPi[[r]]))
  }
  for(a in auxs){
    r <- norn$auxs[[a]]$parents
    ownPi[[a]] <- piToFrom[[a]][[r]]
  }
  
  # Now, initialize the lambda messages from auxs to labels
  for(a in auxs){
    if(!is.null(input[[a]])){
      ownLambda[[a]] <- unlist(input[[a]])
    } else {
      ownLambda[[a]] <- c(1, 1)
    }
  }
  for(r in rules){
    lambdaToFrom[[r]] <- list()
    a <- norn$rules[[r]]$children
    lambdaToFrom[[r]][[a]] <- apply(norn$auxs[[a]]$cpt, 1, function(x) max(x * ownLambda[[a]]))
  }
  for(r in rules){
    a <- norn$rules[[r]]$children
    ownLambda[[r]] <- lambdaToFrom[[r]][[a]]
  }
  for(l in c(labels, classifs)){
    ownLambda[[l]] <- rep(1, length(ownPi[[l]]))
  }
  for(l in labels){
    lambdaToFrom[[l]] <- list()
    for(r in norn$labels[[l]]$children){
      lambdaToFrom[[l]][[r]] <- rep(1, length(ownPi[[l]]))
    }
  }
  for(c in classifs){
    lambdaToFrom[[c]] <- list()
    l <- norn$classifs[[c]]$children
    lambdaToFrom[[c]][[l]] <- apply(norn$classifs[[c]]$conf, 1, function(x) max(x * ownLambda[[l]]))
  }
  
  # Compute initial own beliefs
  for(n in c(labels)){
    ownBel[[n]] <- ownLambda[[n]] * ownPi[[n]]
    ownBel[[n]] <- ownBel[[n]] / sum(ownBel[[n]])
  }
  
  # Do the loopy belief propagation
  converged <- FALSE
  for(i in seq(maxit)){
    # Do not need to send pi messages from classifiers to labels, because remains
    # unchanged
    
    # Do not need to compute ownPi of labels, because remains unchanged
    
    # Send pi messages from label nodes to rule nodes
    for(l in labels){
      for(r in norn$labels[[l]]$children){
        piToFrom[[r]][[l]] <- ownPi[[l]]
        for(r2 in norn$labels[[l]]$children){
          if(r2 != r){
            piToFrom[[r]][[l]] <- piToFrom[[r]][[l]] * lambdaToFrom[[l]][[r2]]
          }
        }
        piToFrom[[r]][[l]] <- piToFrom[[r]][[l]] / sum(piToFrom[[r]][[l]])
      }
    }
    
    # Compute helping values for the rule node
    f0 <- list()
    f1 <- list()
    for(r in rules){
      if(!is.null(norn$rules[[r]]$parents)){
        piProd <- array(apply(expand.grid(piToFrom[[r]]), 1, prod), dim = sapply(piToFrom[[r]], length))
        triggerProd <- array(apply(expand.grid(norn$rules[[r]]$probs), 1, prod), dim = sapply(norn$rules[[r]]$probs, length))
        f0[[r]] <- piProd * triggerProd
        f1[[r]] <- piProd - f0[[r]]
      } else {
        f0[[r]] <- 1
        f1[[r]] <- 1
      }
    }
    
    # compute ownPi beliefs of rule nodes
    for(r in rules){
      ownPi[[r]] <- c(max(f0[[r]]), max(f1[[r]]))
    }
    
    # send pi messages from rules to auxs (which are equal to the auxs ownPi, 
    # because it only has one parent)
    for(r in rules){
      a <- norn$rules[[r]]$children
      ownPi[[a]] <- apply(norn$auxs[[a]]$cpt, 1, function(x) max(x * ownPi[[r]]))
    }
    
    # No need to compute ownlambda of auxs or to send lambda from aux to rules
    # because evidence on rules does not change
    
    # send lambda messages from rule nodes to label nodes
    for(r in rules){
      for(lInd in seq(along = norn$rules[[r]]$parents)){
        l <- norn$rules[[r]]$parents[lInd]
        lambdaToFrom[[l]][[r]] <- 
          pmax(sapply(seq(along = norn$rules[[r]]$probs[[l]]), function(x) 
            max(gRbase::tabSlice2(f0[[r]], slice = list(x), margin.idx = lInd))) * ownLambda[[r]][1], 
            sapply(seq(along = norn$rules[[r]]$probs[[l]]), function(x) 
              max(gRbase::tabSlice2(f1[[r]], slice = list(x), margin.idx = lInd))) * ownLambda[[r]][2])
        # Prevent division by 0 (it is ok to "ignore" this division, because if
        # piToFrom was 0, then ownPi would be 0, thus when calculating ownBel 
        # it would stay 0 regardless)
        is0 <- piToFrom[[r]][[l]] == 0
        lambdaToFrom[[l]][[r]][!is0] <- lambdaToFrom[[l]][[r]][!is0] / piToFrom[[r]][[l]][!is0]
        
        normConst <- sum(lambdaToFrom[[l]][[r]])
        if(normConst != 0){
          lambdaToFrom[[l]][[r]] <- lambdaToFrom[[l]][[r]] / normConst
        }
      }
    }
    
    # compute ownLambda beliefs of label nodes
    for(l in labels){
      ownLambda[[l]] <- rep(1, length(ownLambda[[l]]))
      for(r in norn$labels[[l]]$children){
        ownLambda[[l]] <- ownLambda[[l]] * lambdaToFrom[[l]][[r]]
      }
    }
    
    # send lambda messages from label nodes to classif nodes
    for(c in classifs){
      l <- norn$classifs[[c]]$children
      lambdaToFrom[[c]][[l]] <- apply(norn$classifs[[c]]$conf, 1, function(x) max(x * ownLambda[[l]]))
    }
    
    # compute ownLambda of classifier nodes
    for(c in classifs){
      l <- norn$classifs[[c]]$children
      ownLambda[[c]] <- lambdaToFrom[[c]][[l]]
    }
    
    # compute ownBel beliefs of all nodes
    newOwnBel <- ownBel
    for(n in c(labels)){
      newOwnBel[[n]] <- ownPi[[n]] * ownLambda[[n]]
      newOwnBel[[n]] <- newOwnBel[[n]] / sum(newOwnBel[[n]])
    }
    
    # Check if beliefs are similar to previous iteration (then it has converged)
    converged <- TRUE
    for(i in seq(along = ownBel)){
      if(any(abs(ownBel[[i]] - newOwnBel[[i]]) > convThresh)){
        converged <- FALSE
        break
      }
    }
    
    ownBel <- newOwnBel
    if(converged) break
  }
  
  # Compute own beliefs
  for(n in c(classifs, labels, rules, auxs)){
    ownBel[[n]] <- ownLambda[[n]] * ownPi[[n]]
    # cast to crisp values
    maxInd <- which.max(ownBel[[n]])
    ownBel[[n]][maxInd] <- 1
    ownBel[[n]][-maxInd] <- 0
  }
  
  if(!converged) warning("Loopy belief propagation has not converged.")
  return(ownBel[outNodes])
}


# predict.norn - approximate predictions with the same API as predict.rsl
predict.norn <- function(norn, rsl, data, type = "marginal", cluster = NULL, showProgress = FALSE){
  if(showProgress) cat("Preprocessing data...")
  dataList <- .preprocessData(rsl, data)
  
  auxs <- .getAllAuxNodes(rsl)
  auxObs <- lapply(auxs, function(x) c(0, 1))
  names(auxObs) <- auxs
  
  # compute a-posteriori probabilities
  labelNodes <- getLabels(rsl)
  labels <- unlist(labelNodes)
  post <- matrix(NA_real_, ncol = length(labels), nrow = nrow(data))
  colnames(post) <- labels
  rownames(post) <- rownames(data)
  post <- as.data.frame(post)
  
  clusterExport(cluster, c("norn", "dataList", "auxObs", "type", "labelNodes", "showProgress"), 
                envir = environment())
  if(showProgress) cat("Predicting...\n")
  post[] <- t(parSapply(cluster, seq(nrow(data)), function(i){
    if(showProgress) cat(i)
    observation <- lapply(dataList, "[", i, , drop = FALSE) # argument left blank on purpose
    observation <- c(observation, auxObs)
    
    if(type == "marginal"){
      est <- .beliefPropagation(norn, observation, outNodes = names(labelNodes))
    } else if(type == "joint"){
      est <- .beliefPropagationJoint(norn, observation, outNodes = names(labelNodes))
    }
    est <- unlist(est)
    
    return(est)
  }))
  
  return(post)
}
