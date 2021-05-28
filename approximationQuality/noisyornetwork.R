# Implementation of Pearl's noisy-or message passing algorithm in linear time.
# Only proven to work exactly in polytrees, but appears to have some approximate
# power also in general bayesian networks.
# Author: michael.kirchhof@udo.edu
# Created: 14.12.2020
# version: 0.1.0 "Makes me feel good"

# as.noisyornetwork - takes the labels of an rsl and casts them into a 
#                     noisyornetwork (does not copy the rules, because they are
#                     not saved in a proper structure in an rsl)
#                     Note: We don't provide direct constructors for noisyornet,
#                           because we want to have the rsl code take care of
#                           various error checking etc
# Input:
#  rsl - an rsl object
# Output:
#  a norn (noisyornetwork) object
as.norn.rsl <- function(rsl){
  norn <- list()
  class(norn) <- "norn"
  
  # Add label nodes
  norn$labels <- rsl$labels
  for(i in seq(along = norn$labels)){
    norn$labels[[i]]$children <- character(0)
  }
  
  # Add (empty) rule nodes
  norn$rules <- list()
  
  return(norn)
}


# .getNewRuleID.norn - returns an unused ID for a rule of the form "R1", "R2", ...
.getNewRuleID.norn <- function(norn){
  ruleIDs <- gsub("^R([[:digit:]]+)$", "\\1", names(norn$rules))
  maxID <- max(c(0, as.integer(ruleIDs)))
  return(paste0("R", maxID + 1))
}


# addRule.norn - adds a noisy or rule to the network
# Input:
#  norn - a noisyornetwork
#  inhProbs - a list of numeric vectors, each containing the inhibition probs of
#             one label node to be included in the noisy-or rule. Names of the
#             list elements must be the names of the label nodes
# Output:
#  the updated norn object
addRule.norn <- function(norn, inhProbs){
  rID <- .getNewRuleID.norn(norn)
  
  # Add new rule
  norn$rules[[rID]] <- list()
  norn$rules[[rID]]$probs <- inhProbs
  norn$rules[[rID]]$parents <- names(inhProbs)
  
  # Connect to the label nodes
  for(p in norn$rules[[rID]]$parents){
    norn$labels[[p]]$children <- c(norn$labels[[p]]$children, rID)
  }
  
  return(norn)
}


# .predict.norn - given soft or crisp input on (some) label and rule nodes, 
#                 returns posteriors of all nodes in the network using
#                 loopy belief propagation
#                 CAUTION: Does not do any error checking and expects complete input
#                 CAUTION: May not converge and produce bad approximations
# Inputs:
#  norn - a noisyornetwork object
#  input - a list where each element is named after a node in the network and
#          contains a numeric vector for the initial beliefs of that node.
#          For rule node, it is assumed that they come in the order
#          c(not_fulfilled, fulfilled).
#  maxit - maximum number of iterations to run the loopy belief propagation
#  convThresh - if the beliefs of all variables change by at most this much 
#               the algorithm is said to have converged
# Output:
#  a list just like input, but with the a-posteriori estimates
.predict.norn <- function(norn, input, maxit = 20, convThresh = 1e-4){
  rules <- names(norn$rules)
  labels <- names(norn$labels)
  
  # Initialize messages
  ownPi <- list()
  for(r in rules){
    ownPi[[r]] <- rep(1, 2)
  }
  for(l in labels){
    ownPi[[l]] <- input[[l]]
  }
  
  ownLambda <- list()
  for(r in rules){
    ownLambda[[r]] <- input[[r]]
  }
  for(l in labels){
    ownLambda[[l]] <- rep(1, length(ownPi[[l]]))
  }
  
  ownBel <- list()
  for(n in c(labels, rules)){
    ownBel[[n]] <- ownLambda[[n]] * ownPi[[n]]
    ownBel[[n]] <- ownBel[[n]] / sum(ownBel[[n]])
  }
  
  piToFrom <- list()
  for(r in rules){
    piToFrom[[r]] <- list()
    for(l in norn$rules[[r]]$parents){
      piToFrom[[r]][[l]] <- ownPi[[l]]
    }
  }
  
  lambdaToFrom <- list()
  for(l in labels){
    lambdaToFrom[[l]] <- list()
    for(r in norn$labels[[l]]$children){
      lambdaToFrom[[l]][[r]] <- rep(1, length(ownPi[[l]]))
    }
  }
  
  # Do the loopy belief propagation
  converged <- FALSE
  for(i in seq(maxit)){
    # Send pi messages from label nodes to rule nodes
    for(l in labels){
      for(r in norn$labels[[l]]$children){
        piToFrom[[r]][[l]] <- ownPi[[l]] * Reduce("*", c(lambdaToFrom[[l]][setdiff(norn$labels[[l]]$children, r)], list(rep(1, length(ownPi[[l]])))))
        piToFrom[[r]][[l]] <- piToFrom[[r]][[l]] / sum(piToFrom[[r]][[l]])
      }
    }
    
    # send lambda messages from rule nodes to label nodes
    prodParts <- matrix(1, nrow = length(rules), ncol = length(labels))
    colnames(prodParts) <- labels
    rownames(prodParts) <- rules
    for(r in rules){
      for(l in norn$rules[[r]]$parents){
        prodParts[r, l] <- sum(norn$rules[[r]]$probs[[l]] * piToFrom[[r]][[l]])
      }
    }
    for(r in rules){
      for(l in norn$rules[[r]]$parents){
        prod <- prod(prodParts[r, setdiff(norn$rules[[r]]$parents, l)])
        lambdaToFrom[[l]][[r]] <- ownLambda[[r]][1] * norn$rules[[r]]$probs[[l]] * prod +
          ownLambda[[r]][2] * (1 - norn$rules[[r]]$probs[[l]] * prod)
        lambdaToFrom[[l]][[r]] <- lambdaToFrom[[l]][[r]] / sum(lambdaToFrom[[l]][[r]])
      }
    }
    
    # compute ownPi beliefs
    for(r in rules){
      prod <- prod(prodParts[r, ])
      ownPi[[r]] <- c(prod, 1 - prod)
    }
    
    # compute ownLambda beliefs
    for(l in labels){
      ownLambda[[l]] <- Reduce("*", c(lambdaToFrom[[l]], list(rep(1, length(ownLambda[[l]])))))
    }
    
    # compute ownBel beliefs
    newOwnBel <- ownBel
    for(n in c(labels, rules)){
      newOwnBel[[n]] <- ownPi[[n]] * ownLambda[[n]]
      newOwnBel[[n]] <- newOwnBel[[n]] / sum(newOwnBel[[n]])
    }
    
    # Check if beliefs are similar to previous iteration (then it has converged)
    # TODO
    converged <- isTRUE(all.equal(ownBel, newOwnBel, tolerance = convThresh))
    ownBel <- newOwnBel
    if(converged) break
  }
  
  if(!converged) warning("Loopy belief propagation has not converged.")
  return(ownBel)
}


# predict.norn - approximate predictions with the same API as predict.rsl
predict.norn <- function(norn, rsl, data, showProgress = FALSE){
  if(showProgress) cat("Preprocessing data...")
  dataList <- .preprocessData(rsl, data)
  for(i in seq(along = dataList)){
    names(dataList)[i] <- .classifierIDtoLabelID(rsl, names(dataList)[i])
  }
  
  rules <- .getAllRules(rsl)
  ruleObs <- lapply(rules, function(x) c(0, 1))
  names(ruleObs) <- rules
  
  # compute a-posteriori probabilities
  labelNodes <- getLabels(rsl)
  labels <- unlist(labelNodes)
  post <- matrix(NA_real_, ncol = length(labels), nrow = nrow(data))
  colnames(post) <- labels
  rownames(post) <- rownames(data)
  post <- as.data.frame(post)
  if(showProgress) cat("Predicting...\n")
  for(i in seq(nrow(data))){
    if(showProgress) cat(i, "/", nrow(data), "\n")
    observation <- lapply(dataList, "[", i, , drop = FALSE) # argument left blank on purpose
    observation <- c(observation, ruleObs)
    
    est <- .predict.norn(norn, observation)
    est[rules] <- NULL
    est <- unlist(est)
    
    post[i, ] <- est
  }
  
  return(post)
}
