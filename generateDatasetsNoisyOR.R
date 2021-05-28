# This script contains all functions necessary to create simulation datasets
# Author: michael.kirchhof@udo.edu
# Date: 07.10.2020


# Dependencies:
source("rsl.R")
# library(MCMCpack)


# .generateLabelSet - builds a list of nLabels labels containing a random number
#                     of labels (at most maxCategories) each
.generateLabelSet <- function(nLabels, maxCategories){
  if(maxCategories < 2){
    stop("maxCategories must be >= 2")
  }
  if(nLabels < 1){
    stop("nLabels must be >= 1")
  }
  
  nCategories <- replicate(nLabels, sample(2:maxCategories, 1))
  
  labelSet <- list()
  for(i in seq(nLabels)){
    labelName <- paste0("L", i)
    categories <- paste0(labelName, "_", seq(nCategories[i]))
    labelSet[[labelName]] <- categories
  }
  
  return(labelSet)
}


# .rtriangle - generates n random values from a triangle distribution
.rtriangle <- function(n, min = 0, max = 1){
  return(min + sqrt(runif(n)) * (max - min))
}


# .generateRuleSet - given a list of labels, generates random rules in the form
#                    "B, C <- A, D"
.generateRuleSet <- function(labels, nRules, maxLabelsPerRule){
  .generateRule <- function(labels, maxLabelsPerRule){
    nLabels <- sample(2:maxLabelsPerRule, 1)
    
    # We want to sample such that each category (the things inside the label groups)
    # has the same chance of being in the rule, so we sample on categories
    categories <- sample(unlist(labels), nLabels)
    
    # Sample a rule probability from a triangle distribution
    p <- .rtriangle(length(categories), min = 1, max = 0)
    
    return(list(labels = categories, p = p))
  }
  
  ruleSet <- lapply(seq(nRules), function(x)
    .generateRule(labels, maxLabelsPerRule))
  
  return(ruleSet)
}


# .generateDirichletPriors - generate prior probabilities for the labels via
#                            dirichlet distributions
.generateDirichletPriors <- function(rsl, n, alpha){
  # .genPriorForLabel - generate dirichlet for a single label group
  .genPriorForLabel <- function(labels){
    nLabels <- length(labels)
    rndDir <- MCMCpack::rdirichlet(n, c(alpha, rep(1, nLabels - 1)))
    
    # right now rndDir is always biased towards the first entry. 
    # Sample it so that each label has high probability sometimes
    biasTowardsLabel <- sample(nLabels, n, replace = TRUE)
    for(i in seq(nLabels)){
      # order will be for e.g. for i = 1: (1, 2, 3), for i = 2: (2, 3, 1) etc.
      order <- seq(i,  i + nLabels - 1) %% nLabels
      order[order == 0] <- nLabels
      rndDir[i == biasTowardsLabel, order] <- rndDir[i == biasTowardsLabel, ]
    }
    
    colnames(rndDir) <- labels
    
    return(rndDir)
  }
  
  probs <- do.call(cbind, lapply(getLabels(rsl), .genPriorForLabel))
  probs <- as.data.frame(probs)
  
  return(probs)
}


# evalPerformance - evaluates hamming loss, accuracy and log-likelihood on 
#                   test, validation and train data
.evalPerformance <- function(rsl, train, trainActual, val, valActual, test, testActual){
  cat("Predicting on train...\n")
  predTrainRSL <- predict(rsl, train)
  predTrainRSL <- .probabilisticToCrispData(rsl, predTrainRSL)
  accTrainRSL <- accuracy(predTrainRSL, trainActual)
  hamTrainRSL <- hammingLoss(predTrainRSL, trainActual)
  likTrainRSL <- .avgLogLikelihood(rsl, train, trainActual)
  # cat("Predicting train2...\n")
  # predTrainMAP <- predict(rsl, train, type = "joint")
  # predTrainMAP <- .probabilisticToCrispData(rsl, predTrainMAP)
  # accTrainMAP <- accuracy(predTrainMAP, trainActual)
  
  cat("Predicting on val...\n")
  predValRSL <- predict(rsl, val)
  predValRSL <- .probabilisticToCrispData(rsl, predValRSL)
  accValRSL <- accuracy(predValRSL, valActual)
  hamValRSL <- hammingLoss(predValRSL, valActual)
  likValRSL <- .avgLogLikelihood(rsl, val, valActual)
  # predValMAP <- predict(rsl, val, type = "joint")
  # predValMAP <- .probabilisticToCrispData(rsl, predValMAP)
  # accValMAP <- accuracy(predValMAP, valActual)
  
  cat("Predicting on test...\n")
  predTestRSL <- predict(rsl, test)
  predTestRSL <- .probabilisticToCrispData(rsl, predTestRSL)
  accTestRSL <- accuracy(predTestRSL, testActual)
  hamTestRSL <- hammingLoss(predTestRSL, testActual)
  likTestRSL <- .avgLogLikelihood(rsl, test, testActual)
  # predTestMAP <- predict(rsl, test, type = "joint")
  # predTestMAP <- .probabilisticToCrispData(rsl, predTestMAP)
  # accTestMAP <- accuracy(predTestMAP, testActual)
  
  return(list(accTrain = accTrainRSL, hamTrain = hamTrainRSL, logLikTrain = likTrainRSL,
              accVal = accValRSL, hamVal = hamValRSL, logLikVal = likValRSL,
              accTest = accTestRSL, hamTest = hamTestRSL, logLikTest = likTestRSL))
}


# .buildRSL - builds an rsl given a ruleset and a labelset
.buildRSL <- function(labels, rules){
  rsl <- createRSL()
  for(i in seq(along = labels)){
    rsl <- addClassifier(rsl, names(labels)[i], labels[[i]], accuracy = 1)
  }
  for(i in seq(along = rules)){
    probs <- rules[[i]]$p
    names(probs) <- rules[[i]]$labels
    rsl <- .addNoisyOR(rsl, probs)
  }
  
  return(rsl)
}


# generateDataset - creates an rsl and samples a dataset from it
# Input:
#  nTrain - integer, number of training observations
#  nVal - integer, number of validation observations
#  nTest - integer, number of test obserations
#  nLabels - integer, how many labelsets should be built
#  maxCategories - integer, how many labels may at most be in each labelset
#  nRules - integer, how many rules should be applied on the dataset
#  maxLabelsPerRule - integer, how long may the rules be (how many labels)
#  onlyPositiveRules - boolean, should rules always have p > 0.5
#  alpha - parameter for dirichlet distribution, used to generate the probabilities
#          per label. Increase to get them closer to actual label, decrease to
#          move them further away from the true label (keep > 1, otherwise 
#          the true label will become less likely than the false ones)
# Output:
#  list with 8 elements:
#    train - dataframe, the train dataset, probabilities per label
#    trainActual - dataframe, the true labels per train observation
#    val - dataframe, the validation dataset, probabilities per label
#    valActual - dataframe, the true labels per validation observation
#    test - dataframe, the test dataset, probabilities per label
#    testActual - dataframe, the true labels per test observation
#    labels - list of , each entry gives the labels in a labelset
#    rules - list of lists, each entry gives a rule and its probability
generateDataset<- function(nTrain = 4000, nVal = 2000, nTest = 3000,
                           nLabels = 10, maxCategories = 3, nRules = 20,
                           maxLabelsPerRule = 5,
                           alpha = 3){
  # TODO: Do we also want to simulate priors for the labels?
  
  # Globally for the dataset, generate labels and rules
  labels <- .generateLabelSet(nLabels, maxCategories)
  rules <- .generateRuleSet(labels, nRules, maxLabelsPerRule)
  rsl <- .buildRSL(labels, rules)
  
  # Generate priors
  priors <- .generateDirichletPriors(rsl, nTrain + nVal + nTest, alpha)
    
  # Generate fitting actuals from the true model
  trueLabels <- predict.rsl(rsl, priors)
  trueLabels <- .probabilisticToCrispData(rsl, trueLabels)
  
  # split into train, val, test
  cur <- 1
  train <- NULL
  trainActual <- NULL
  if(nTrain > 0){
    train <- priors[seq(cur, nTrain), ]
    trainActual <- trueLabels[seq(cur, nTrain), ]
    cur <- cur + nTrain
  }
  
  val <- NULL
  valActual <- NULL
  if(nVal > 0){
    val <- priors[seq(cur, nTrain + nVal), ]
    valActual <- trueLabels[seq(cur, nTrain + nVal), ]
    cur <- cur + nVal
  }
  
  test <- NULL
  testActual <- NULL
  if(nTest > 0){
    test <- priors[seq(cur, nTrain + nVal + nTest), ]
    testActual <- trueLabels[seq(cur, nTrain + nVal + nTest), ]
  }
  
  return(list(train = train, trainActual = trainActual, 
              val = val, valActual = valActual, 
              test = test, testActual = testActual, 
              labels = labels, rules = rules))
}

# evalTrueRSL - evaluates hamming loss and accuracy of the data generating rsl
# Input:
#  dataset - object generated by generateDataset
# Output:
#  list of 15 elements:
#    accTrain, accVal, accTest - accuracy on train, val and test when using the
#                                raw probabilities as predictions
#    hamTrain, hamVal, hamTest - hamming loss on train, val and test when using
#                                the raw probabilities as predictions
#    accTrainRSL, etc. - accuracy on train etc. when running the probabilities
#                        through the rsl and using the marginals as prediction
#    acc2TrainRSL, etc. - accuracy on train etc. when running the probabilities
#                         through the rsl and using the MAP as prediction
#    hamTrainRSL, etc. - hamming loss on train etc. when running the probabilities
#                        through the rsl and using the marginals as prediction
evalTrueRSL <- function(dataset){
  rsl <- .buildRSL(dataset$labels, dataset$rules)
  return(.evalPerformance(rsl, dataset$train, dataset$trainActual,
                    dataset$val, dataset$valActual,
                    dataset$test, dataset$testActual))
}

# evalEmptyRSL - evaluates hamming loss and accuracy of the baseline (without 
#                rules; using only the priors to predict)
evalEmptyRSL <- function(dataset){
  rsl <- .buildRSL(dataset$labels, list())
  return(.evalPerformance(rsl, dataset$train, dataset$trainActual,
                    dataset$val, dataset$valActual,
                    dataset$test, dataset$testActual))
}

# evalLearnedRSL - learns an RSL on the data and evaluates it
evalLearnedRSL <- function(dataset, nRules = 20, method = "noisyor", maxIter = 50, 
                           batchsize = 20, alpha = 0.002, beta1 = 0.9, 
                           beta2 = 0.999, eps = 1e-8){
  rsl <- .buildRSL(dataset$labels, list())
  rsl <- learnRules(rsl, dataset$train, dataset$trainActual, nRules, maxIter = maxIter,
                    batchsize = batchsize, alpha = alpha, beta1 = beta1,
                    beta2 = beta2, eps = eps, method = method)
  return(.evalPerformance(rsl, dataset$train, dataset$trainActual,
                    dataset$val, dataset$valActual,
                    dataset$test, dataset$testActual))
}

# Generate small simulation datasets
set.seed(25112020)
for(i in 1:10){
  cat(i, "...\n")
  data <- generateDataset(5000, 2000, 3000, nLabels = 10, maxCategories = 4, nRules = 10, maxLabelsPerRule = 5, alpha = 1)
  save(data, file = paste0("simData_", i, ".RData"))
}

set.seed(25112020)
for(i in 1:10){
  cat(i, "...\n")
  load(paste0("simData_", i, ".RData"))
  res <- evalEmptyRSL(data)
  save(res, file = paste0("res_emptyRSL_simData", i, ".RData"))
}

set.seed(25112020)
for(i in 1:10){
  cat(i, "...\n")
  load(paste0("simData_", i, ".RData"))
  res <- evalTrueRSL(data)
  save(res, file = paste0("res_trueRSL_simData", i, ".RData"))
}

# 
# resEmptyRSL <- evalEmptyRSL(data)
# resTrueRSL <- evalTrueRSL(data)
# resLearnedRSL <- evalLearnedRSL(data)
# cbind(resEmptyRSL, resTrueRSL, resLearnedRSL)

# Hyperparameter tuning:
# params <- expand.grid(list(alpha = c(0.01, 0.05, 0.001, 0.0005, 0.0001),
#                            beta1 = c(0.95, 0.9, 0.85, 0.8),
#                            beta2 = c(0.99, 0.999, 0.9999),
#                            eps = c(1e-6, 1e-8, 1e-10),
#                            maxIter = 50))
# allRes <- matrix(NA, nrow = nrow(params), ncol = 6)
# colnames(allRes) <- c("accTrain", "hamTrain", "accVal", "hamVal", "accTest", "hamTest")
# for(i in seq(nrow(params))){
#   cat(i, "\n", unlist(params[i, ]))
#   allRes[i, ] <- unlist(evalLearnedRSL(data, 20, "noisyor", maxIter = params$maxIter,
#                                        batchsize = 20, alpha = params$alpha,
#                                        beta1 = params$beta1, beta2 = params$beta2,
#                                        eps = params$eps))
# }
# 
# paramsTransformed <- params
# paramsTransformed$alpha <- 2 * (paramsTransformed$alpha - 0.0001) / (0.05 - 0.0001) - 1
# paramsTransformed$beta1 <- 2 * (paramsTransformed$beta1 - 0.8) / (0.95 - 0.8) - 1
# paramsTransformed$beta2[paramsTransformed$beta2 == 0.99] <- -1
# paramsTransformed$beta2[paramsTransformed$beta2 == 0.999] <- 0
# paramsTransformed$beta2[paramsTransformed$beta2 == 0.9999] <- 1
# paramsTransformed$eps[paramsTransformed$eps == 1e-6] <- 1
# paramsTransformed$eps[paramsTransformed$eps == 1e-8] <- 0
# paramsTransformed$eps[paramsTransformed$eps == 1e-10] <- -1
# dat <- cbind(paramsTransformed, y = allRes[, "hamTrain"])
# summary(lm(y ~ alpha^2 + beta1^2 + beta2 + eps + 1, data = dat))
