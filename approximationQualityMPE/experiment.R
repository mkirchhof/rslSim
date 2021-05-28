# Get the ID that SLURM hands over to the script as argument
folds <- as.integer(Sys.getenv("PBS_ARRAYID"))
if(is.na(folds)){
  folds <- 1:10
}
cat("ID = ", paste(folds, collapse = ", "), "\n")


# Dependencies:
source("rsl.R")

# And the simData_k.RData have to be located in the directory

# evalPerformance - evaluates hamming loss, accuracy and log-likelihood on 
#                   test, validation and train data and the time it took per 
#                   sample for a prediction on the test dataset
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
  predTime <- microbenchmark(predTestRSL <- predict(rsl, test), times = 1)$time / nrow(test)
  predTestRSL <- .probabilisticToCrispData(rsl, predTestRSL)
  accTestRSL <- accuracy(predTestRSL, testActual)
  hamTestRSL <- hammingLoss(predTestRSL, testActual)
  likTestRSL <- .avgLogLikelihood(rsl, test, testActual)
  # predTestMAP <- predict(rsl, test, type = "joint")
  # predTestMAP <- .probabilisticToCrispData(rsl, predTestMAP)
  # accTestMAP <- accuracy(predTestMAP, testActual)
  
  return(list(accTrain = accTrainRSL, hamTrain = hamTrainRSL, logLikTrain = likTrainRSL,
              accVal = accValRSL, hamVal = hamValRSL, logLikVal = likValRSL,
              accTest = accTestRSL, hamTest = hamTestRSL, logLikTest = likTestRSL,
              avgPredTime = predTime))
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

# Main:
for(k in folds){
  cat(k, "\n")
  # Learn an rsl
  # load an object called data containing all important data
  load(paste0("../medium/simData/simData_", k, ".RData"))
  
  # build the true RSL
  rsl <- .buildRSL(data$labels, data$rules)
  
  # Compute the marginals in correct, marginal and approximate way
  dat <- data$test[1:100, ]
  truePred <- predict.rsl(rsl, dat, method = "exact", type = "joint", showProgress = TRUE)
  truePred <- .probabilisticToCrispData(rsl, truePred)
  margPred <- predict.rsl(rsl, dat, method = "exact", type = "marginal")
  margPred <- .probabilisticToCrispData(rsl, margPred)
  approxPred <- predict.rsl(rsl, dat, method = "approximate", type = "joint", showProgress = TRUE)
  approxPred <- .probabilisticToCrispData(rsl, approxPred)
  
  # Calculate % correct
  accMarg <- accuracy(margPred, truePred)
  accApprox <- accuracy(approxPred, truePred)
  res <- c(marg = accMarg, approx = accApprox)
  
  save(res, file = paste0(k, "_res.RData"))
}
