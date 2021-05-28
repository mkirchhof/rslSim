# This script compares the exact with the approximate estimates

# Dependencies:
source("rsl.R")
library("parallel")
library(microbenchmark)
# And the simData_k.RData have to be located in the directory

# .buildRSL - builds an rsl given a ruleset and a labelset
.buildRSL <- function(labels, rules){
  rsl <- createRSL()
  for(i in seq(along = labels)){
    rsl <- addClassifier(rsl, names(labels)[i], labels[[i]], accuracy = 1)
  }
  for(i in seq(along = rules)){
    probs <- rules[[i]]$p
    names(probs) <- rules[[i]]$labels
    probList <- .preprocessInhProbs(rsl, probs)
    rsl <- .addNoisyOR(rsl, probs)
  }
  
  return(rsl)
}

ntasks <- 10
cl <- makeCluster(ntasks)
clusterEvalQ(cl, source("rsl.R"))

cat("Small:\n")
for(i in 1:10){
  cat(i, "...\n")
  load(paste0("../small/simData/simData_", i, ".RData"))
  rsl <- .buildRSL(data$labels, data$rules)
  data$train <- data$train[1:100, ]
  
  clusterSetRNGStream(cl, iseed=22122020)
  set.seed(22122020)
  exactMargTime <- microbenchmark(exactMarg <- predict.rsl(rsl, data$train, method = "exact", type = "marginal", cluster = cl), times = 1)$time / 100 * ntasks
  approxMargTime <- microbenchmark(approxMarg <- predict.rsl(rsl, data$train, method = "approximate", type = "marginal", cluster = cl), times = 1)$time / 100 * ntasks
  exactJointTime <- microbenchmark(exactJoint <- predict.rsl(rsl, data$train, method = "exact", type = "joint", cluster = cl), times = 1)$time / 100 * ntasks
  approxJointTime <- microbenchmark(approxJoint <- predict.rsl(rsl, data$train, method = "approximate", type = "joint", cluster = cl), times = 1)$time / 100 * ntasks
  
  res <- list(exactMarg = exactMarg, approxMarg = approxMarg, exactJoint = exactJoint, approxJoint = approxJoint,
              exactMargTime = exactMargTime, approxMargTime = approxMargTime, exactJointTime = exactJointTime, approxJointTime = approxJointTime)
  save(res, file = paste0("res_small_", i, ".RData"))
}

cat("Medium:\n")
for(i in 1:10){
  cat(i, "...\n")
  load(paste0("../medium/simData/simData_", i, ".RData"))
  rsl <- .buildRSL(data$labels, data$rules)
  data$train <- data$train[1:100, ]
  
  clusterSetRNGStream(cl, iseed=22122020)
  set.seed(22122020)
  exactMargTime <- microbenchmark(exactMarg <- predict.rsl(rsl, data$train, method = "exact", type = "marginal", cluster = cl), times = 1)$time / 100 * ntasks
  approxMargTime <- microbenchmark(approxMarg <- predict.rsl(rsl, data$train, method = "approximate", type = "marginal", cluster = cl), times = 1)$time / 100 * ntasks
  exactJointTime <- microbenchmark(exactJoint <- predict.rsl(rsl, data$train, method = "exact", type = "joint", cluster = cl), times = 1)$time / 100 * ntasks
  approxJointTime <- microbenchmark(approxJoint <- predict.rsl(rsl, data$train, method = "approximate", type = "joint", cluster = cl), times = 1)$time / 100 * ntasks
  
  res <- list(exactMarg = exactMarg, approxMarg = approxMarg, exactJoint = exactJoint, approxJoint = approxJoint,
              exactMargTime = exactMargTime, approxMargTime = approxMargTime, exactJointTime = exactJointTime, approxJointTime = approxJointTime)
  save(res, file = paste0("res_medium_", i, ".RData"))
}

cat("almostLarge:\n")
for(i in 1:10){
  cat(i, "...\n")
  load(paste0("../almostLarge/simData/simData_", i, ".RData"))
  rsl <- .buildRSL(data$labels, data$rules)
  data$train <- data$train[1:100, ]
  
  clusterSetRNGStream(cl, iseed=22122020)
  set.seed(22122020)
  exactMargTime <- microbenchmark(exactMarg <- predict.rsl(rsl, data$train, method = "exact", type = "marginal", cluster = cl), times = 1)$time / 100 * ntasks
  approxMargTime <- microbenchmark(approxMarg <- predict.rsl(rsl, data$train, method = "approximate", type = "marginal", cluster = cl), times = 1)$time / 100 * ntasks

  res <- list(exactMarg = exactMarg, approxMarg = approxMarg, 
              exactMargTime = exactMargTime, approxMargTime = approxMargTime)
  save(res, file = paste0("res_almostLarge_", i, ".RData"))
}

stopCluster(cl)
