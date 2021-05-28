# This script compares the exact with the approximate estimates

# Dependencies:
source("rsl.R")
source("noisyornetwork.R")
# And the simData_k.RData have to be located in the directory

# .buildRSL - builds an rsl given a ruleset and a labelset
.buildRSL <- function(labels, rules){
  rsl <- createRSL()
  for(i in seq(along = labels)){
    rsl <- addClassifier(rsl, names(labels)[i], labels[[i]], accuracy = 1)
  }
  norn <- as.norn.rsl(rsl)
  for(i in seq(along = rules)){
    probs <- rules[[i]]$p
    names(probs) <- rules[[i]]$labels
    probList <- .preprocessInhProbs(rsl, probs)
    rsl <- .addNoisyOR(rsl, probs)
    norn <- addRule.norn(norn, probList)
  }
  
  return(list(rsl = rsl, norn = norn))
}

set.seed(25112020)
pdf("small.pdf", width = 8, height = 16)
par(mfrow = c(5, 2), mar = c(5, 4, 1, 2))
for(i in 1:10){
  cat(i, "...\n")
  load(paste0("../small/simData/simData_", i, ".RData"))
  net <- .buildRSL(data$labels, data$rules)
  
  data$train <- data$train[1:100, ]
  predExact <- predict(net$rsl, data$train)
  predApprox <- predict.norn(net$norn, net$rsl, data$train, showProgress = TRUE)
  
  plot(unlist(predExact), unlist(predApprox), col = adjustcolor("black", alpha.f = 0.05),
       xlab = "Exact Marginal", ylab = "Approximate Marginal", sub = paste("Dataset", i))
}
dev.off()