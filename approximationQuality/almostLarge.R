# This script compares the exact with the approximate estimates

# Dependencies:
source("rsl.R")
source("noisyornetwork.R")
library(tikzDevice)
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

pdf("almostLarge.pdf", width = 8, height = 16)
par(mfrow = c(5, 2), mar = c(5, 4, 1, 2))
for(i in 4){
  cat(i, "...\n")
  load(paste0("../almostLarge/simData/simData_", i, ".RData"))
  set.seed(25112020)
  net <- .buildRSL(data$labels, data$rules)
  
  data$train <- data$train[1:100, ]
  predExact <- predict(net$rsl, data$train, showProgress = TRUE)
  predApprox <- predict.norn(net$norn, net$rsl, data$train, showProgress = TRUE)
  
  save(predExact, predApprox, file = paste0("pred_", i, ".RData"))
  
  plot(unlist(predExact), unlist(predApprox), col = adjustcolor("black", alpha.f = 0.05),
       xlab = "Exact Marginal", ylab = "Approximate Marginal", sub = paste("Dataset", i))
}
dev.off()

# Plot for paper:
load("pred_4.RData")
tikz('predQuality.tex', standAlone = FALSE, width=3.3, height=2.7)
par(mar = c(4, 5.7, 0, 0))
plot(unlist(predExact), unlist(predApprox), col = adjustcolor("black", alpha.f = 0.01), 
     ylab = "Approximate Marginal", xlab = "",
     pch = 16)
mtext("Exact Marginal", side = 1, line = 2.5)
dev.off()