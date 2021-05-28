# Summarizes the results of the approximate loopy belief propagation experiment

# Dependencies:
source("rsl.R")
library(tikzDevice)


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


# Small
corMarg <- numeric(10)
precJointSmall <- numeric(10)
precMargSmall <- numeric(10)
for(i in 1:10){
  load(paste0("res_small_", i,".RData"))
  
  # Marginals
  plot(unlist(res$exactMarg), unlist(res$approxMarg), col = adjustcolor("black", alpha.f = 0.01), 
       ylab = "Approximate Marginal", xlab = "Exact Marginal",
       pch = 16, sub = i)
  
  corMarg[i] <- cor(unlist(res$exactMarg), unlist(res$approxMarg))
  
  # MPEs
  load(paste0("../small/simData/simData_", i,".RData"))
  rsl <- .buildRSL(data$labels, data$rules)
  truePred <- .probabilisticToCrispData(rsl, res$exactJoint)
  margPred <- .probabilisticToCrispData(rsl, res$approxMarg)
  approxPred <- .probabilisticToCrispData(rsl, res$approxJoint)
  
  # Calculate % correct
  precJointSmall[i] <- accuracy(approxPred, truePred)
  precMargSmall[i] <- accuracy(margPred, truePred)
}

mean(corMarg) # 0.9949422
precJointSmall # 0.98 0.89 0.98 1.00 0.94 0.94 0.98 0.94 0.95 0.98
precMargSmall # 0.71 0.47 0.90 0.80 0.70 0.58 0.66 0.63 0.72 0.41

# Medium
corMarg <- numeric(10)
precJointMedium <- numeric(10)
precMargMedium <- numeric(10)
for(i in 1:10){
  load(paste0("res_medium_", i,".RData"))
  
  # Marginals
  plot(unlist(res$exactMarg), unlist(res$approxMarg), col = adjustcolor("black", alpha.f = 0.01), 
       ylab = "Approximate Marginal", xlab = "Exact Marginal",
       pch = 16, sub = i)
  
  corMarg[i] <- cor(unlist(res$exactMarg), unlist(res$approxMarg))
  
  # MPEs
  load(paste0("../medium/simData/simData_", i,".RData"))
  rsl <- .buildRSL(data$labels, data$rules)
  truePred <- .probabilisticToCrispData(rsl, res$exactJoint)
  margPred <- .probabilisticToCrispData(rsl, res$approxMarg)
  approxPred <- .probabilisticToCrispData(rsl, res$approxJoint)
  
  # Calculate % correct
  precJointMedium[i] <- accuracy(approxPred, truePred)
  precMargMedium[i] <- accuracy(margPred, truePred)
}

mean(corMarg) # 0.9955726
precJointMedium # 0.90 0.91 0.98 0.89 0.95 0.99 0.96 0.86 0.94 0.94
precMargMedium # 0.42 0.37 0.69 0.36 0.35 0.42 0.59 0.50 0.51 0.28

# almostLarge
corMarg <- numeric(10)
precJointAlmostLarge <- numeric(10)
for(i in 1:10){
  load(paste0("res_almostLarge_", i,".RData"))
  
  plot(unlist(res$exactMarg), unlist(res$approxMarg), col = adjustcolor("black", alpha.f = 0.01), 
       ylab = "Approximate Marginal", xlab = "Exact Marginal",
       pch = 16, sub = i)
  
  corMarg[i] <- cor(unlist(res$exactMarg), unlist(res$approxMarg))
}

mean(corMarg) # 0.9950552


# Plot MPE for small and medium
tikz('boxplots.tex', standAlone = FALSE, width=3, height=2)
par(mar = c(2.2, 2.8, 1.1, 0), mfrow = c(2, 1))
boxplot(precMargSmall, precJointSmall,
        ylim = c(0, 1), yaxt = "n", xaxt = "n", horizontal = TRUE)
axis(2, at = 1:2, labels = c("Marg", "MPE"), las = 2)
axis(1, at = seq(0, 1, 0.25), labels = c("", "", "", "", ""))
mtext(side = 3, line = 0.2, "Small")
par(mar = c(3.3, 2.8, 0.0, 0))
boxplot(precMargMedium, precJointMedium,
        ylim = c(0, 1), yaxt = "n", xaxt = "n", horizontal = TRUE)
axis(2, at = 1:2, labels = c("Marg", "MPE"), las = 2)
axis(1, at = seq(0, 1, 0.25), labels = paste0(seq(0, 1, 0.25)))
mtext(text = "Predictions equal to exact MPE", side = 1, line = 2.2)
mtext(side = 3, line = 0.2, "Medium")
dev.off()
