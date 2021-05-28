# collects and visualizes the performances of marginal and approx joint predictions

library(tikzDevice)

allRes <- list()
for(i in 1:10){
  load(paste0(i, "_res.RData"))
  allRes[[i]] <- res
}
allRes <- do.call(rbind, allRes)

pdf("performance.pdf", width = 8, height = 5)
par(mar = c(4, 4, 1, 1))
plot(NA, xlim = c(1, 10), ylim = c(0, 1), xlab = "Dataset", ylab = "Accuracy")
points(y = rep(1, 10), x = 1:10, col = "black", pch = 16)
points(y = allRes[, 1], x = 1:10, col = "red", pch = 16)
points(y = allRes[, 2], x = 1:10, col = "blue", pch = 16)
legend("bottomright", legend = c("Exact MPE", "Approx. MPE", "Marginals"), col = c("black", "blue", "red"),
       pch = 16)
dev.off()

allRes <- allRes * 100
tikz('boxplots.tex', standAlone = FALSE, width=3.3, height=1.35)
par(mar = c(3.5, 5.7, 0, 0))
colnames(allRes) <- c("Marginals", "Approx. MPE")
boxplot(allRes, horizontal = TRUE, ylim = c(0, 100), las = 1, xaxt = "n")
axis(1, at = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "\\%"))
mtext(text = "Predictions equal to exact MPE", side = 1, line = 2.5)
dev.off()