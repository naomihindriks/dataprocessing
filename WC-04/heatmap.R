args = commandArgs(trailingOnly=TRUE)
data <-read.csv(file = args[1], header=FALSE, sep=",")


d <- as.matrix(data[-1,-1])
rownames(d) <- data[-1,1]
colnames(d) <- data[1,-1]

jpeg(args[2], width = 1500, height = 1500)

par(cex.main=3, oma = c(6,4,7,2) + 0.1,  mar = c(7, 4, 4, 2) + 0.1)

heatmap(
  d, 
  Colv = NA, Rowv = NA, 
  xlab = "", ylab = "", main = "Expression data of yeast",
  cexCol = 2.5
)

mtext("Time", side=1, line=9, cex = 2.5)
mtext("Gene", side=4, cex = 2.5)

dev.off()