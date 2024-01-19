y <- DGEList(counts = counts, group = Tissue)

cpm <- cpm(y)
lcpm <- cpm(y, log=TRUE)

L <- mean(y$samples$lib.size) * 1e-6
M <- median(y$samples$lib.size) * 1e-6
c(L, M)

table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=Tissue)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

samplenames <- substring(colnames(y), 1, nchar(colnames(y)))
samplenames

lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(y)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(y, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
