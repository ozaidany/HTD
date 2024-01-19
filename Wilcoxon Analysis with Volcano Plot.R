library("edgeR")
my_data <- Book1_HTD
head(my_data, 6)

counts <- my_data[,-1]

thresh <- my_data[, 2:26] > 200 & my_data[, 27:51] > 200
head(thresh)
keep <- rowSums(thresh) >= 5
summary(keep)

my_data1 <- my_data[keep,]
View(my_data1)

thresh1 <- my_data1[,2:51] <= 100
head(thresh1)
keep1 <- rowSums(thresh1) == 0
summary(keep1)

my_data2 <- my_data1[keep1,]
View(my_data2)

counts <- my_data2[,-1]

Tissue <- factor(c("N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","N","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D","D"))

y <- DGEList(counts = counts, group = Tissue)

keep <- filterByExpr(y)

y <- y[keep, keep.lib.sizes = FALSE]

y <- calcNormFactors(y, method = "TMM")

count_norm <- cpm(y)

count_norm <- as.data.frame(count_norm)

my_data3 = cbind(my_data2[,1], count_norm)
my_data3 = as.data.frame(my_data3)
View(my_data3)

count_norm <- my_data3[,-1]
count_norm <- as.data.frame(count_norm)

countCV <- my_data3[,-1]
gene.mean <- apply(countCV,1,mean)
gene.sd <- apply(countCV,1,sd)
gene.cv <- gene.sd/gene.mean
summary(gene.cv)

thresh3 <- gene.cv >= 0.4
head(thresh3)
summary(thresh3)

my_data4 <- my_data3[thresh3,]
View(my_data4)

count_norm <- my_data4[,-1]
count_norm <- as.data.frame(count_norm)

pvalues <- sapply(1:nrow(count_norm), function(i){
data <- cbind.data.frame(gene = as.numeric(t(count_norm[i,])), Tissue)
p <- wilcox.test(gene~Tissue, data)$p.value
return(p)
})

fdr <- p.adjust(pvalues, method = "fdr")

conditionsLevel <- levels(Tissue)

dataCon1 <- count_norm[,c(which(Tissue==conditionsLevel[1]))]

dataCon2 <- count_norm[,c(which(Tissue==conditionsLevel[2]))]

foldChanges <- log2(rowMeans(dataCon2)/rowMeans(dataCon1))

outRst <- data.frame(log2foldChange = foldChanges, pValues = pvalues, FDR = fdr)

rownames(outRst) <- rownames(count_norm)

outRst <- na.omit(outRst)

fdrThres <- 0.001

write.table(outRst[outRst$FDR<fdrThres,], file = "C:/Users/31683/OneDrive/Desktop/HTD/sth.xlsx", sep="\t", quote = F, row.names = T, col.names = T)

results = cbind(my_data4[,1,drop=F], foldChanges, pvalues, fdr)
results = as.data.frame(results)
View(results)

library("ggplot2")
p1 <- ggplot(results, aes(foldChanges, -log(fdr,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

library("magrittr") # needs to be run every time you start R and want to use %>%
library("dplyr")    # alternatively, this also loads %>%
library("hms")
library("kableExtra")
library("tidyverse")
library("ggrepel")

# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

results <- results %>% 
  mutate(
    Expression = case_when(foldChanges >= log(2) & fdr <= 0.05 ~ "Up-regulated",
                           foldChanges <= -log(2) & fdr <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(results) %>% 
  knitr_table()
View(results)

p2 <- ggplot(results, aes(foldChanges, -log(fdr,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

results %>% 
  dplyr::count(Expression) %>% 
  knitr_table()

results <- results %>% 
  mutate(
    Significance = case_when(
      abs(foldChanges) >= log(2) & fdr <= 0.05 & fdr > 0.01 ~ "FDR 0.05", 
      abs(foldChanges) >= log(2) & fdr <= 0.01 & fdr > 0.001 ~ "FDR 0.01",
      abs(foldChanges) >= log(2) & fdr <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )
head(results) %>% 
  knitr_table()

p3 <- ggplot(results, aes(foldChanges, -log(fdr,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

results %>% 
  dplyr::count(Expression, Significance) %>% 
  knitr_table()

top <- 5
top_genes <- bind_rows(
  results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(fdr, desc(abs(foldChanges))) %>% 
    head(top),
  results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(fdr, desc(abs(foldChanges))) %>% 
    head(top)
)
top_genes %>% 
  knitr_table()

p3 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(foldChanges, -log(fdr,10), label = GeneID),
                   size = 2)
p3

sessionInfo()