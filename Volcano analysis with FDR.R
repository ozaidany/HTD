my_data <- Book1_HTD
counts <- my_data[,-1]

thresh <- my_data[, 2:26] > 200 & my_data[, 27:51] > 200
head(thresh)
keep <- rowSums(thresh) >= 5
summary(keep)

my_data1 <- my_data[keep,]
View(my_data1)

counts <- my_data1[,-1]

# The packages can be installed within R with:
# intstall.packages(c("tidyverse", "ggrepel"))

# Load packages
library(tidyverse)
library(ggrepel)

# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

# Import data
dim(my_data1)

install.packages("magrittr") # package installations are only needed the first time you use it
install.packages("dplyr")    # alternative installation of the %>%
library(magrittr) # needs to be run every time you start R and want to use %>%
library(dplyr)    # alternatively, this also loads %>%
install.packages("hms")
library("hms")
library("kableExtra")

head(my_data1) %>% 
  knitr_table()

ttestRat <- function(df, grp1, grp2) {
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)  
  results = t.test(x, y)
  results$p.value
}
rawpvalue = apply(counts, 1, ttestRat, grp1 = c(1:25), grp2 = c(26:50))
rawpvalue

hist(rawpvalue)

##transform our data into log2 base.
logdata = log2(counts)

#calculate the mean of each gene per control group
control = apply(counts[,1:25], 1, mean)

#calcuate the mean of each gene per test group
test = apply(counts[, 26:50], 1, mean) 

#confirming that we have a vector of numbers
class(control)
#confirming that we have a vector of numbers
class(test)

foldchange <- control - test 

hist(foldchange, xlab = "log2 Fold Change (Control vs Test)")

FDR <- p.adjust(rawpvalue, method="BH")

results = cbind(my_data1[,1], foldchange, rawpvalue, FDR)
results = as.data.frame(results)

View(results)

p1 <- ggplot(results, aes(foldchange, -log(FDR,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR"))
p1

results <- results %>% 
  mutate(
    Expression = case_when(foldchange >= log(2) & FDR <= 0.05 ~ "Up-regulated",
                           foldchange <= -log(2) & FDR <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )

head(results) %>% 
  knitr_table()
View(results)

p2 <- ggplot(results, aes(foldchange, -log(FDR,10))) +
  geom_point(aes(color = Expression), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 
p2

results %>% 
  count(Expression) %>% 
  knitr_table()

results <- results %>% 
  mutate(
    Significance = case_when(
      abs(foldchange) >= log(2) & FDR <= 0.05 & FDR > 0.01 ~ "FDR 0.05", 
      abs(foldchange) >= log(2) & FDR <= 0.01 & FDR > 0.001 ~ "FDR 0.01",
      abs(foldchange) >= log(2) & FDR <= 0.001 ~ "FDR 0.001", 
      TRUE ~ "Unchanged")
  )
head(results) %>% 
  knitr_table()

p3 <- ggplot(results, aes(foldchange, -log(FDR,10))) +
  geom_point(aes(color = Significance), size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"FDR")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3

results %>% 
  count(Expression, Significance) %>% 
  knitr_table()

top <- 5
top_genes <- bind_rows(
  results %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(FDR, desc(abs(foldchange))) %>% 
    head(top),
  results %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(FDR, desc(abs(foldchange))) %>% 
    head(top)
)
top_genes %>% 
  knitr_table()

p3 <-  p3 +
  geom_label_repel(data = top_genes,
                   mapping = aes(foldchange, -log(FDR,10), label = GeneID),
                   size = 2)
p3

sessionInfo()