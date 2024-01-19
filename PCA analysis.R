my_data <- Book1_HTD
head(my_data, 6)

library("factoextra")
library("FactoMineR")
pca.data <- PCA(my_data[,-1], scale.unit = TRUE, graph = FALSE)
fviz_eig(pca.data, addlabels = TRUE, ylim = c(0, 70))
fviz_pca_var(pca.data, col.var = "cos2",
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"),
             repel = TRUE)

pca.data <- PCA(t(my_data[,-1]), scale.unit = TRUE, graph = FALSE)

fviz_pca_ind(pca.data, col.ind = "cos2", 
             gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
             repel = TRUE)

library(ggpubr)
a <- fviz_pca_ind(pca.data, col.ind = "cos2", 
                  gradient.cols = c("#FFCC00", "#CC9933", "#660033", "#330033"), 
                  repel = TRUE)

ggpar(a,
      title = "Principal Component Analysis",
      xlab = "PC1", ylab = "PC2",
      legend.title = "Cos2", legend.position = "top",
      ggtheme = theme_minimal())

fviz_cos2(pca.data, choice = "ind")
fviz_contrib(pca.data, choice = "ind", axes = 1:2)

fviz_pca_ind(pca.data, axes = c(2, 3))