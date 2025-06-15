require("poolfstat")
require("mclust")
require("ggplot2")


load("all_pops_data.RData")
poolcsv <- read.csv("pools1.csv")
png("randomallele.png", width=1000, height=1000)
all.pops.pca <- randomallele.pca(all.pops.data)
pc1_2 <- all.pops.pca$pop.loadings[,c(1,2)]
rownames(pc1_2) <- poolcsv$NAME
perc.pc1 <- all.pops.pca$perc.var[[1]]
perc.pc2 <- all.pops.pca$perc.var[[2]]
perc.pc1.rounded <- round(perc.pc1, digits=2)
perc.pc1.str <- paste0(perc.pc1.rounded)
xlabel <- paste("PC1", "(", perc.pc1.str, "%)")
perc.pc2.rounded <- round(perc.pc2, digits=2)
perc.pc2.str <- paste0(perc.pc2.rounded)
ylabel <- paste("PC2", "(", perc.pc2.str, "%)")
plot(pc1_2[,1], pc1_2[,2], xlab = xlabel, ylab = ylabel, main = "Random Allele PCA", pch = 19, col = "red")
text(pc1_2[,1], pc1_2[,2], labels = rownames(pc1_2), cex = 0.8, pos = 4)
dev.off()
loadings <- all.pops.pca$pop.loadings
pc_variances <- apply(loadings, 2, var)
prop_variance <- pc_variances / sum(pc_variances)
cum_var <- cumsum(prop_variance)
plot(cum_var, type = "b", pch = 19, xlab = "Principal Component", ylab = "Cumulative Proportion of Variance Explained", main = "Scree Plot (Cumulative)")
abline(h = 0.8, col = "red", lty = 2)
# Limit to first 15 principal components
loadings_pc15 <- all.pops.pca$pop.loadings[, 1:15]
rownames(loadings_pc15) <- poolcsv$NAME
mclust_result <- Mclust(loadings_pc15)
plot(mclust_result, what = "classification")
# Limit to first 2 PCs
mclust_result1 <- Mclust(pc1_2)
plot(mclust_result1, what = "classification")
clusters <- mclust_result$classification
cluster_df <- data.frame(Pool = rownames(loadings_pc15),
                         Cluster = clusters)
write.csv(cluster_df, "pops_clusters_from_15.csv")
df_plot <- as.data.frame(loadings_pc15)
df_plot$Cluster <- as.factor(clusters)
df_plot$Pool <- rownames(df_plot)
# Plot nicely
p <- ggplot(df_plot, aes(x = V1, y = V2, color = Cluster, label = Pool)) +
  geom_point(size = 4) +
  geom_text(vjust = -1) +
  theme_minimal() +
  labs(title = "Clustering of Pools Based on First 15 PCs",
       x = "PC1", y = "PC2")
ggsave("clustering_plot.png", plot = p, width = 8, height = 6, dpi = 300)
# Plot PC1,2 clusters
plot(mclust_result1, what = "classification")
clusters <- mclust_result1$classification
cluster_df <- data.frame(Pool = rownames(loadings_pc15),
                         Cluster = clusters)
write.csv(cluster_df, "pops_clusters_from12.csv")
