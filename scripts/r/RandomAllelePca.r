load("all_pops_data.RData") # -> all.pops.data

poolcsv <- read.csv("pools.csv")

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