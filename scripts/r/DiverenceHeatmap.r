require(poolfstat)

png("divergence_fst_heatmap.png", width=1000, height=1000)

require(ComplexHeatmap)
load("results/all_pops_data_wosim.RData") # -> all.pops.fstats 
all.pops.fstats <- compute.fstats(all.pops.data, nsnp.per.bjack.block=10000)
save(all.pops.fstats, file="results/all_pops_fstats_wosim.RData")
poolcsv <- read.csv("pools1.csv")
hm.fst <- all.pops.fstats@pairwise.fst
rownames(hm.fst) <- poolcsv$NAME
colnames(hm.fst) <- poolcsv$NAME

hm.div <- all.pops.fstats@pairwise.div
rownames(hm.div) <- poolcsv$NAME
colnames(hm.div) <- poolcsv$NAME

div.hm <- Heatmap(hm.div, cluster_rows = TRUE, cluster_columns=TRUE, name="Divergence (pairwise)", show_heatmap_legend=FALSE, column_title="Divergence (1-Q2)")
fst.hm <- Heatmap(hm.fst, cluster_rows=TRUE, cluster_columns=TRUE, name="values", column_title = "fst=(Q1-Q2)/(1-Q2)")
div.hm+fst.hm
dev.off()
