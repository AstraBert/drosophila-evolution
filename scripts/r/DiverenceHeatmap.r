png("divergence_fst_heatmap.png", width=1000, height=1000)

require(ComplexHeatmap)
load("/gatk_modified/userdata/abertelli/drosophila-evolution/results/all_pops_data.RData") # -> all.pops.fstats 
poolcsv <- read.csv("pools.csv")
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