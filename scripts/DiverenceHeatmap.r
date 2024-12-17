png("divergence_fst_heatmap.png", width=1000, height=1000)

require(ComplexHeatmap)
load("Fstats_all.RData") # -> all.pops.fstats 
div.hm <- Heatmap(all.pops.fstats@pairwise.div, cluster_rows = TRUE, cluster_columns=TRUE, name="Divergence (pairwise)", show_heatmap_legend=FALSE, column_title="Divergence (1-Q2)")
fst.hm <- Heatmap(all.pops.fstats@pairwise.fst, cluster_rows=TRUE, cluster_columns=TRUE, name="values", column_title = "fst=(Q1-Q2)/(1-Q2)")
div.hm+fst.hm
dev.off()