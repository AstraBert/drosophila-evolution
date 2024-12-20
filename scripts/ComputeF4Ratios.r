all.popswsim.data <- vcf2pooldata(vcf.file="/gatk_modified/userdata/abertelli/drosophila-evolution/results/drosophila_evolution.bcftools_fakepools_withsim.vcf.gz", min.maf="0.05", poolsizes=rep(100000,21))
save(all.popswsim.data, file="all_popswsim_data.RData")
all.popswsim.fstats <- compute.fstats(all.popswsim.data, nsnp.per.bjack.block = 1000, computeDstat = TRUE,verbose=FALSE)
save(all.popswsim.fstats, file="all_popswsim_fstats.RData")

# Load the CSV file into a data frame (assuming tab-separated file)
df <- read.table("F4_to_compute.tsv", header = TRUE, sep = "\t")

# Loop through each row of the data frame
for (i in 1:nrow(df)) {
  # Extract NUM and DEN values for the current row
  num_quadruplets <- df$NUM[i]
  den_quadruplets <- df$DEN[i]
  
  # Construct the compute.f4ration command
  command <- paste0("compute.f4ratio(all.popswsim.fstats, num.quadruplets='", num_quadruplets, "', den.quadruplets='", den_quadruplets, "')")
  
  # Print the command (optional, to check what is being run)
  print(command)
  
  # Run the command (eval() will evaluate the constructed string as R code)
  eval(parse(text = command))
}