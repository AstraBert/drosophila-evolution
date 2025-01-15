# Read the text file after the first two lines
dfpool <- read.table("POOL_PopData.snp", skip = 2, header = FALSE, nrows = 19998)

# Check the dimension of dfpool
cat("Initial dimension of dfpool:", dim(dfpool), "\n")

# Define the minimum number of reads for each SNP in each of the 4 populations
n.reads.min = 20
n.reads.max = 500

# Compute the sums for each row
dfpool$sum12 <- dfpool$V1 + dfpool$V2
dfpool$sum34 <- dfpool$V3 + dfpool$V4
dfpool$sum56 <- dfpool$V5 + dfpool$V6
dfpool$sum78 <- dfpool$V7 + dfpool$V8

# Rename to dfpoolsum (optional, the dataframe is already modified)
dfpoolsum <- dfpool

# Remove rows based on the criterion (we need at least 50X, so n.reads.min = 50 reads for each SNP in all 4 populations)
dfpoolsel <- dfpoolsum[!(dfpoolsum$sum12 < n.reads.min | dfpoolsum$sum34 < n.reads.min | dfpoolsum$sum56 < n.reads.min | dfpoolsum$sum78 < n.reads.min | dfpoolsum$sum12 > n.reads.max | dfpoolsum$sum34 > n.reads.max | dfpoolsum$sum56 > n.reads.max | dfpoolsum$sum78 > n.reads.max),]

# Check the dimension of dfpoolsel
cat("Dimension after selection of dfpoolsel:", dim(dfpoolsel), "\n")

# Display the first 10 rows of dfpoolsel
cat("The first 10 rows of dfpoolsel:\n")
print(head(dfpoolsel, 10))

# Select the first 6 columns of dfpoolsel
df_poolseq_readmin <- dfpoolsel[, 1:8]

# Check the dimension of df_poolseq_readmin
cat("Dimension of df_poolseq_readmin:", dim(df_poolseq_readmin), "\n")

# Output file name
output.file.name = paste0("POOL_PopData_df_poolseq_readmin_",n.reads.min,".txt")

# Write the dataframe to a text file
write.table(
  df_poolseq_readmin,
  file = output.file.name,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Check the dimension of dfpoolsel
cat("Dimension after selection of dfpoolsel:", dim(df_poolseq_readmin), "\n")
cat("Number of selected SNPs =", nrow(df_poolseq_readmin), "\n")

# Confirmation
cat("The file",output.file.name, "was successfully created", "\n")
