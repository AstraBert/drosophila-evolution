import polars as pl

# Read CSV files
csv = pl.read_csv("f2_nohet.csv")
csv = csv.rename({"PoolName": "Pools"})
df = pl.read_csv("dest_samples.csv")

# Create mappings from pools to populations
pools = df["PoolName"].to_list()
poops = df["PopName"].to_list()
pools2pops = {pools[i]: poops[i] for i in range(len(pools))}

# Extract pools and map to populations
poolsf4 = csv["Pools"].to_list()
pops1 = [
    pools2pops[poolf4.split(",")[0]] if poolf4.split(",")[0] not in ["DGN", "CNXJ", "CnOther", "CnQTP", "ISR", "DrosSim"]
    else poolf4.split(",")[0]
    for poolf4 in poolsf4
]
pops2 = [
    pools2pops[poolf4.split(",")[1]] if poolf4.split(",")[1] not in ["DGN", "CNXJ", "CnOther", "CnQTP", "ISR", "DrosSim"]
    else poolf4.split(",")[1]
    for poolf4 in poolsf4
]

# Insert new columns into the CSV
csv.insert_column(1, pl.Series("Pop1", pops1))
csv.insert_column(2, pl.Series("Pop2", pops2))

# Define constants
A1 = "CnOther"
B3 = "DGN"
D2 = "DrosSim"
C4 = "whatever"
expression = f"0.5*({A1},{C4} + {D2},{B3} - {A1},{B3} - {D2},{C4})"

# Filter data based on conditions
drossim_whatever = csv.filter((pl.col("Pop1") == A1) | (pl.col("Pop2") == A1))
pt_whatever = csv.filter((pl.col("Pop1") == D2) | (pl.col("Pop2") == D2))
drossim_cnxj = csv.filter(((pl.col("Pop2") == A1) | (pl.col("Pop1") == A1)) & ((pl.col("Pop1") == B3) | (pl.col("Pop2") == B3)))["Estimate"].to_list()[0]
pt_cnxj = csv.filter(((pl.col("Pop2") == B3) & (pl.col("Pop1") == D2)) | ((pl.col("Pop1") == B3) & (pl.col("Pop2") == D2)))["Estimate"].to_list()[0]

# Update population list
poops += ["CnOther", "CnQTP", "CNXJ", "ISR"]
if D2 in poops:
	poops.remove(D2)

# Calculate f4 statistics
f4stats = []
for pop in poops:
    drossim_pop = drossim_whatever.filter((pl.col("Pop1") == pop) | (pl.col("Pop2") == pop))["Estimate"].to_list()[0]
    pt_pop = pt_whatever.filter((pl.col("Pop1") == pop) | (pl.col("Pop2") == pop))["Estimate"].to_list()[0]
    f4 = 0.5 * (drossim_pop + pt_cnxj - drossim_cnxj - pt_pop)
    f4stats.append(f4)

# Create DataFrame and write to CSV
dff = pl.DataFrame({"Pop": poops, "f4": f4stats})
dff.write_csv("f4ratios_dgn_cnother_nohet.csv")

