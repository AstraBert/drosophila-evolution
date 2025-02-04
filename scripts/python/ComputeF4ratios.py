import polars as pl

# Read CSV files
csv = pl.read_csv("F2stats_noinv.csv")
csv = csv.rename({"": "Pools"})
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
C4 = "TR_Ank_Yes_1_2021-10-16"
expression = f"0.5*({A1},{C4} + {D2},{B3} - {A1},{B3} - {D2},{C4})"

# Filter data based on conditions
drossim_trk = csv.filter(((pl.col("Pop1") == A1) & (pl.col("Pop2") == C4)) | ((pl.col("Pop2") == A1) & (pl.col("Pop1") == C4)))["Estimate"].to_list()[0]
pt_trk = csv.filter(((pl.col("Pop1") == D2) & (pl.col("Pop2") == C4)) | ((pl.col("Pop2") == D2) & (pl.col("Pop1") == C4)))["Estimate"].to_list()[0]
drossim_cnxj = csv.filter(((pl.col("Pop2") == A1) & (pl.col("Pop1") == B3)) | ((pl.col("Pop1") == A1) & (pl.col("Pop2") == B3)))["Estimate"].to_list()[0]
pt_cnxj = csv.filter(((pl.col("Pop2") == B3) & (pl.col("Pop1") == D2)) | ((pl.col("Pop1") == B3) & (pl.col("Pop2") == D2)))["Estimate"].to_list()[0]

f4 = 0.5*(drossim_trk + pt_cnxj - drossim_cnxj - pt_trk)
print(f"F4({A1},{D2};{B3},{C4}):", f4)

f4_csv = pl.read_csv("f4ratios_dgn_cnother.csv")
other_f4s = f4_csv["f4"].to_list()
f4_ratios=[]
for el in other_f4s:
	ratio = abs(el/f4)
	f4_ratios.append(abs(1-ratio))
pseudo_df = {"Pop": f4_csv["Pop"].to_list(),"f4ratio":f4_ratios}
df = pl.DataFrame(pseudo_df)
df.write_csv("f4ratios_dgn_cnother_trk.csv")
