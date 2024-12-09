import polars as pl  
import matplotlib.pyplot as plt
import numpy as np
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging

def most_fre_in_col(df: pl.DataFrame, col: str, col1: str):
    most_frequent = (
        df.group_by(col, col1)
        .agg(pl.count().alias("count"))   # Count occurrences of each locality
        .sort("count", descending=True)      # Sort by count in descending order
        .limit(10)                         # Take the first (most frequent) row
    )
    return most_frequent

# Read and process data
DEST_SAMPLES = pl.read_csv("data/dest_v2.samps_3May2024.csv", null_values=["NA"])
DEST_SUBSET = DEST_SAMPLES.select(["locality", "lat", "long", "country", "continent"])
DEST_EUROPE = DEST_SUBSET.filter((pl.col("continent")=="Europe"))
eu_pops = DEST_EUROPE["locality"].to_list() 
DEST_FILTERED = DEST_SUBSET.filter(pl.col("continent")=="Europe")
DEST_F3 = pl.read_csv("data/DEST_f3.tsv", separator="\t", n_threads=16)

# Process unique locations
unique_lat_long = DEST_FILTERED.select(["locality", "lat", "long"]).unique().drop_nulls()
pops = unique_lat_long["locality"].to_list()
lats = [round(lat,2) for lat in unique_lat_long["lat"].to_list()]
longs = [round(lon,2) for lon in unique_lat_long["long"].to_list()]

MAX_LAT = max(lats)
MIN_LAT = min(lats)
MAX_LONG = max(longs)
MIN_LONG = min(longs)

# Create population to lat-long mapping
pops2latlongs = {pops[k]: [lats[k], longs[k]]  for k in range(len(pops))}

# Filter F3 data of interest
DEST_F3_INTEREST = DEST_F3.filter(
    (pl.col("focal").is_in(pops)) &
    (pl.col("p1")=="ET_Oro_Fic") &
    (pl.col("p2")=="US_Flo_Hom")
)

# ET_Oro_Fic US_Flo_Hom


# f3_pops = DEST_F3_INTEREST["focal"].to_list()
# print(f3_pops)
# f3_f3s = DEST_F3_INTEREST["f3"].to_list()

# # Prepare data for Kriging
# data = np.array([[pops2latlongs[f3_pops[i]][0], pops2latlongs[f3_pops[i]][1], f3_f3s[i]] for i in range(len(f3_pops))])

# # Create grid
# gridx = np.arange(MIN_LAT, MAX_LAT, 1)
# gridy = np.arange(MIN_LONG, MAX_LONG, 1)

# # Perform Ordinary Kriging
# OK = OrdinaryKriging(
#     data[:, 0],
#     data[:, 1],
#     data[:, 2],
#     variogram_model="linear",
#     verbose=False,
#     enable_plotting=False,
# )

# # Execute Kriging
# z, ss = OK.execute("grid", gridx, gridy)

# # Write ASC grid
# kt.write_asc_grid(gridx, gridy, z, filename="output.asc")

# # Create a more detailed visualization
# plt.figure(figsize=(12, 10))

# # Create the main heatmap with improved color mapping
# im = plt.imshow(
#     z,
#     extent=[MIN_LONG, MAX_LONG, MIN_LAT, MAX_LAT],  # Set correct geographical extent
#     origin='lower',  # Ensure correct orientation
#     cmap='viridis',  # Choose a perceptually uniform colormap
#     aspect='auto'    # Adjust aspect ratio
# )

# # Add a colorbar
# plt.colorbar(im, label='F3 Value')

# # Add labels and title
# plt.xlabel('Longitude')
# plt.ylabel('Latitude')
# plt.title('F3 Kriging Interpolation')

# # Scatter plot of original data points with annotations
# for i in range(len(f3_pops)):
#     plt.scatter(
#         pops2latlongs[f3_pops[i]][1],  # Longitude
#         pops2latlongs[f3_pops[i]][0],  # Latitude
#         c='red',
#         marker='x',
#         s=100  # Increased marker size
#     )
#     # Annotate with population name and F3 value
#     plt.annotate(
#         f'{f3_pops[i]}\nF3: {f3_f3s[i]:.4f}',
#         (pops2latlongs[f3_pops[i]][1], pops2latlongs[f3_pops[i]][0]),
#         xytext=(5, 5),  # 5 points offset
#         textcoords='offset points',
#         fontsize=8,
#         bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7),
#         arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2')
#     )

# plt.legend(['Data Points'], loc='best')

# # Save and show the plot
# plt.tight_layout()
# plt.savefig("./results/F3_kriging_with_annotations.png", dpi=300)
# plt.show()
