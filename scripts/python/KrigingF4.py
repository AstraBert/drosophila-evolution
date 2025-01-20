import pandas as pd
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from mpl_toolkits.basemap import Basemap  # For basemap functionality
from matplotlib.colors import Normalize
from matplotlib import cm


continent = "EU"

df_pools = pd.read_csv("data/f4_stats_all/dest_drosevol_latlong.csv")
df_pools = df_pools[(df_pools["Continent"] == "EU") & (df_pools["sampleId"] != "FI_Pir_Aka_1_2021-09-17") & (df_pools["sampleId"] != "CNXJ")]
pops = df_pools["sampleId"].to_list()
lats = df_pools["lat"].to_list()
longs = df_pools["long"].to_list()
df_stats = pl.read_csv("data/f4_stats_all/f4_cnxj_fin.csv")
df_stats = df_stats.filter(pl.col("Pop").is_in(pops))
pops1 = df_stats["Pop"].to_list()
stats = df_stats["f4"].to_list()
pops2stats = {pops1[i]: [stats[i]] for i in range(len(pops1))}
pops2coord = {pops[i]: [lats[i], longs[i]] for i in range(len(pops))}
pops2every = {}
for k in pops2coord:
    l1 = pops2coord[k]
    l2 = pops2stats[k]
    pops2every.update({k: l1+l2})

# Define plot boundaries
MINY = min(lats) - 10
MAXY = max(lats) + 10
MINX = min(longs) - 10
MAXX = max(longs) + 10

data = np.array([pops2every[k] for k in pops2every])
print(data)
gridx = np.arange(MINX, MAXX, 0.1)
gridy = np.arange(MINY, MAXY, 0.1)

# Perform Ordinary Kriging
OK = OrdinaryKriging(
    data[:, 0],
    data[:, 1],
    data[:, 2],
    variogram_model="linear",
    verbose=False,
    enable_plotting=False,
)

# Execute Kriging
z, ss = OK.execute("grid", gridx, gridy)

# Write ASC grid
kt.write_asc_grid(gridx, gridy, z, filename="output.asc")

# Create the map with Basemap
plt.figure(figsize=(12, 10))

# Create Basemap instance for Europe, Africa, and Asia
m = Basemap(projection='merc', llcrnrlat=MINY, urcrnrlat=MAXY, llcrnrlon=MINX, urcrnrlon=MAXX, resolution='i')

# Draw coastlines, countries, and other map elements (without background color)
m.drawcoastlines()
m.drawcountries()
m.drawmapboundary(fill_color='none')  # Make the map boundary transparent
m.drawrivers()
m.drawstates()

# Convert grid points to map projection coordinates
xx, yy = np.meshgrid(gridx, gridy)
x, y = m(xx, yy)

vmin, vmax = min(data[:,2]), max(data[:,2])
abs_max = max(abs(vmin), abs(vmax))
norm = Normalize(vmin=-abs_max, vmax=abs_max)
# Overlay the heatmap onto the map (make sure the heatmap is not covered)
im = m.imshow(z, extent=[MINX, MAXX, MINY, MAXY], origin='lower', cmap="RdBu_r", alpha=0.7)

# Add colorbar for the heatmap
plt.colorbar(im, label='F4 Value')

# Scatter plot the original data points on the map
for k in pops2every:
    xpt, ypt = m(pops2every[k][1], pops2every[k][0])
    plt.scatter(xpt, ypt, c='orange', marker='x', s=100)
    # Annotate with population name and F3 value
    # plt.annotate(
    #     f'{k}\nF4: {pops2every[k][2]:.4f}',
    #     (xpt, ypt),
    #     xytext=(5, 5),
    #     textcoords='offset points',
    #     fontsize=8,
    #     bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7),
    #     arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2')
    # )

# Title and labels
plt.title(f'Outgroup: D. simulans, PopC: DGN - Zambia, PopA: China - Other')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save and show the plot
plt.tight_layout()
plt.savefig(f"imgs/F4_wonorm_cnxjfin_EU.png", dpi=300, transparent=True)  # Save with transparent background
plt.show()
