import pandas as pd
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import Normalize
from matplotlib import cm

# Load data
df_pools = pd.read_csv("data/f4_stats_all/dest_drosevol_latlong.csv")
df_pools = df_pools[(df_pools["Continent"] == "EU") & (df_pools["sampleId"] != "ISR") & (df_pools["sampleId"] != "DGN")]
pops = df_pools["sampleId"].to_list()
lats = df_pools["lat"].to_list()
longs = df_pools["long"].to_list()
df_stats = pl.read_csv("data/f4_stats_all/f4_dgn_isr.csv")
df_stats = df_stats.filter(pl.col("Pop").is_in(pops))
pops1 = df_stats["Pop"].to_list()
print(pops1)
stats = df_stats["f4"].to_list()

# Create dictionaries for mapping
pops2stats = {pops1[i]: [stats[i]] for i in range(len(pops1))}
pops2coord = {pops1[i]: [lats[i], longs[i]] for i in range(len(pops1))}
pops2every = {}
for k in pops2coord:
    l1 = pops2coord[k]
    l2 = pops2stats[k]
    pops2every.update({k: l1 + l2})
print(pops2every)
# Define plot boundaries
MINY = min(lats) - 10
MAXY = max(lats) + 10
MINX = min(longs) - 10
MAXX = max(longs) + 10

data = np.array([pops2every[k] for k in pops2every])
print(data)


# Create the map with Basemap
plt.figure(figsize=(12, 10))
ax = plt.gca()  # Get the current axes
m = Basemap(ax=ax, projection='merc', llcrnrlat=MINY, urcrnrlat=MAXY, llcrnrlon=MINX, urcrnrlon=MAXX, resolution='i')
m.drawcoastlines()
m.drawmapboundary(fill_color='none')


# Scatter plot of the data points
vmin, vmax = min(data[:, 2]), max(data[:, 2])
abs_max = max(abs(vmin), abs(vmax))
norm = Normalize(vmin=-abs_max, vmax=abs_max)


# Prepare scatter plot data
x_coords = []
y_coords = []
f4_values = []
for k in pops2every:
    xpt, ypt = m(pops2every[k][1], pops2every[k][0])
    x_coords.append(xpt)
    y_coords.append(ypt)
    f4_values.append(pops2every[k][2]) # Store the f4 value to use in color mapping

# Create scalar mappable
sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)


# Create scatter plot, map the values to colors, and set the colormap
scatter = plt.scatter(x_coords, y_coords, c=f4_values, cmap=sm.cmap, marker='o', s=200, alpha=0.7, norm=sm.norm)

# Add colorbar for the F4 values
plt.colorbar(scatter, label='F4 Value')


# Title and labels
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save and show the plot
plt.tight_layout()
plt.savefig(f"imgs/F4_scatterplot_dgnisr.png", dpi=300, transparent=True)  # Save with transparent background
plt.show()
