import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from mpl_toolkits.basemap import Basemap  # For basemap functionality
from matplotlib.colors import Normalize
from matplotlib import cm

p1 = "DrosSim"
p2 = "DGN"
p3 = "CYP_1"
continent = "AS"

# Read your data
df = pd.read_csv("../data/all_f4stats_wsim.tsv", sep="\t")
df_pools = pd.read_csv("../results/pools.csv")
pops = df_pools["NAME"].to_list() 
conts = df_pools["CONTINENT"].to_list() 
pops2cont = {pops[i]: conts[i] for i in range(len(pops))}
df_subs = df[(df["PopO"] == p2) & (df["PopX"] == p1) & (df["PopA"] == p3) & (df["PopC"].map(pops2cont) == continent)]

# Calculate coordinates and prepare the data
longs = df_pools['LONG'] 
lats = df_pools['LAT']  
pops2dist = {pops[i]: [longs[i],lats[i]] for i in range(len(pops))}
focals = df_subs["PopC"].to_list() 
f3s = df_subs["Estimate"].to_list() 
focal2zscoref3 = {focals[i]: [f3s[i]] for i in range(len(focals))}


latsused = [] 
longsused = [] 

# Prepare data for Kriging
for k in focal2zscoref3:
    focal2zscoref3[k].insert(0, pops2dist[k][0])
    focal2zscoref3[k].insert(1, pops2dist[k][1])
    longsused.append(pops2dist[k][0])
    latsused.append(pops2dist[k][1])

# Define plot boundaries
MINY = min(latsused) - 10
MAXY = max(latsused) + 10
MINX = min(longsused) - 10
MAXX = max(longsused) + 10

data = np.array([focal2zscoref3[k] for k in focal2zscoref3])
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
im = m.imshow(z, extent=[MINX, MAXX, MINY, MAXY], origin='lower', cmap="RdBu_r", alpha=0.7, norm=norm)

# Add colorbar for the heatmap
plt.colorbar(im, label='F4 Value')

# Scatter plot the original data points on the map
for k in focal2zscoref3:
    xpt, ypt = m(focal2zscoref3[k][0], focal2zscoref3[k][1])
    plt.scatter(xpt, ypt, c='yellow', marker='x', s=100)
    # Annotate with population name and F3 value
    plt.annotate(
        f'{k}\nF4: {focal2zscoref3[k][2]:.4f}',
        (xpt, ypt),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=8,
        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2')
    )

# Title and labels
plt.title(f'Outgroup: {p1}, PopC: {p2}, PopA: {p3}')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save and show the plot
plt.tight_layout()
plt.savefig(f"./F4_{p1}_{p2}_{p3}_{continent}.png", dpi=300, transparent=True)  # Save with transparent background
plt.show()
