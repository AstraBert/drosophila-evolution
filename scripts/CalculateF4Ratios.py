import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging
from mpl_toolkits.basemap import Basemap  # For basemap functionality

pO = "DrosSim"
pA = "EERU_1"
pC = "ISR"


df = pd.read_csv("F4stats.tsv", sep="\t")

df_sel = df[(df["PopC"] == pC) & (df["PopO"] == pO) & (df["PopA"] == pA)]
print(df_sel)
df_pools = pd.read_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/results/pools.csv")

Xs= df_sel["X"].to_list()
Bs= df_sel["A"].to_list()

longs = df_pools['LONG'] 
lats = df_pools['LAT']
pops = df_pools["NAME"].to_list() 
pops2dist = {pops[i]: [longs[i],lats[i]] for i in range(len(pops))}
f4s = df_sel["Estimate"].to_list() 
focal2f4 = {Bs[i]: [f4s[i]] for i in range(len(Bs))}

# Define plot boundaries
MINY = min(lats) - 10
MAXY = max(lats) + 10
MINX = min(longs) - 10
MAXX = max(longs) + 10

# Prepare data for Kriging
for k in focal2f4:
    focal2f4[k].insert(0, pops2dist[k][0])
    focal2f4[k].insert(1, pops2dist[k][1])

data = np.array([focal2f4[k] for k in focal2f4])
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

# Overlay the heatmap onto the map (make sure the heatmap is not covered)
im = m.imshow(z, extent=[MINX, MAXX, MINY, MAXY], origin='lower', cmap='coolwarm', alpha=0.7)

# Add colorbar for the heatmap
plt.colorbar(im, label='F3 Value')

# Scatter plot the original data points on the map
for k in focal2f4:
    xpt, ypt = m(focal2f4[k][0], focal2f4[k][1])
    plt.scatter(xpt, ypt, c='yellow', marker='x', s=100)
    # Annotate with population name and F3 value
    plt.annotate(
        f'{k}\nF3: {focal2f4[k][2]:.4f}',
        (xpt, ypt),
        xytext=(5, 5),
        textcoords='offset points',
        fontsize=8,
        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2')
    )

# Title and labels
plt.title(f'Outgroup: {pO}, InnerOutgroup: {pC}, InnerNode: {pA}')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Save and show the plot
plt.tight_layout()
plt.savefig(f"./F4_{pO}_{pC}_{pA}.png", dpi=300, transparent=True)  # Save with transparent background
plt.show()


