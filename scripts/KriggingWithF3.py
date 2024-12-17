import pandas as pd  
import numpy as np
import matplotlib.pyplot as plt
import pykrige.kriging_tools as kt
from pykrige.ok import OrdinaryKriging


df = pd.read_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/results/f3stats.tsv", sep="\t")

df_subs = df[(df["Pop1"] == "DGN") & (df["Pop2"] == "CnOther")]

# Haversine function to calculate distance
def haversine(lat1, lon1, lat2, lon2):
    # Radius of Earth in kilometers
    R = 6371.0
    
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    
    # Differences in coordinates
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    
    # Haversine formula
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    # Distance in kilometers
    distance = R * c
    return distance


# Define a reference point (e.g., the equator or a specific location)
# Example: Coordinates of the equator (latitude = 0, longitude = 0)
ref_lat = 0.0
ref_lon = 0.0

df_pools = pd.read_csv("/gatk_modified/userdata/abertelli/drosophila-evolution/results/pools.csv")
# Calculate distance from each point to the reference point
df_pools['distance_to_ref'] = df_pools.apply(lambda row: haversine(row['LAT'], row['LONG'], ref_lat, ref_lon), axis=1)
pops = df_pools["NAME"].to_list() 
dists = df_pools['distance_to_ref'].to_list()
# Output the DataFrame with distances
pops2dist = {pops[i]: dists[i] for i in range(len(pops))}
focals = df_subs["Focal"].to_list() 
zscores = df_subs["Z-score"].to_list() 
f3s = df_subs["Estimate"].to_list() 
focal2zscoref3 = {focals[i]: [zscores[i], f3s[i]] for i in range(len(focals))}

MINY= min(zscores) - 1
MAXY = max(zscores) + 1
MINX = min(dists) - 1000
MAXX = max(dists) + 1000

for k in focal2zscoref3:
    focal2zscoref3[k].insert(0, pops2dist[k]) 


data = np.array([focal2zscoref3[k] for k in focal2zscoref3])

gridx = np.arange(MINX, MAXX, 0.1)
gridy = np.arange(MINY, MAXY, 0.1)

# # Perform Ordinary Kriging
OK = OrdinaryKriging(
    data[:, 0],
    data[:, 1],
    data[:, 2],
    variogram_model="linear",
    verbose=False,
    enable_plotting=False,
)

# # Execute Kriging
z, ss = OK.execute("grid", gridx, gridy)

# # Write ASC grid
kt.write_asc_grid(gridx, gridy, z, filename="output.asc")

# # Create a more detailed visualization
plt.figure(figsize=(12, 10))

# # Create the main heatmap with improved color mapping
im = plt.imshow(
    z,
    extent=[MINX, MAXX, MINY, MAXY],  # Set correct geographical extent
    origin='lower',  # Ensure correct orientation
    cmap='viridis',  # Choose a perceptually uniform colormap
    aspect='auto'    # Adjust aspect ratio
)

# Add a colorbar
plt.colorbar(im, label='F3 Value')

# # Add labels and title
plt.xlabel('Distance from the Equator (km)')
plt.ylabel('Z-score')
plt.title('P1: DGN, P2: CnOther')

# # Scatter plot of original data points with annotations
for k in focal2zscoref3:
    plt.scatter(
        focal2zscoref3[k][0],  # Longitude
        focal2zscoref3[k][1],  # Latitude
        c='red',
        marker='x',
        s=100  # Increased marker size
    )
    # Annotate with population name and F3 value
    plt.annotate(
        f'{k}\nF3: {focal2zscoref3[k][2]:.4f}',
        (focal2zscoref3[k][0], focal2zscoref3[k][1]),
        xytext=(5, 5),  # 5 points offset
        textcoords='offset points',
        fontsize=8,
        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.7),
        arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.2')
    )

plt.legend(['Data Points'], loc='best')

# # Save and show the plot
plt.tight_layout()
plt.savefig("./results/F3_all_DGN_CnOther.png", dpi=300)
plt.show()
