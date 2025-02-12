from statistics import median, stdev
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("seaborn-v0_8-paper")

groups = {
    "Western": ["ES", "FR", "GB", "GR", "DK", "CH", "NL", "SE"],
    "Western Border": ["DE", "IT"],
    "Eastern Border": ["HU", "PL", "RS", "AT"],
    "Eastern": ["UA", "BY", "RU", "FI"],
    "Middle East": ["TR", "ISR"]
}

df_stats = pd.read_csv("data/f4_stats_all/f4ratios_dgn_cnother_fin.csv")
pops = df_stats["Pop"].to_list()
stats = df_stats["f4ratio"].to_list()

data = {key: [] for key in groups}

for i in range(len(pops)):
    for k in groups:
        if pops[i].split("_")[0] in groups[k]:
            data[k].append(stats[i])

xlabs = list(data.keys())
yvals = list(data.values())

fig, ax = plt.subplots(figsize=(10, 5))

# Create a boxplot
box = ax.boxplot(yvals, patch_artist=True, labels=xlabs, showfliers=True)

# Customize colors
colors = ["#5cd65c", "#85e085", "#ffb84d", "#ff9900", "#e4080a"]
for patch, color in zip(box["boxes"], colors):
    patch.set_facecolor(color)

for med in box["medians"]:
	med.set(color="black", linewidth=1)

medians = [round(med, 5) for med in [median(data[key]) for key in data]]
for i, median_val in enumerate(medians):
	ax.text(i + 1, 0, f"{median_val:.5f}", ha="center", va="bottom", fontsize=8, fontweight="bold")

ax.set_ylabel("F4 Ratio")
ax.set_title("F4 Ratio Distribution per European Clusters (Finland)")

fig.savefig("imgs/f4ratios/F4_boxplot_dgncnotherfin.png")
fig.show()

