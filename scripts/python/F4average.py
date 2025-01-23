from statistics import mean, stdev
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use("seaborn-v0_8-paper")

groups = {"Western": ["ES", "FR", "GB", "GR", "DK", "CH", "NL", "SE"], "Western Border": ["DE", "IT"], "Eastern Border": ["HU", "PL", "RS", "AT"], "Eastern": ["UA", "BY", "RU", "FI"]}

df_stats = pd.read_csv("data/f4_stats_all/Dest_F4_stats.csv")
pops = df_stats["Pop"].to_list()
stats = df_stats["F4"].to_list()

data = {key: [] for key in groups}

for i in range(len(pops)):
    for k in groups:
        if pops[i].split("_")[0] in groups[k]:
            data[k].append(stats[i])
        else:
            continue

xlabs = [key for key in data]
yvals = [mean(data[key]) for key in data]
yerrs = [stdev(data[key]) for key in data]
colors = ["#5cd65c", "#85e085", "#ffb84d", "#ff9900"]

fig, ax = plt.subplots(figsize=(10,5))
bars = ax.bar(xlabs, yvals, yerr=yerrs, color=colors)
ax.set_xlabel("Cluster")
ax.set_ylabel("F4")
for bar in bars:
    height = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        height,
        f"{height:.5f}",
        ha="right",
        va="top",
    )

ax.set_title("Average F4 per European Clusters")
fig.savefig("imgs/F4_barplot.png")
fig.show()