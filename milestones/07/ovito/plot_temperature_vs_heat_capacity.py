import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotx

# Prepare data
labels = [r"$E_\mathrm{total}$", r"Temperature}"]
headers = ["Etotal", "Temperature"]
mycolors = ['tab:blue', 'tab:green', 'tab:orange',
            'tab:brown', 'tab:grey', 'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

# Import Data
# numbers = ["55"]
# numbers = ["55", "147", "309", "561", "923", "1415", "2057"]
# numbers = ["2869", "3871", "5083", "6525", "8217", "10179", "12431", "14993", "17885", "21127", "24739", "28741"]
# numbers = ["55", "147", "309", "561", "923", "1415", "2057", "2869", "3871",
#            "5083", "6525", "8217", "10179", "12431", "14993", "17885",
#            "21127", "24739", "28741"]
numbers = ["1415", "2057", "2869", "3871",
           "5083", "6525", "8217", "10179", "12431", "14993", "17885",
           "21127", "24739", "28741"]

dfs = [pd.read_csv('energy_' + num + '.csv', header=None, usecols=[0, 1], sep="\t") for num in numbers]
for df in dfs:
    df.columns = headers

# Calculate Heat Capacities
dfs_heat_capacities = []
for i, df in enumerate(dfs):
    heat_capacities = []
    total_energy_diff = max(df[headers[0]]) - min(df[headers[0]])
    energy0, energy1 = 0, 0
    temperature0, temperature1 = 0, 0

    for p in np.linspace(0, len(df) - 1,  13, dtype=int):
        energy1, temperature1 = df[headers[0]].iloc[p], df[headers[1]].iloc[p]
        heat_capacity = (energy1 - energy0) / (temperature1 - temperature0)
        heat_capacities.append((temperature1, heat_capacity))
        energy0, temperature0 = energy1, temperature1
    dfs_heat_capacities.append(heat_capacities)

# Plot
plt.figure(figsize=(12, 7), dpi=100)
for i, df in enumerate(dfs):
    t, hc = zip(*dfs_heat_capacities[i])
    plt.plot(t, hc, color=mycolors[i % len(mycolors)], label=numbers[i])

# Decoration
plt.xlabel(r'Temperature in $K$')
plt.ylabel(r'Heat Capacity in $\frac{eV}{K}$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Temperature versus Heat Capacity", fontsize=22)
matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

path = r"T_vs_C"
plt.savefig(path + ".png")

plt.show()  # clears the plot
