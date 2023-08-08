import pandas as pd
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib


# Prepare data
filenames = [r'energy_0.001000', r'energy_0.005000', r'energy_0.010000', r'energy_0.020000', r'energy_0.025000']
headers = ["Epot", "Ekin", "Etotal", "Temperature"]

# Import Data
dfs = [pd.read_csv(f'ovito/{filename}.csv', header=None, sep="\t") for filename in filenames]
for df in dfs:
    df.columns = headers

# Draw Plot
mycolors = ['tab:pink', 'tab:blue', 'tab:green', 'tab:orange',
            'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(12, 7), dpi=100)

plottable_tuples = []
for i, df in enumerate(dfs):
    plottable_tuples.append([])
    for time in range(len(dfs[0])):
        # print(len(df), time, time * len(df) // len(dfs[0]))
        e = df["Etotal"].iloc[time * len(df) // len(dfs[0])]
        plottable_tuples[i].append((time, e))

for i, points in enumerate(plottable_tuples):
    e, time = zip(*points)
    timestep = float(filenames[i][len("energy_"):])
    print(i, timestep)
    plt.plot(e, time, color=mycolors[i], label=str(timestep), linewidth=1)

# Decoration
# plt.ylim(-200, 60)
# plt.xlim(0, 100)
plt.xlabel(r'Time in $\sqrt{m \sigma^2 / \varepsilon}$')
plt.ylabel(r'Energy in $eV$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Total Energy for Different Timesteps", fontsize=22)
matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

tikzplotlib.save(r"ovito/different_timesteps.tex")
plt.savefig(r"ovito/different_timesteps.png")

plt.show()  # clears the plot
