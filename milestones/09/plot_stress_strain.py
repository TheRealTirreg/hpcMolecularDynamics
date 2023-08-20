import pandas as pd
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib

# Prepare data
labels = [r"$E_\mathrm{total}$", r"Temperature}"]
headers = ["Etotal", "Temperature", "xStress", "yStress", "zStress", "zDomainSize"]

# Import Data
df = pd.read_csv(r'ovito/energy_whisker_small.csv', header=None, sep="\t")
df.columns = headers

# Draw Plot
mycolors = ['tab:pink', 'tab:blue', 'tab:green', 'tab:orange',
            'tab:brown', 'tab:grey', 'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(12, 7), dpi=100)

plt.plot(df["zDomainSize"], df["zStress"], color=mycolors[0])

# Decoration
plt.xlabel(r'Domain length in $z$-direction in $\AA$', fontsize=16)
plt.ylabel(r'Stress in $\frac{eV}{\AA^2}$', fontsize=14)
plt.yticks(fontsize=12, alpha=.7)
plt.title("Stress versus Strain", fontsize=22)
# matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

path = r"ovito/Stress_vs_Strain"
plt.savefig(path + ".png")

plt.show()  # clears the plot
