import pandas as pd
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib

# Prepare data
labels = [r"$E_\mathrm{total}$", r"Temperature}"]
headers = ["Epot", "Ekin", "Etotal", "Temperature"]

# Import Data
num = "147"
df = pd.read_csv('energy_' + num + '.csv', header=None, sep="\t")
df.columns = headers

# Draw Plot
mycolors = ['tab:pink', 'tab:blue', 'tab:green']
# ['tab:orange', 'tab:brown', 'tab:grey', 'tab:red',
# 'tab:olive', 'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(8, 5), dpi=80)

plt.plot(df[headers[2]], df[headers[3]], color=mycolors[1], label=labels[0])

# Decoration
plt.xlim(-550, -300)
plt.ylim(0, 4000)
plt.xlabel(r'Energy in $eV$')
plt.ylabel(r'Temperature in $K$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Total Energy versus Temperature " + num + " Atoms", fontsize=22)
matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

path = r"MP_vs_N"
tikzplotlib.save(path + ".tex")
plt.savefig(path + ".png")

plt.show()  # clears the plot
