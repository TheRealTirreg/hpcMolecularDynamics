import pandas as pd
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib

# Prepare data
labels = [r"$E_\mathrm{pot}$", r"$E_\mathrm{kin}$", r"$E_\mathrm{total}$"]
headers = ["Epot", "Ekin", "Etotal"]

# Import Data
df = pd.read_csv('ovito/energy.csv', header=None, sep="\t")
df.columns = headers

# Draw Plot
mycolors = ['tab:pink', 'tab:blue', 'tab:green']
# ['tab:orange', 'tab:brown', 'tab:grey', 'tab:red',
# 'tab:olive', 'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(8, 5), dpi=80)

for i in range(3):
    print(df.shape, "\n")
    plt.plot(df[headers[i]], color=mycolors[i], label=labels[i])

# Decoration
plt.ylim(-200, 60)
plt.xlim(0, 100)
plt.xlabel(r'Time in $\sqrt{m \sigma^2 / \varepsilon}$')
plt.ylabel(r'Energy in $eV$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Energy Conservation", fontsize=22)
matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

tikzplotlib.save(r"ovito/plot.tex")
plt.savefig(r"ovito/plot.png")

plt.show()  # clears the plot
