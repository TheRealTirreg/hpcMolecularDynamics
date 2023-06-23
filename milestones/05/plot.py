import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib

# Prepare data
labels = ["With Neighborlist", "No Neighborlist"]
headers = ["n", "With Neighborlist", "No Neighborlist"]

# Import Data
df = pd.read_csv('ovito/perf.csv', header=None, sep=",")
df.columns = headers
print(df)

# Draw Plot
mycolors = ['tab:pink', 'tab:blue']
# ['tab:orange', 'tab:brown', 'tab:grey', 'tab:red',
# 'tab:olive', 'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(8, 5), dpi=80)

for i in range(2):
    print(df.shape, "\n")
    plt.plot(df[headers[0]], df[headers[i+1]], color=mycolors[i], label=labels[i])

# Decoration
plt.ylim(0, 1000)
plt.xlim(0, 1000)
plt.xticks(df[headers[0]])
plt.xticks(fontsize=12, alpha=.7)
plt.xlabel(r'Number of Atoms')
plt.ylabel(r'Execution time in $s$')
plt.title("Execution time MS05", fontsize=22)
plt.legend(loc="upper left")
# matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

# tikzplotlib.save(r"ovito/plot.tex")
plt.savefig(r"ovito/plot.png")

plt.show()  # clears the plot
