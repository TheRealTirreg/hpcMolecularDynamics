import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotx

mycolors = ['tab:green', 'tab:orange', 'tab:pink', 'tab:blue',
            'tab:brown', 'tab:grey', 'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']
plt.figure(figsize=(8, 5), dpi=80)

# Import and draw data
nums_and_points = []
with open('intersection points all', "r") as f:
    for i, line in enumerate(f):
        num, *pts = line.split("\t")
        for p in pts:
            x, y = p.strip("()\n").split(", ")
            plt.plot(float(x), float(y), "o", color=mycolors[i % len(mycolors)])
            plt.text(float(x), float(y) + 20, num, horizontalalignment='center', color=mycolors[i % len(mycolors)])


# Decoration
plt.xlabel(r'Energy in $eV$')
plt.ylabel(r'Temperature in $K$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Melting Points", fontsize=22)

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

path = r"Melting_points"
plt.savefig(path + ".png")

plt.show()  # clears the plot
