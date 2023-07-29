import pandas as pd
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib

# Prepare data
labels = [r"$E_\mathrm{total}$", r"Temperature}"]
headers = ["Etotal", "Temperature"]

# Import Data
# numbers = ["923"]
numbers = ["55", "147", "309", "561", "923", "1415", "2057", "2869"]
# numbers = ["3871", "5083", "6525", "8217", "10179", "12431", "14993", "17885", "21127", "24739", "28741"]
# numbers = ["55", "147", "309", "561", "923", "1415", "2057", "2869", "3871",
#            "5083", "6525", "8217", "10179", "12431", "14993", "17885",
#            "21127", "24739", "28741"]

dfs = [pd.read_csv('energy_' + num + '.csv', header=None, usecols=[0, 1], sep="\t") for num in numbers]
for df in dfs:
    df.columns = headers

# Draw Plot
mycolors = ['tab:pink', 'tab:blue', 'tab:green', 'tab:orange',
            'tab:brown', 'tab:grey', 'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(8, 5), dpi=80)

for i, df in enumerate(dfs):
    plt.plot(df[headers[0]], df[headers[1]], color=mycolors[i % len(mycolors)], label=numbers[i])
    plt.text(df[headers[0]].iloc[-1], max(df[headers[1]]) + 20, numbers[i], horizontalalignment='center', color=mycolors[i % len(mycolors)])

# Decoration
plt.ylabel(r'Temperature in $K$')
plt.xlabel(r'Energy in $eV$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Total Energy versus Temperature", fontsize=22)
# matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

path = r"E_vs_T"
# tikzplotlib.save(path + ".tex")
plt.savefig(path + ".png")

plt.show()  # clears the plot
