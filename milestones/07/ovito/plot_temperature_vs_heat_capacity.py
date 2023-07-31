import pandas as pd
import matplotlib.pyplot as plt
import matplotx
import tikzplotlib

# Prepare data
labels = [r"$E_\mathrm{total}$", r"Temperature}"]
headers = ["Etotal", "Temperature", "Heat Capacity"]

# Import Data
# numbers = ["55"]
# numbers = ["55", "147", "309", "561", "923", "1415", "2057"]
# numbers = ["2869", "3871", "5083", "6525", "8217", "10179", "12431", "14993", "17885", "21127", "24739", "28741"]
numbers = ["55", "147", "309", "561", "923", "1415", "2057", "2869", "3871",
           "5083", "6525", "8217", "10179", "12431", "14993", "17885",
           "21127", "24739", "28741"]

dfs = [pd.read_csv('energy_' + num + '.csv', header=None, sep="\t") for num in numbers]
for df in dfs:
    df.columns = headers

# Set first element to 0
for df in dfs:
    df.iat[0, 2] = None

# Draw Plot
mycolors = ['tab:pink', 'tab:blue', 'tab:green', 'tab:orange',
            'tab:brown', 'tab:grey', 'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(12, 7), dpi=100)

for i, df in enumerate(dfs):
    plt.plot(df[headers[1]], df[headers[2]], "o", markersize=2, color=mycolors[i % len(mycolors)], label=numbers[i])
    # plt.text(df[headers[2]].iloc[-1], max(df[headers[2]]) + 50, numbers[i], horizontalalignment='center', color=mycolors[i % len(mycolors)])

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
# tikzplotlib.save(path + ".tex")
plt.savefig(path + ".png")

plt.show()  # clears the plot
