import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotx

# Prepare data
labels = [r"$E_\mathrm{total}$", r"Temperature}"]
headers = ["Etotal", "Temperature"]

# numbers = ["8217"]
# numbers = ["55", "147", "309", "561", "923", "1415", "2057", "2869", "3871", "5083", "6525", "8217"]
numbers = ["55", "147", "309", "561", "923", "1415", "2057", "2869", "3871",
           "5083", "6525", "8217", "10179", "12431", "14993", "17885", "21127",
           "24739", "28741"]

# Import Data
dfs = [pd.read_csv('energy_' + num + '.csv', header=None, usecols=[0, 1], sep="\t") for num in numbers]
for df in dfs:
    df.columns = headers

# Draw Plot
mycolors = ['tab:green', 'tab:orange', 'tab:pink', 'tab:blue',
            'tab:brown', 'tab:grey', 'tab:red', 'tab:olive',
            'deeppink', 'steelblue', 'firebrick', 'mediumseagreen']

plt.figure(figsize=(8, 5), dpi=80)

# Clear intersections and latened heat files
with open("intersection points", "w") as f:
    f.write("")
with open("latened heat", "w") as f:
    f.write("")

for i, df in enumerate(dfs):
    upper_part = 3
    if i > 6:
        upper_part = 4

    plt.plot(df[headers[0]], df[headers[1]], color=mycolors[i % len(mycolors)], label=numbers[i])

    # Get lines for linear regression
    b, a = np.polyfit(df[headers[0]].array[:len(df)//upper_part], df[headers[1]].array[:len(df)//upper_part], deg=1)
    xseq = np.linspace(0, max(df[headers[0]].array), num=len(df))
    plt.plot(xseq, a + b * xseq, color=mycolors[i % len(mycolors)], lw=0.5)
    b2, a2 = np.polyfit(df[headers[0]].array[2*len(df)//3:], df[headers[1]].array[2*len(df)//3:], deg=1)
    plt.plot(xseq, a2 + b2 * xseq, color=mycolors[i % len(mycolors)], lw=0.5)

    # Draw bisecting line
    plt.plot(xseq, (a2 - (a2-a)/2) + (b2 - (b2-b)/2) * xseq, color=mycolors[i % len(mycolors)], lw=0.5)

    # Calculate intersection point
    intersections = []
    t = np.linspace(0, 50, 1000)
    prev_dif = 0
    t0, prev_ydata, prev_yline = None, None, None
    for j, ZIPPEDY in enumerate(zip(t, df[headers[1]].array, (a2 - (a2-a)/2) + (b2 - (b2-b)/2) * xseq)):
        t1, ydata, yline = ZIPPEDY
        new_dif = ydata - yline
        if np.abs(new_dif) < 1e-12:
            intersections.append((t1, ydata))
        elif new_dif * prev_dif < 0:  # interpolate intersection point
            denom = prev_dif - new_dif
            p = (df[headers[0]].iloc[j], -(ydata * prev_yline - yline * prev_ydata) / denom)
            if p[1] > 500:
                intersections.append(p)
        t0, prev_ydata, prev_yline, prev_dif = t1, ydata, yline, new_dif

    # Save intersection points
    with open("intersection points", "a") as f:
        f.write(f"{numbers[i]}\t" + "\t".join([str(p) for p in intersections]) + "\n")

    # Calculate vertical lines on intersection point from left line to right line
    py = intersections[len(intersections)//2][1]
    vertical_line = np.array([py for _ in range(len(df))])
    left_idx = np.argwhere(np.diff(np.sign(a2 + b2 * xseq - vertical_line))).flatten()[0]
    right_idx = np.argwhere(np.diff(np.sign(a + b * xseq - vertical_line))).flatten()[0]

    # Save latened heat
    with open("latened heat", "a") as f:
        f.write(f"{numbers[i]}\t{xseq[left_idx] - xseq[right_idx]}\n")

    # Draw latent heat line and melting point
    plt.plot(np.linspace(xseq[left_idx], xseq[right_idx], num=100), np.linspace((a2 + b2 * xseq)[left_idx], (a + b * xseq)[right_idx], num=100), color="red")
    plt.plot(*zip(*intersections), 'o', ms=7, color="black")

    # Draw cluster size
    plt.text(df[headers[0]].iloc[-1], max(df[headers[1]]) + 20, numbers[i], horizontalalignment='center', color=mycolors[i % len(mycolors)])

# Decoration
plt.ylabel(r'Temperature in $K$')
plt.xlabel(r'Energy in $eV$')
plt.yticks(fontsize=12, alpha=.7)
plt.title("Derivation of Latent Heat and Melting Point", fontsize=22)
# matplotx.line_labels()

# Plot grid behind graphs
plt.grid(axis='y', alpha=.3)

# Remove borders
plt.gca().spines["top"].set_alpha(0.0)
plt.gca().spines["bottom"].set_alpha(0.5)
plt.gca().spines["right"].set_alpha(0.0)
plt.gca().spines["left"].set_alpha(0.5)

path = r"LH_vs_N"
plt.savefig(path + ".png")

plt.show()  # clears the plot
