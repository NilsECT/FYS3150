import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

sns.set_context({'axes.linewidth': 0.9,
 'grid.linewidth': 0.6,
 'lines.linewidth': 1.5,
 'lines.markersize': 9.0,
 'patch.linewidth': 1.5,
 'xtick.major.width': 1.875,
 'ytick.major.width': 1.875,
 'xtick.minor.width': 1.5,
 'ytick.minor.width': 1.5,
 'xtick.major.size': 9.0,
 'ytick.major.size': 9.0,
 'xtick.minor.size': 6.0,
 'ytick.minor.size': 6.0,
 'font.size': 25.0,
 'axes.labelsize': 25.0,
 'axes.titlesize': 25.0,
 'xtick.labelsize': 25.0,
 'ytick.labelsize': 25.0,
 'legend.fontsize': 25.0,
'legend.title_fontsize': 18.0})

sns.set_style({'axes.facecolor': 'white',
 'axes.edgecolor': 'black',
 'axes.grid': True,
 'axes.axisbelow': 'line',
 'axes.labelcolor': 'black',
 'figure.facecolor': 'white',
 'grid.color': '#b9f2f0',
 'grid.linestyle': '-',
 'text.color': 'black',
 'xtick.color': 'black',
 'ytick.color': 'black',
 'xtick.direction': 'out',
 'ytick.direction': 'out',
 'patch.edgecolor': 'black',
 'patch.force_edgecolor': False,
 'image.cmap': 'viridis',
 'font.family': ['sans-serif'],
 'font.sans-serif':
    ['DejaVu Sans',
    'Bitstream Vera Sans',
    'Computer Modern Sans Serif',
    'Lucida Grande',
    'Verdana',
    'Geneva',
    'Lucid',
    'Arial',
    'Helvetica',
    'Avant Garde',
    'sans-serif'],
 'xtick.bottom': True,
 'xtick.top': False,
 'ytick.left': True,
 'ytick.right': False,
 'axes.spines.left': True,
 'axes.spines.bottom': True,
 'axes.spines.right': True,
 'axes.spines.top': True})

palette = sns.color_palette("colorblind")
sns.set_palette(palette=palette, desat=0.9, color_codes=True)

cm = 1/2.54
figsize=(12*cm, 10*cm)

def get_data(filename, delimiter=", ", col_names_ind=0):
    # henter ut data fra fil og legger det i en pandas DataFrame

    name = filename + ".txt"

    data = pd.read_csv(name, sep=delimiter, header=col_names_ind, index_col=False, engine='python')
    
    return data

def load_multiple(filenames):
    temp = []
    for i in filenames:
        temp.append(get_data(i))
    
    return pd.concat(temp, ignore_index=True)


def sortdata(data, col_name):
    '''sorts the dataframe data by column names, can be a list, the first will have biggest importance.
    Example sort with ["other", "age"]:
        name  age  other
        0  james   13      3
        1  jimmy   98      5
        2    mad   43      4
        3   blog   32      3
        4   pray   92      7
        5  boing   31      3
        6  extra   55      5
        7  annen   55      3
        Sorting
            name  age  other
        0  james   13      3
        1  boing   31      3
        2   blog   32      3
        3  annen   55      3
        4    mad   43      4
        5  extra   55      5
        6  jimmy   98      5
        7   pray   92      7

    '''
    
    return data.sort_values(by=col_name, ignore_index=True)


# eps_dist = get_data("epsilon_distribution")

# sns.displot(data=eps_dist, x="Energy [J]", col="Temperature [J/kb]", stat="probability", bins=100, kde=True)
# plt.show()

# ordered_1 = get_data("varying_cycles_ordered_1.000000")
# ordered_2 = get_data("varying_cycles_ordered_2.400000")
# unordered_1 = get_data("varying_cycles_unordered_1.000000")
# unordered_2 = get_data("varying_cycles_unordered_2.400000")

var_N = get_data("varying_cycles_1")

sns.lineplot(data=var_N.loc[var_N["Order"] == "Ordered"], x="MC cycles", y="Energy [J]", hue="Temperature [J/kb]")
plt.xscale("log")
plt.show()

sns.lineplot(data=var_N.loc[var_N["Order"] == "Unordered"], x="MC cycles", y="Energy [J]", hue="Temperature [J/kb]")
plt.xscale("log")
plt.show()