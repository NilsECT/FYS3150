import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

sns.set(rc={'figure.figsize':(8,7)})

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
 'font.size': 16.,
 'axes.labelsize': 16.,
 'axes.titlesize': 16.,
 'xtick.labelsize': 16.,
 'ytick.labelsize': 16.,
 'legend.fontsize': 16.,
'legend.title_fontsize': 16.0})

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

cycles = "Cycle"
temperature = "Temperature [J/kb]"
lattice = "Lattice"
energy = "Energy [J]"
sample = "Sample"
magnet = "Magnetisation"
burn = "Burn in [cycles]"
cv = "Heat capacity"
chi = "Susceptibility"


###########  LOOKING AT BURN IN ORDERED VS UNORDERED #######

var_N = get_data("varying_cycles_no_burn")

lattices = [20]    # [20, 40, 60, 100]

for lat in lattices:

    # # using fontsizes 16

    # # plt.figure(figsize=(12, 10))
    # ordered = var_N.loc[var_N["Order"] == "Ordered"]
    # sns.lineplot(data=ordered.loc[ordered[lattice] == lat], x=cycles, y=energy, hue=temperature)
    # plt.xscale("log")
    # plt.savefig("ordered_var_MC_%d.pdf" %lat)
    # plt.show()

    # # plt.figure(figsize=(12, 10))
    # unordered = var_N.loc[var_N["Order"] == "Unordered"]
    # sns.lineplot(data=unordered.loc[unordered[lattice] == lat], x=cycles, y=energy, hue=temperature)
    # plt.xscale("log")
    # plt.savefig("unordered_var_MC_%d.pdf" %lat)
    # plt.show()

    # plt.figure(figsize=(12, 10))
    sns.lineplot(data=var_N, x=cycles, y=energy, hue=temperature, style="Order")
    plt.xscale("log")
    plt.savefig("var_MC_%d.pdf" %lat)
    plt.show()

###########################################################

################### LOOKING AT ENERGY DISTRINUTION #######

# eps_dist = get_data("epsilon_distribution")

# using fontsizes 12

# temp = [1., 1.5, 2., 2.4]

# for i in temp:
#     data = eps_dist.loc[eps_dist[temperature] == i]
#     for lat in lattices:

#         # plt.figure(figsize=(15, 13))
#         sns.displot(data=data.loc[data[lattice] == lat], x="Energy [J]", stat="probability", bins=75, kde=True)
#         # title_temp_lattice_#MC_#samples
#         plt.savefig("energy_dist_%.d_%d.pdf" %(i*10, lat))
#         plt.show()

#######################################################

############# LOOKING AT PHASE TRANSITIONS ###############

# phase = get_data("phase_transition_varL")

# sns.lineplot(data=phase, x=temperature, y=cv)
# plt.savefig("cv_20.pdf")
# plt.show()

# sns.lineplot(data=phase, x=temperature, y=chi)
# plt.savefig("chi_20.pdf")
# plt.show()