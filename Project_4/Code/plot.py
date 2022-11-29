import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from scipy.stats import linregress

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

def load_multiple_sample(filenames, n_samples):
    temp = []
    for i in range(len(filenames)):
        data = get_data(filenames(i))
        data["Sample"] += n_samples*i
        temp.append(data)

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
lattice = "Lattice size"
energy = "Average energy [J]"
sample = "Sample"
magnet = "Average magnetisation"
burn = "Burn in [cycles]"
cv = "Heat capacity"
chi = "Susceptibility"
average_energy = "Average energy [J]"

###########  COMPARISON WITH ANALYTICAL RESULTS #######

data = get_data('analytical_comparison')

analytical_labels = ['Analytical energy per spin [J]', 'Analytical absolute magnetisation per spin', 'Analytical heat capacity', 'Analytical susceptibility']
labels = ['Energy per spin [J]', 'Absolute magnetisation per spin', cv, chi]

for i, axis_label, anal_axis_label in zip(range(4), labels, analytical_labels):
    plt.figure(figsize=(12, 7))

    graph = sns.lineplot(data=data, x="MC cycles",  y=axis_label, color='red', label='Computed')#, color=palette[0])
    graph.axhline(data[anal_axis_label][0], linestyle='--', label='Analytical', color='orange')

    plt.legend()
    plt.xscale("log")
    plt.savefig(f"analytical_{axis_label.split(' ')[0]}.pdf")
    plt.show()

###########  LOOKING AT BURN IN ORDERED VS UNORDERED #######

var_N = get_data("varying_cycles_no_burn_lattices_true")

lattice = 20


plt.figure(figsize=(12, 10))
sns.lineplot(data=var_N, x=cycles, y=average_energy, hue=temperature, style="Order")
plt.xscale("log")
plt.savefig("var_MC_e_%d_true.pdf" %lattice)

plt.figure(figsize=(12, 10))
sns.lineplot(data=var_N, x=cycles, y=magnet, hue=temperature, style="Order")
plt.xscale("log")
plt.savefig("var_MC_m_%d_true.pdf" %lattice)

###########################################################

################### LOOKING AT ENERGY DISTRIBUTION #######

eps_dist = get_data("epsilon_distribution")

temp = [1., 1.5, 2., 2.4]

for i in temp:
    data = eps_dist.loc[eps_dist[temperature] == i]
    

    plt.figure(figsize=(15, 13))
    sns.displot(data=data, x="Energy [J]", stat="probability", bins=75, kde=True)
    plt.savefig("energy_dist_%.d_%d.pdf" %(i*10, lattice))


#######################################################

############# LOOKING AT PHASE TRANSITIONS ###############

lattices = [20,40,60,80,100]
steps = [23,24,25,100]
T_c_max = np.zeros((len(steps),len(lattices)))
datas = []
for i,step in enumerate(steps):
    data = get_data(f"phase_transition_varL_step_{step}")
    datas.append(data)
    for j,lat in enumerate(lattices):
        idxmax = data.loc[data[lattice] == lat][chi].idxmax()
        T_c_max[i,j] = data[temperature].to_numpy()[idxmax]


T_c_mean = np.zeros(len(T_c_max[0]))
for i in range(len(T_c_mean)):
    T_c_mean[i] = np.mean(T_c_max[:,i])

L_frac = 1/lattices

lin = linregress(L_frac, T_c_mean)
a = lin.slope
T_c_infty = lin.intercept

da = lin.stderr
dT_c = lin.intercept_stderr
T_c_onsager = 2.269
print("For all results:")
print(f"T_c(L->infinity) = {T_c_infty:.4f} pm {dT_c:.4f}")
print(f"Absolut error: {abs(T_c_infty-T_c_onsager)}")
print(f"Relativ error: {abs(T_c_infty-T_c_onsager)/T_c_onsager}")
print()

plt.figure()
plt.plot(L_frac,T_c_mean,"o",label="Computed $T_c(L)$")
plt.plot(L_frac,a*L_frac+T_c_infty,label="Linear fit")
plt.xlabel("1/L")
plt.ylabel("T_c [J/k_B]")
plt.legend()
plt.savefig("linear_fit_T_c_1.pdf")

T_c_mean = np.zeros(len(T_c1))
for i in range(len(T_c1)):
    T_ci = np.array([T_c1[i], T_c3[i], T_c4[i]])
    T_c_mean[i] = np.mean(T_ci)

L_frac = 1/lattices

lin = linregress(L_frac, T_c_mean)
a = lin.slope
T_c_infty = lin.intercept

da = lin.stderr
dT_c = lin.intercept_stderr
T_c_onsager = 2.269
print("Neglecting Fridas results:")
print(f"T_c(L->infinity) = {T_c_infty:.4f} pm {dT_c:.4f}")
print(f"Absolute error: {abs(T_c_infty-T_c_onsager)}")
print(f"Relative error: {abs(T_c_infty-T_c_onsager)/T_c_onsager}")
plt.figure()
plt.plot(L_frac,T_c_mean,"o",label="Computed $T_c(L)$")
plt.plot(L_frac,a*L_frac+T_c_infty,label="Linear fit")
plt.xlabel("1/L")
plt.ylabel("T_c [J/k_B]")
plt.legend()
plt.savefig("linear_fit_T_c_2.pdf")

# Plotting av C_V og chi og sånt


phase_all = pd.concat(datas, ignore_index=True)

plt.figure()
sns.lineplot(data=phase_all, x=temperature, y=cv, hue=lattice)
plt.savefig("specific_heat_capacity.pdf")

plt.figure()
sns.lineplot(data=phase_all, x=temperature, y=chi, hue=lattice)
plt.savefig("susceptibility.pdf")

plt.figure()
sns.lineplot(data=phase_all, x=temperature, y=average_energy, hue=lattice)
plt.savefig("epsilon.pdf")

plt.figure()
sns.lineplot(data=phase_all, x=temperature, y=magnet, hue=lattice)
plt.savefig("magnetization.pdf")
