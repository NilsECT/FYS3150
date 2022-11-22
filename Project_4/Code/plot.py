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
energy = "Energy [J]"
sample = "Sample"
magnet = "Magnetisation"
burn = "Burn in [cycles]"
cv = "Heat capacity"
chi = "Susceptibility"
average_energy = "Average energy [J]"


###########  LOOKING AT BURN IN ORDERED VS UNORDERED #######

var_N = get_data("varying_cycles_no_burn_lattices_true")

lattices = [40]# np.array([20, 40, 60, 80, 100])
# var_N = var_N.loc[var_N[cycles] > 3000]

for lat in lattices:

    # using fontsizes 16

    # plt.figure(figsize=(12, 10))
    # ordered = var_N.loc[var_N["Order"] == "Ordered"]
    # sns.lineplot(data=ordered.loc[ordered[lattice] == lat], x=cycles, y=average_energy, hue=temperature)
    # plt.xscale("log")
    # plt.savefig("ordered_var_MC_%d.pdf" %lat)
    # # plt.show()

    # plt.figure(figsize=(12, 10))
    # unordered = var_N.loc[var_N["Order"] == "Unordered"]
    # sns.lineplot(data=unordered.loc[unordered[lattice] == lat], x=cycles, y=average_energy, hue=temperature)
    # plt.xscale("log")
    # plt.savefig("unordered_var_MC_%d.pdf" %lat)
    # plt.show()

    plt.figure(figsize=(12, 10))
    sns.lineplot(data=var_N, x=cycles, y=average_energy, hue=temperature, style="Order")
    plt.xscale("log")
    plt.savefig("var_MC_e_%d_true.pdf" %lat)
    # plt.show()

###########################################################

################### LOOKING AT ENERGY DISTRINUTION #######

eps_dist = get_data("epsilon_distribution")

#using fontsizes 12

temp = [1., 1.5, 2., 2.4]

for i in temp:
    data = eps_dist.loc[eps_dist[temperature] == i]
    for lat in lattices:

        plt.figure(figsize=(15, 13))
        sns.displot(data=data.loc[data[lattice] == lat], x="Energy [J]", stat="probability", bins=75, kde=True)
        # title_temp_lattice_#MC_#samples
        plt.savefig("energy_dist_%.d_%d.pdf" %(i*10, lat))
        plt.show()

#######################################################

############# LOOKING AT PHASE TRANSITIONS ###############

phase = get_data("phase_transition_varL_100")
phase_2 = get_data("phase_transition_varL_alt_frem_til_100_feil_størrelseorden")
phase_2[cv] = phase_2[cv]*phase_2[lattice]**4
phase_2[chi] = phase_2[chi]*phase_2[lattice]**4

# Specific heat capacity
plt.figure()
sns.lineplot(data=phase, x=temperature, y=cv)
sns.lineplot(data=phase_2, x=temperature, y=cv, hue=lattice)
plt.savefig("cv.pdf")

# Susceptibility 
plt.figure()
sns.lineplot(data=phase, x=temperature, y=chi)
sns.lineplot(data=phase_2, x=temperature, y=chi, hue=lattice)
plt.savefig("chi.pdf")



# Problem 9:
data_9 = get_data("21_26_1sample_phase_transition_varL")

T_c1 = []
for lat in lattices:
    i = data_9.loc[data_9[lattice] == lat][chi].idxmax()
    T_c1.append(data_9[temperature].to_numpy()[i])


# print(T_c1)

data_10 = get_data("phase_transition_varL_FRIDAS_2_kveld (1)")
T_c2 = []
for lat in lattices:
    i = data_10.loc[data_10[lattice] == lat][chi].idxmax()
    T_c2.append(data_9[temperature].to_numpy()[i])

# print(T_c2)

lattices_2468 = np.array([20, 40, 60, 80])
T_c3 = []
for lat in lattices_2468:
    i = phase_2.loc[phase_2[lattice] == lat][chi].idxmax()
    T_c3.append(phase_2[temperature].to_numpy()[i])
i = phase[chi].idxmax()
T_c3.append(phase[temperature].to_numpy()[i])
# print(T_c3)

data_nils = get_data("phase_transition_varL_nils")
T_c4 = []
for lat in lattices:
    i = data_nils.loc[data_nils[lattice] == lat][chi].idxmax()
    T_c4.append(data_nils[temperature].to_numpy()[i])

T_c_mean = np.zeros(len(T_c1))
for i in range(len(T_c1)):
    T_ci = np.array([T_c1[i], T_c2[i], T_c3[i], T_c4[i]])
    T_c_mean[i] = np.mean(T_ci)

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

phase_all_list = [phase, phase_2, data_9, data_10, data_nils]
phase_all = pd.concat(phase_all_list, ignore_index=True)

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

# sns.lineplot(data=var_N, x=cycles, y=energy, hue=temperature, style="Order")