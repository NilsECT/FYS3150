import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
plt.rcParams.update({"font.size": 16})
sns.set_theme()

data = pd.read_fwf("Problem_10.txt")

plt.errorbar(np.log10(data["n"]), data["dt_mean"], yerr=data["dt_stddev"], fmt="o", label="General algorithm")
plt.errorbar(np.log10(data["n"]), data["dt_mean_spec"], yerr=data["dt_stddev_spec"], fmt="o", label="Special algorithm")
plt.legend()
plt.xlabel("log_10(n steps)")
plt.ylabel("Time (s)")
plt.show()