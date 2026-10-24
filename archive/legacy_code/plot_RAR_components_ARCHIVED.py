import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm
import os
self_path = os.path.dirname(__file__)
from utils import load_sparc_table, vel2acc

testing = False

table = load_sparc_table()
galaxy_count = 3 if testing else len(table["Galaxy"])
columns = [ "Rad", "Vobs", "errV", "Vgas",
            "Vdisk", "Vbul", "SBdisk", "SBbul" ]

for i in tqdm(range(galaxy_count), desc="SPARC galaxies"):
    g = table["Galaxy"][i]

    data_path = "./data/Newtonian_models/"+g+"_rotmod.dat"
    data_full_path = os.path.join(self_path, data_path)
    rawdata = np.loadtxt(data_full_path)
    data = pd.DataFrame(rawdata, columns=columns)

    r = data["Rad"].to_numpy()
    Vobs = data["Vobs"].to_numpy()
    Vgas, Vdisk, Vbul = data["Vgas"].to_numpy(), data["Vdisk"].to_numpy(), data["Vbul"].to_numpy()
    bulged = np.any(Vbul != 0.0)

    g_obs = vel2acc(Vobs**2, r)
    g_gas, g_disk, g_bul = vel2acc(Vgas*np.abs(Vgas) * 1.0, r), vel2acc(Vdisk**2 * 0.5, r), 0.0
    if bulged: g_bul = vel2acc(Vbul**2 * 0.7, r)
    g_bar_SPARC = g_gas + g_disk + g_bul

    plt.figure()
    plt.plot(r, g_obs, marker='.', color='black', label='g_obs', alpha=0.6)
    plt.plot(r, g_bar_SPARC, marker='.', color='red', label='g_bar', alpha=0.6)
    plt.plot(r, g_gas, marker='.', color='green', label='g_gas', alpha=0.6)
    plt.plot(r, g_disk, marker='.', color='blue', label='g_disk', alpha=0.6)
    if bulged: plt.plot(r, g_bul, marker='.', color='orange', label='g_bul', alpha=0.6)
    plt.xlabel('Radius (kpc)')
    plt.ylabel('Acceleration (m/s^2)')
    plt.yscale('log')
    plt.title(g)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(self_path, f'./plots/acc_profiles/{g}.png'), dpi=150, bbox_inches='tight')
    plt.close()
