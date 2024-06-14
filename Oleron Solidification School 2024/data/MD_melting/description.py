from dscribe.descriptors import MBTR, LMBTR
from ase.io import read
from glob import glob
import numpy as np
import subprocess

# Setup
mbtr = MBTR(
    species=["Cu"],
    geometry={"function": "distance"},
    grid={"min": 1.5, "max": 5, "n": 100, "sigma": 0.05},
    weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
    periodic=True,
    normalization="l2")

struct=read('Cu_10-10-10-melting.dump', format='lammps-dump-text',index=':')

for S in struct:
    S.symbols=['Cu' for i in range(1000)]

discripts=mbtr.create(struct, n_jobs=14)

print(len(discripts), len(discripts[0]))

np.save('mbtr_copper.npy', discripts)


#--------------------------------------------------------------------------

import pandas as pd
import numpy as np

def read_lammps_dump(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()
    data = []
    read_data = False
    for line in lines:
        if line.startswith("ITEM: ATOMS"):
            read_data = True
            continue
        if line.startswith("ITEM: TIMESTEP"):
            read_data = False
        if read_data:
            data.append(line.split())
    column_names = ["id", "type", "x", "y", "z", "pe"]
    df = pd.DataFrame(data, columns=column_names, dtype=float)
    pe_array = df["pe"].to_numpy()
    return pe_array

file_name = "dump.pe"
pe_array = read_lammps_dump(file_name)
a=pe_array.reshape((1001,1000))
np.save('all_atomic_energies.npy', a)
aux=[a[i][i] for i in range(1000)]
np.save('pe_Cuatoms_i_i.npy', aux)

#---------------------------------------------------------------------------

mbtr = LMBTR(
    species=["Cu"],
    geometry={"function": "distance"},
    grid={"min": 1.5, "max": 5, "n": 100, "sigma": 0.05},
    weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
    periodic=True,
    normalization="l2")

discripts=[]
for i in range(1000):
    print(float(i)/10.0)
    list_atom=mbtr.create(struct[i], n_jobs=14, centers=[i])
    discripts.append(list_atom[0][100:200])

print(len(discripts), len(discripts[0]))

np.save('lmbtr_Cuatoms_i_i.npy', discripts)

