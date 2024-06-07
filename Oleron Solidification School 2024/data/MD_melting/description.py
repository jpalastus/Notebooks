from dscribe.descriptors import MBTR
from ase.io import read
from glob import glob
import numpy as np
import subprocess

# Setup descriptor
mbtr = MBTR(
    species=["Cu"],
    geometry={"function": "distance"},
    grid={"min": 1.5, "max": 5, "n": 100, "sigma": 0.05},
    weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
    periodic=True,
    normalization="l2")
discripts=[]

# Readinf trajectory with ASE
struct=read('Cu_10-10-10-melting.dump', format='lammps-dump-text',index=':')

# Changing LAMMPS labels '1' to 'Cu' 
for S in struct:
    S.symbols=['Cu' for i in range(1000)]

# Computing descriptors
discripts=mbtr.create(struct, n_jobs=14)

# Saving results
np.save('mbtr_copper.npy', discripts)
