from dscribe.descriptors import MBTR, CoulombMatrix
from ase.io import read
from glob import glob
import numpy as np
import subprocess

# Setup the descriptors
cm = CoulombMatrix(n_atoms_max=147,permutation='eigenspectrum')
mbtr = MBTR(
    species=["Cu"],
    geometry={"function": "inverse_distance"},
    grid={"min": 0, "max": 1, "n": 100, "sigma": 0.1},
    weighting={"function": "exp", "scale": 0.5, "threshold": 1e-3},
    periodic=False,
    normalization="l2")

list_atom=[]
list_atom2=[]
energies=[]
for file in glob("Cu147-LM/*.xyz"):
  print(file)
  # Load geometries 
  struct=read(file, format='xyz')
  # Assures correct label for atoms
  struct.symbols=['Cu' for i in range(147)]
  # Compute descriptors
  list_atom2.append(mbtr.create(struct))
  list_atom.append(cm.create(struct))
  # Read energy from commentary line of xyz
  energies.append(float(subprocess.check_output("head -n 2 "+file+" | tail -n 1", shell=True)))

# Saving results
np.save('cm.npy', list_atom)
np.save('mbtr.npy', list_atom2)
np.save('energies.npy', energies)


