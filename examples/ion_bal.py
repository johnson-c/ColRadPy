import sys
sys.path.append('../')
from colradpy_class import colradpy
import numpy as np
import matplotlib.pyplot as plt
from ionization_balance_class import ionization_balance


fils = np.array(['cpb03_ls#be0.dat','cpb03_ls#be1.dat','be2_adf04','be3_adf04'])
temp = np.linspace(3,100,20)
dens = np.array([1.e9,1.e10,1.e11,1.e12,1.e13,1.e14,1.e15,1.e16])
metas = [np.array([0,1]),np.array([0]),np.array([0,1]),np.array([0])]

ion = ionization_balance(fils, metas, temp, dens, keep_species_data = True)

ion.populate_ion_matrix()

ion.solve_no_source(np.array([1,0,0,0,0,0,0]),np.linspace(0,10,1e5))












#ion.solve_source(np.array([1,0,0,0,0,0,0]),np.array([1,0,0,0,0,0,0]),np.linspace(0,10,1e5))
