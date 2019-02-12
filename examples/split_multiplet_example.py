import sys
sys.path.append('../')
from colradpy_class import *
import numpy as np

he = colradpy('./mom97_ls#he1.dat',[0],np.array([20]),np.array([1.e13]),use_recombination=False,use_recombination_three_body = False,use_ionization=True)

he.solve_cr()
he.split_pec_multiplet()

wave_8_3 = np.array([468.5376849,468.5757974,468.5704380])
ind_8_3 = np.where( (he.data['processed']['pec_levels'][:,0] == 8) & (he.data['processed']['pec_levels'][:,1] == 3))[0]

wave_6_5 = np.array([468.5407225,468.5568006])
ind_6_5 = np.where( (he.data['processed']['pec_levels'][:,0] == 6) & (he.data['processed']['pec_levels'][:,1] == 5))[0]

wave_7_3 = np.array([468.5524404,468.5905553])
ind_7_3 = np.where( (he.data['processed']['pec_levels'][:,0] == 7) & (he.data['processed']['pec_levels'][:,1] == 3))[0]

wave_9_4 = np.array([468.5703849, 468.5830890, 468.5804092])
ind_9_4 = np.where( (he.data['processed']['pec_levels'][:,0] == 9) & (he.data['processed']['pec_levels'][:,1] == 4))[0]

wave_6_4 = np.array([ 468.5917884, 468.5757080, 468.5884123])
ind_6_4 = np.where( (he.data['processed']['pec_levels'][:,0] == 6) & (he.data['processed']['pec_levels'][:,1] == 4))[0]


wave_468 = np.hstack((wave_8_3,wave_6_5,wave_7_3,wave_9_4,wave_6_4))
pecs_468 = np.vstack((he.data['processed']['split']['pecs'][ind_8_3[0]],
                      he.data['processed']['split']['pecs'][ind_6_5[0]],
                      he.data['processed']['split']['pecs'][ind_7_3[0]],
                      he.data['processed']['split']['pecs'][ind_9_4[0]],
                      he.data['processed']['split']['pecs'][ind_6_4[0]]))[np.argsort(wave_468)]
wave_468 = wave_468[np.argsort(wave_468)]



plt.figure()
plt.vlines(wave_468,np.zeros_like(wave_468),pecs_468[:,0,0,0])
