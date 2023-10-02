import sys
import numpy as np
import h5py as h5
import pickle
from multiprocessing import Pool
import matplotlib.pyplot as plt
import glob
import os
sys.path.insert(0, '/home/tnguyen/Tools')
from _Distance import *
# Distance calculation
os.chdir('/home/tnguyen/FAH_MSM_4')

# Load traj, top and population
pdb = '/home/tnguyen/fah/prot_masses.pdb'
projects_paths = '/home/tnguyen/fah_cut_data/'
projects_names = ['PROJ14234', 'PROJ14542', 'PROJ14543']
# Add residue index
res1 = 143 #Replace selected residue number 
res2 = 145 #Replace selected redigue number
info = []
for i, name in enumerate(os.listdir(projects_paths+projects_names[0])):
    for j in range(5):
        traj = projects_paths + \
            projects_names[0] + '/' + name + '/' + f'CLONE{j}.xtc'
        info.append(np.array([pdb, str(traj), res1, res2]))

for r, cl1 in enumerate(glob.glob(projects_paths+projects_names[1]+'/*.xtc')):
    info.append(
        np.array([pdb, str(cl1), res1, res2]))

for k, cl2 in enumerate(glob.glob(projects_paths+projects_names[2]+'/*.xtc')):
    info.append(
        np.array([pdb, str(cl2), res1, res2]))

print (info)

with h5.File(f'data/test_dis_oxyanion_NH{res1}_NH{res2}_full.h5', 'a') as f:
    for l, c in enumerate(chunked_iterable(info, size=10)):
        with Pool() as p:
            # expect to have a map of 10 each time for calculation
            dis_run = p.starmap(oxyanion_dis, c)
        for k, t in enumerate(dis_run):
            f.create_dataset(f'{l*10+k:04d}', data=t, dtype=float)
print('Distance calculation is done!')
