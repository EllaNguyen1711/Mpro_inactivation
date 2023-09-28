import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os
import glob
import enspara
import h5py
import pandas as pd
import pyemma.plots as pyemma_plots
from pyemma.util.contexts import settings

os.chdir('/home/tnguyen/FAH_MSM_4')

pdb = '/home/tnguyen/fah/prot_masses.pdb'
projects_paths = '/home/tnguyen/fah_cut_data/'
projects_names = ['PROJ14234', 'PROJ14542', 'PROJ14543']

info = []
for i, name in enumerate(os.listdir(projects_paths+projects_names[0])):
    for j in range(5):
        traj = projects_paths + projects_names[0] + '/' + name + '/' + f'CLONE{j}.xtc'
        info.append(str(traj))

for r, cl1 in enumerate(glob.glob(projects_paths+projects_names[1]+'/*.xtc')):
    info.append(str(cl1))

for k, cl2 in enumerate(glob.glob(projects_paths+projects_names[2]+'/*.xtc')):
    info.append(str(cl2))

#print (info)

feat = coor.featurizer(pdb)
feat.add_backbone_torsions(cossin=True)
feat.add_sidechain_torsions(which=['chi1','chi2'], cossin=True)
stride = 1
data =coor.source(info, features=feat, stride=stride)

#data.write_to_hdf5('feature_data/chi1_2.h5')
print ('Featurization has been done')
print ('Tica calculation starts ...')
var_cutoff = 0.9 # adjust to find elbow of of cumulative kinetic variance
stride_tica = 1
tica = coor.tica(data = data, lag= 10, kinetic_map=False, commute_map=True, var_cutoff = 0.9, stride =stride_tica)

#tica.var_cutoff = var_cutoff

print('Number of dimensions saved is: ', tica.dimension())

with open(f'tica_data/chi1_2_cumvar_stride_{stride_tica}.npy','wb') as handle:
    np.save(handle, tica.cumvar)
print ('Tica calculation is done, data is being saved')
tica.write_to_hdf5(f'tica_data/chi1_2_stride_{stride_tica}.h5')
eig = tica.eigenvectors
np.save(f'tica_data/eigenvectors_stride_{stride_tica}.npy', eig)
