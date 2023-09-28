import numpy as np
import pickle
import os
import h5py as h5
from multiprocessing import Pool
import itertools
import gc


os.chdir('/home/tnguyen/FAH_MSM_4')
num_fr = 50
if os.path.isdir(f'Dihedral_{num_fr}') == False:
    os.mkdir(f'Dihedral_{num_fr}')
if os.path.isdir(f'Dihedral_{num_fr}/jobs') == False:
    os.mkdir(f'Dihedral_{num_fr}/jobs')
group_clusters = np.load(
    '/home/tnguyen/FAH_MSM_4/MSM/pyemma/pathways_clusters_coarse_grained.npy', allow_pickle=True)
states = ['A-I', 'Inter1', 'Inter2', 'Inter3', 'I-A']
for l, gr in enumerate(group_clusters):
    slurm_script = f'''import sys
import numpy as np
import os 
import MDAnalysis as mda
import h5py as h5
from MDAnalysis.analysis.dihedrals import Ramachandran, Janin
sys.path.insert(0, \'/home/tnguyen/Tools\')
from _Dihedral import *


os.chdir(\'/home/tnguyen/FAH_MSM_4\')

states = [\'A-I\', \'Inter1\', \'Inter2\', \'Inter3\', \'I-A\']
if os.path.isdir(\'Dihedral_{num_fr}/%s\'%states[{l}]) == False:
	os.mkdir(\'Dihedral_{num_fr}/%s\'%states[{l}])

pdb = \'/home/tnguyen/fah/prot_masses.pdb\'
path_traj = \'/home/tnguyen/FAH_MSM_4\'
group_clusters = np.load(\'/home/tnguyen/FAH_MSM_4/MSM/pyemma/pathways_clusters_coarse_grained.npy\', allow_pickle = True)
order_clusters = group_clusters[{l}]
n_state = \'%s\'%states[{l}]
trajs = [path_traj +\'/Samples_{num_fr}/cluster_%s.xtc\'%i for i in order_clusters]
cal = Dihedral(pdb, trajs)
cal.write_resdues_h5(n_residues=\'all\', path = \'Dihedral_{num_fr}/%s\'%states[{l}])
print(\'Done!\') 
'''
    print(slurm_script)
    F = open(f'Dihedral_{num_fr}/jobs/{states[l]}.py', 'w')
    F.write(slurm_script)
    F.close()
    slurm_script1 = f'''#!/bin/bash
#PBS -S /bin/bash
#PBS -o /home/tnguyen/FAH_MSM_4/Dihedral_{num_fr}/jobs/{states[l]}.log
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=5gb,walltime=200:00:00

module load miniconda/3
source activate bayes
date
free -h
python /home/tnguyen/FAH_MSM_4/Dihedral_{num_fr}/jobs/{states[l]}.py'''
    print(slurm_script1)
    F1 = open(f'Dihedral_{num_fr}/jobs/{states[l]}.job', 'w')
    F1.write(slurm_script1)
    F1.close()
    os.system(f'qsub Dihedral_{num_fr}/jobs/{states[l]}.job')
