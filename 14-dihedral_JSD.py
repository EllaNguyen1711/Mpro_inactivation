import os
import numpy as np

os.chdir('/home/tnguyen/FAH_MSM_4')
chain = ['A', 'B']
num_frs = 2000
if os.path.isdir(f'Jensen_Shannon_{num_frs}') == False:
    os.mkdir(f'Jensen_Shannon_{num_frs}')

if os.path.isdir(f'Jensen_Shannon_{num_frs}/data') == False:
    os.mkdir(f'Jensen_Shannon_{num_frs}/data')

if os.path.isdir(f'Jensen_Shannon_{num_frs}/jobs') == False:
    os.mkdir(f'Jensen_Shannon_{num_frs}/jobs')

name_resids = ['resid%d' % u for u in range(1, 307)]

for k in chain:
    if k == 'A':
        for i, vl in enumerate(name_resids):
            slurm_script = f'''import numpy as np
import os
from scipy.stats import *
import h5py as h5
from scipy.stats import entropy
from numpy.linalg import norm
from multiprocessing import Pool
import MDAnalysis as mda
import itertools
import gc
import glob

os.chdir(\'/home/tnguyen/FAH_MSM_4\')

def norm_const(f):
    N, error = integrate.quad(lambda x: f(x), -np.pi, np.pi)
    return N

def js(P, Q):
    _P = P
    _Q = Q
    _M = 0.5 * (_P + _Q)
    return entropy(_M) - 0.5*(entropy(_P)+entropy(_Q))

def JSD(x1, x2, weight1, weight2):
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    kde_x1 = gaussian_kde(X1, weights=np.concatenate([weight1]*3), bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights=np.concatenate([weight2]*3), bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 10000)
    N1 = norm_const(kde_x1)
    N2 = norm_const(kde_x2)
    print (N1, N2)
    z1 = kde_x1(points)/N1
    z2 = kde_x2(points)/N2
    jsd = js(z1, z2)
    return jsd, N1, N2

def weight_samples(ls_clusters, f_probs, cluster_indexes, path):
	cl_probs = f_probs[ls_clusters]
	x = []
	for j, vl in enumerate(ls_clusters):
		cl_ind = cluster_indexes[vl]
		traj = path + \'/cluster_%s.xtc\'%vl
		u = mda.Universe(traj)
		x.append([cl_probs[j]/len(cl_ind)]*len(u.trajectory))
	x = np.concatenate(x)
	return x		
			
states = [\'A-I\', \'Inter1\', \'Inter2\', \'Inter3\', \'I-A\']
group_clusters = np.load(\'/home/tnguyen/FAH_MSM_4/MSM/pyemma/pathways_clusters_coarse_grained.npy\', allow_pickle = True)
population = np.load(\'MSM/pyemma/populations_40.npy\')
cluster_indexes = np.load(\'/home/tnguyen/FAH_MSM_4/clustering/2nd_clustering/final_cluster_indexes.pickle\', allow_pickle = True)
data = []
for st in states:
	dt = h5.File(\'Dihedral_{num_frs}/%s/protomer_{k}.h5\'%st, 'r')
	data.append(np.array(dt[\'{vl}\']))
list_torsion = [\'phi\', \'psi\', \'chi1\', \'chi2\']
path = \'/home/tnguyen/FAH_MSM_4/Samples_{num_frs}\'
with h5.File(\'Jensen_Shannon_{num_frs}/data/{vl}.h5\', \'a\') as handle:
	for t in range(len(data[0])):
		x = []
		for r in range(len(data)-1):
			weight1 = weight_samples(group_clusters[r], population, cluster_indexes, path)
			weight2 = weight_samples(group_clusters[r+1], population, cluster_indexes, path)
			x.append(JSD(data[r][t], data[r+1][t], weight1, weight2))
		handle.create_dataset(\'%s\'%list_torsion[t], data = x, dtype = float)
print(\'Done\')'''
            print(slurm_script)
            F = open(f'Jensen_Shannon_{num_frs}/jobs/{vl}.py', 'w')
            F.write(slurm_script)
            F.close()
            slurm_script1 = f'''#!/bin/bash
#PBS -S /bin/bash
#PBS -o /home/tnguyen/FAH_MSM_4/Jensen_Shannon_{num_frs}/jobs/{vl}.log
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=1gb,walltime=100:00:00

module load miniconda/3
source activate bayes
date
free -h
autopep8 -i /home/tnguyen/FAH_MSM_4/Jensen_Shannon_{num_frs}/jobs/{vl}.py
python /home/tnguyen/FAH_MSM_4/Jensen_Shannon_{num_frs}/jobs/{vl}.py'''
            print(slurm_script1)
            F1 = open(f'Jensen_Shannon_{num_frs}/jobs/{vl}.job', 'w')
            F1.write(slurm_script1)
            F1.close()
            os.system(f'qsub Jensen_Shannon_{num_frs}/jobs/{vl}.job')
    if k == 'B':
        for i1, vl1 in enumerate(range(307, 613)):
            slurm_script = f'''import numpy as np
import os
from scipy.stats import *
from scipy.stats import entropy
from numpy.linalg import norm
import h5py as h5
from multiprocessing import Pool
import itertools
import MDAnalysis as mda
import gc
import glob

os.chdir(\'/home/tnguyen/FAH_MSM_4\')

def norm_const(f):
    N, error = integrate.quad(lambda x: f(x), -np.pi, np.pi)
    return N

def js(P, Q):
    _P = P
    _Q = Q
    _M = 0.5 * (_P + _Q)
    return entropy(_M) - 0.5*(entropy(_P)+entropy(_Q))

def JSD(x1, x2, weight1, weight2):
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    kde_x1 = gaussian_kde(X1, weights=np.concatenate([weight1]*3), bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights=np.concatenate([weight2]*3), bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 10000)
    N1 = norm_const(kde_x1)
    N2 = norm_const(kde_x2)
    print (N1, N2)
    z1 = kde_x1(points)/N1
    z2 = kde_x2(points)/N2
    jsd = js(z1, z2)
    return jsd, N1, N2

def weight_samples(ls_clusters, f_probs, cluster_indexes, path):
        cl_probs = f_probs[ls_clusters]
        x = []
        for j, vl in enumerate(ls_clusters):
                cl_ind = cluster_indexes[vl]
                traj = path + \'/cluster_%s.xtc\'%vl
                u = mda.Universe(traj)
                x.append([cl_probs[j]/len(cl_ind)]*len(u.trajectory))
        x = np.concatenate(x)
        return x 
                        
states = [\'A-I\', \'Inter1\', \'Inter2\', \'Inter3\', \'I-A\']
group_clusters = np.load(\'/home/tnguyen/FAH_MSM_4/MSM/pyemma/pathways_clusters_coarse_grained.npy\', allow_pickle = True)
population = np.load(\'MSM/pyemma/populations_40.npy\')
cluster_indexes = np.load(\'/home/tnguyen/FAH_MSM_4/clustering/2nd_clustering/final_cluster_indexes.pickle\', allow_pickle = True)
data = []
for st in states:
        dt = h5.File(\'Dihedral_{num_frs}/%s/protomer_{k}.h5\'%st, 'r')
	name = \'resid{int(i1 + 1)}\'
        data.append(np.array(dt[\'%s\'%name]))
list_torsion = [\'phi\', \'psi\', \'chi1\', \'chi2\']
path = \'/home/tnguyen/FAH_MSM_4/Samples_{num_frs}\'
with h5.File(\'Jensen_Shannon_{num_frs}/data/resid{vl1}.h5\', \'a\') as handle:
        for t in range(len(data[0])):
                x = []
                for r in range(len(data)-1):
                        weight1 = weight_samples(group_clusters[r], population, cluster_indexes, path)
                        weight2 = weight_samples(group_clusters[r+1], population, cluster_indexes, path)
                        x.append(JSD(data[r][t], data[r+1][t], weight1, weight2))
               	handle.create_dataset(\'%s\'%list_torsion[t], data = x, dtype = float)
print(\'Done\')'''
            print(slurm_script)
            F = open(f'Jensen_Shannon_{num_frs}/jobs/resid{vl1}.py', 'w')
            F.write(slurm_script)
            F.close()
            slurm_script1 = f'''#!/bin/bash
#PBS -S /bin/bash
#PBS -o /home/tnguyen/FAH_MSM_4/Jensen_Shannon_{num_frs}/jobs/resid{vl1}.log
#PBS -j oe
#PBS -l nodes=1:ppn=1,mem=1gb,walltime=100:00:00

module load miniconda/3
source activate bayes
date
free -h
autopep8 -i /home/tnguyen/FAH_MSM_4/Jensen_Shannon_{num_frs}/jobs/resid{vl1}.py
python /home/tnguyen/FAH_MSM_4/Jensen_Shannon_{num_frs}/jobs/resid{vl1}.py'''
            print(slurm_script1)
            F1 = open(f'Jensen_Shannon_{num_frs}/jobs/resid{vl1}.job', 'w')
            F1.write(slurm_script1)
            F1.close()
            os.system(f'qsub Jensen_Shannon_{num_frs}/jobs/resid{vl1}.job')
