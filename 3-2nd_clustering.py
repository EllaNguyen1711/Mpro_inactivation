import numpy as np
import os

os.chdir('/home/tnguyen/FAH_MSM_4')


# Making directories and storing all the results of clustering 
if os.path.isdir('clustering/2nd_clustering') == False:
        os.mkdir('clustering/2nd_clustering')

for num in range (3):
	slurm_script = f'''
import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import os
import enspara
import h5py as h5
import pandas as pd
import numpy as np
import pyemma.plots as pyemma_plots
from pyemma.util.contexts import settings

os.chdir(\'/home/tnguyen/FAH_MSM_4\')
# Clustering 
# there are 3 clusters of data for the 2nd clustering [0, 1, 2]
cluster_numbers = [2480, 10, 10] 
stride = [5, 1, 1]
feature = coor.source(\'tica_data/cluster_traj/1st_cluster/cluster_{num}.h5\')

cluster = pyemma.coordinates.cluster_kmeans(feature, k=cluster_numbers[{num}], max_iter=50, stride = stride[{num}])

print (\'Clustering is done\')

if os.path.isdir(\'clustering/2nd_clustering/cluster_{num}\') == False:
        os.mkdir(\'clustering/2nd_clustering/cluster_{num}\')

path = \'/home/tnguyen/FAH_MSM_4/clustering/2nd_clustering/cluster_{num}/\' 
cluster.write_to_hdf5(path + \'distance_dtraj.h5\')
clustercenters = cluster.clustercenters
hf_clustercenters = h5.File(path + \'cluster_centers.h5\', \'w\')
hf_clustercenters.create_dataset(\'cluster_centers\', data=clustercenters)
hf_clustercenters.close()

cluster_indexes = cluster.index_clusters

with open(path + \'cluster_indexes.pickle\', \'wb\') as handle:
    pickle.dump(cluster_indexes, handle)

print (\'All clustering results are saved sucessfully\')'''
	print(slurm_script)
	F = open(f'clustering/2nd_clustering/cluster_{num}.py','w')
	F.write(slurm_script)
	F.close()
	slurm_script1 = f'''#!/bin/bash
#PBS -S /bin/bash
#PBS -o /home/tnguyen/FAH_MSM_4/clustering/2nd_clustering/cluster_{num}.log
#PBS -j oe
#PBS -l nodes=1:ppn=8,mem=30gb,walltime=300:00:00

module load miniconda/3
source activate bayes
date
free -h
python /home/tnguyen/FAH_MSM_4/clustering/2nd_clustering/cluster_{num}.py'''
	print(slurm_script1)
	F1 = open(f'clustering/2nd_clustering/cluster_{num}.job','w')
	F1.write(slurm_script1)
	F1.close()
	os.system(f'qsub clustering/2nd_clustering/cluster_{num}.job')

