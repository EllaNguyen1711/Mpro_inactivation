import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
from pyemma.coordinates.clustering import AssignCenters
from pyemma.coordinates import pipeline
import os
import enspara
import h5py as h5
import shutil
import pandas as pd
import numpy as np
from enspara.msm import MSM, builders, transition_matrices
import pyemma.plots as pyemma_plots
from pyemma.util.contexts import settings

os.chdir('/home/tnguyen/FAH_MSM_4')
#Clustering 
centers = [[4.2, 4.2], [4.2, 10.0], [10.0, 4.2]]
cluster_numbers = len(centers)
stride = 1
data = h5.File('data/dis_Phe_His_full.h5', 'r')
feature = []
for i, traj in enumerate(data):
        feature.append(np.array(data[traj]))

#cluster = pyemma.coordinates.cluster_kmeans(feature, k=cluster_numbers, max_iter=50, stride = stride)

if os.path.isdir('clustering') == False:
	os.mkdir('clustering')
if os.path.isdir(f'clustering/1st_clustering_{cluster_numbers}_mstates') == False:
	os.mkdir(f'clustering/1st_clustering_{cluster_numbers}_mstates')
else: 
	shutil.rmtree(f'clustering/1st_clustering_{cluster_numbers}_mstates')
	os.mkdir(f'clustering/1st_clustering_{cluster_numbers}_mstates')


dtrajs = pyemma.coordinates.assign_to_centers(feature, centers, stride=1, return_dtrajs=True, metric='euclidean', n_jobs=None, skip=0)
with h5.File(f'clustering/1st_clustering_{cluster_numbers}_mstates/dtrajs.h5', 'w') as f:
	for i, vl in enumerate(dtrajs):
		f.create_dataset(f'{i:04d}', data = vl, dtype = int)

print ('Saving dtraj is done')
cluster_indexes = []
for i in range (len(centers)):
	cl = []
	for j, trj in enumerate(dtrajs):
		dtrj = []
		for n, vl in enumerate(trj):
			if i == vl:
				dtrj.append([j, n])
		cl.append(np.array(dtrj))
	cluster_indexes.append(np.array(cl))

with open(f'clustering/1st_clustering_{cluster_numbers}_mstates/distance_cluster_indexes.pickle', 'wb') as handle:
    pickle.dump(cluster_indexes, handle)

print ('All clustering results are saved sucessfully')
