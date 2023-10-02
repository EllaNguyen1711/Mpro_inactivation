import mdtraj as md
import pyemma.coordinates as coor
import numpy as np
import pickle
import pyemma
import numbers
import os
import enspara
import h5py
import pandas as pd
import numpy as np
import logging
import csv
import scipy
import scipy.sparse
import scipy.sparse.linalg
from scipy.sparse.csgraph import connected_components
from enspara import exception
from enspara.msm import MSM, builders, transition_matrices
from enspara.msm.transition_matrices import _transitions_helper, assigns_to_counts
import pyemma.plots as pyemma_plots
from pyemma.util.contexts import settings


os.chdir('/home/tnguyen/FAH_MSM_4')

#Building MSM
msm_lag = 40
cluster_numbers = 2500

cluster_assig = coor.load(f'clustering/2nd_clustering/final_dtraj.h5')
cluster_assig = np.array(cluster_assig)
assig = []
for frame in cluster_assig:
	assig.append(frame.astype(int))
assig = np.array(assig)
assig = assig.reshape(1, -1)
print ('MSM is being built ...')
tcounts = assigns_to_counts(assig, lag_time = msm_lag)
prior_counts = 1/tcounts.shape[0]
tcounts = builders._apply_prior_counts(tcounts, prior_counts)
probs = builders._row_normalize(tcounts)
eq_probs_ = builders.eq_probs(probs)

print ('transition maxtrix: ', tcounts)
print ('transition probabilities: ', probs)
print ('equilibrium probabilities: ', eq_probs_)
if os.path.isdir('MSM') == False:
	os.mkdir('MSM')
if os.path.isdir('MSM/enspara') == False:
	os.mkdir('MSM/enspara')
np.save (f'MSM/enspara/tcounts_{msm_lag}.npy', tcounts)
np.save (f'MSM/enspara/tprobs_{msm_lag}.npy', probs)
np.save (f'MSM/enspara/populations_{msm_lag}.npy', eq_probs_)

print('MSM was built successfully')
