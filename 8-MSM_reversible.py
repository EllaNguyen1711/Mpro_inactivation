import pyemma.plots as mplt
import pyemma.msm as msm
import numpy as np
import pyemma

os.chdir('/home/tnguyen/FAH_MSM_4')

P2_nonrev = np.load('MSM/enspara/tprobs_40.npy')
M2_nonrev = msm.markov_model(P2_nonrev)
w_nonrev = M2_nonrev.stationary_distribution
# make reversible
C = np.dot(np.diag(w_nonrev), P2_nonrev)
Csym = C + C.T
P2 = Csym / np.sum(Csym,axis=1)[:,np.newaxis]
M2 = msm.markov_model(P2)
w = M2.stationary_distribution

if os.path.isdir('MSM/pyemma') == False:
	os.mkdir('MSM/pyemma')

pyemma.save('MSM/pyemma/msm_enspara_tprobs.pyemma')
