import pyemma.coordinates as coor
import numpy as np
import pyemma.msm as msm
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt
import h5py 
import os 
os.chdir('/home/tnguyen/FAH_MSM_4/implied_timescales')

dtrajs = coor.load(f'../clustering/2nd_clustering/final_dtraj.h5')

dt2 = [i.astype(np.int_) for i in dtrajs]
dt3 = [i.reshape((i.shape[0])) for i in dt2]
lags = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70]
print('Calculation is starting ...')
its = msm.its(dt3, lags=lags, nits=5)
print ('Caculation is done')
cm = 1/2.54
fig, ax = plt.subplots(figsize=(8.5/cm, 6.5/cm))
ax.set_xlabel('Lag time (ns)', fontsize = 45)
ax.set_ylabel('Implied Timescales (ns)', fontsize = 45)
#ax.set_title(f'Implied timescale max_lag using pyemma', fontsize = 45)
pyemma_plots.plot_implied_timescales(its, units='ns', dt = 0.1, ax=ax)

fig.savefig('implied_timescale_pyemma.pdf')
