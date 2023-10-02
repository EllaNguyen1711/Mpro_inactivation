import pyemma.coordinates as coor
import numpy as np
import enspara.msm as msm
import pyemma.plots as pyemma_plots
import matplotlib.pyplot as plt
import h5py 
import os 
def assigns_to_counts1(
        assigns, lag_time, max_n_states=None, sliding_window=True):
    """Count transitions between states in a single trajectory.
    Parameters
    ----------
    assigns : array, shape=(traj_len, )
        A 2-D array where each row is a trajectory consisting of a
        sequence of state indices.
    lag_time : int
        The lag time (i.e. observation interval) for counting
        transitions.
    max_n_states : int, default=None
        The number of states. This is useful for controlling the
        dimensions of the transition count matrix in cases where the
        input trajectory does not necessarily visit every state.
    sliding_window : bool, default=True
        Whether to use a sliding window for counting transitions or to
        take every lag_time'th state.
    Returns
    -------
    C :  array, shape=(n_states, n_states)
        A transition count matrix.
    """

    if not isinstance(lag_time, numbers.Integral):
        raise exception.DataInvalid(
            "The lag time must be an integer. Got %s type %s." %
            lag_time, type(lag_time))
    if lag_time < 1:
        raise exception.DataInvalid(
            "Lag times must be be strictly greater than 0. Got '%s'." %
            lag_time)
    '''
    # if it's 1d, later stuff will fail
    if len(assigns.shape) == 1:
        raise exception.DataInvalid(
            'The given assignments array has 1-dimensional shape %s. '
            'Two dimensional shapes = (n_trj, n_frames) are expected. '
            'If this is really what you want, try using '
            'assignments.reshape(1, -1) to create a single-row 2d array.')
    '''
    assigns = np.array([a[np.where(a != -1)] for a in assigns])

    if max_n_states is None:
        max_n_states = np.concatenate(assigns).max() + 1

    transitions = [
        _transitions_helper(
            assign, lag_time=lag_time, sliding_window=sliding_window)
        for assign in assigns]
    # generate sparse matrix
    mat_coords = np.hstack(transitions)
    mat_data = np.ones(mat_coords.shape[1], dtype=int)
    C = scipy.sparse.coo_matrix(
        (mat_data, mat_coords), shape=(max_n_states, max_n_states))
    return C

if os.path.isdir('/home/tnguyen/FAH_MSM_4/implied_timescales') == False:
	os.mkdir('/home/tnguyen/FAH_MSM_4/implied_timescales')
os.chdir('/home/tnguyen/FAH_MSM_4/implied_timescales')

dtrajs = coor.load(f'../clustering/2nd_clustering/final_dtraj.h5')

dt2 = [i.astype(np.int_) for i in dtrajs]
dt3 = [i.reshape((i.shape[0])) for i in dt2]
lags = [1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70]
print('Calculation is starting ...')
n_clusters = 2500
def norm_pseudo(C, prior_counts=1/ n_clusters, calculate_eq_probs=True):
    return msm.builders.normalize(C, prior_counts=prior_counts, calculate_eq_probs=calculate_eq_probs)

if os.path.isfile('implied_timescales.npy') == False:
	its = msm.timescales.implied_timescales(np.array(dt3), lags, norm_pseudo, n_times=8)
	np.save('implied_timescales.npy', its)
else:
	its = np.load('implied_timescales.npy')
print ('Caculation is done')

cm = 1/2.54
ft = 50
fig, ax = plt.subplots(figsize=(8.5/cm, 6.5/cm))

for i in range(8):
    ax.plot(np.array(lags)*0.1, np.absolute(its[:,i]), linewidth=10, marker='o', markersize = 15)

ax.plot(np.array(lags)*0.1, np.array(lags)*0.1, linestyle='dashed', linewidth=10, color='k')
ax.set_yscale('log')
ax.set_xlabel('Lag time (ns)', fontsize = ft)
ax.set_ylabel('Implied Timescales', fontsize = ft)
ax.tick_params(axis='both', which='major', labelsize=ft)
fig.savefig(f'implied_timescale_enspara_norm_psuedo.png', bbox_inches ="tight", pad_inches = 0.1, transparent = True, facecolor ="w", edgecolor ='w', orientation ='landscape')
