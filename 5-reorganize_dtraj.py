import numpy as np
import pickle
import h5py as h5
import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
os.chdir('/home/tnguyen/FAH_MSM_4/clustering')


def transf(n_cl, cl_num):  # Transform discrete dtrajs into continous dtrajs assignments
    dt = h5.File(f'2nd_clustering/cluster_{n_cl}/fix_distance_dtraj.h5', 'r')
    if n_cl == 1:
        with h5.File(f'2nd_clustering/cluster_{n_cl}/distance_dtraj_transformed.h5', 'a') as f:
            for i, vl in enumerate(dt):
                dtraj = (np.array(dt[vl]) + cl_num[n_cl - 1])
                f.create_dataset(vl, data=dtraj, dtype=int)
    elif n_cl == 2:
        with h5.File(f'2nd_clustering/cluster_{n_cl}/distance_dtraj_transformed.h5', 'a') as f:
            for i, vl in enumerate(dt):
                dtraj = (np.array(dt[vl]) +
                         cl_num[n_cl - 1] + cl_num[n_cl - 2])
                f.create_dataset(vl, data=dtraj, dtype=int)


def t_ind(data, cl_ind):  # Transform indexes of each data in cluster into traj
    traj = []
    if cl_ind in data[:, 0]:
        ls_ind = data[:, 0] == cl_ind
        pos = np.where(ls_ind)[0]
        for j, vl in enumerate(pos):
            traj.append(data[vl][1])
        traj = np.array(traj)
    return traj


def order(data):  # Get the new order of original data
    full_trajs = []
    for j, traj in enumerate(data):
        order = []
        for i in range(len(traj)):
            order.append(np.where(traj == i)[0])
        order = np.concatenate(order)
        full_trajs.append(order)
    return np.array(full_trajs)


dt_1 = h5.File('2nd_clustering/cluster_0/distance_dtraj.h5', 'r')
dt_2 = h5.File('2nd_clustering/cluster_1/distance_dtraj.h5', 'r')
dt_3 = h5.File('2nd_clustering/cluster_2/distance_dtraj.h5', 'r')

tica_cl2 = h5.File('../tica_data/cluster_traj/1st_cluster/cluster_1.h5', 'r')
tica_cl3 = h5.File('../tica_data/cluster_traj/1st_cluster/cluster_2.h5', 'r')

keys_2 = np.array(tica_cl2)
keys_3 = np.array(tica_cl3)

with h5.File(f'2nd_clustering/cluster_1/fix_distance_dtraj.h5', 'w') as f:
    for i, vl in enumerate(dt_2):
        data = np.array(dt_2[vl])
        f.create_dataset(keys_2[i], data=data, dtype=int)

with h5.File(f'2nd_clustering/cluster_2/fix_distance_dtraj.h5', 'w') as g:
    for j, n in enumerate(dt_3):
        data = np.array(dt_3[n])
        g.create_dataset(keys_3[j], data=data, dtype=int)

a = [2480, 10, 10]
for i in range(1, 3):
    transf(i, a)

ind = pickle.load(
    open('1st_clustering_3_mstates/distance_cluster_indexes.pickle', 'rb'))
dtraj = h5.File('1st_clustering_3_mstates/dtrajs.h5', 'r')

y = []
for i in range(len(dtraj)):
    y.append(np.concatenate([t_ind(ind[0], i)] +
             [t_ind(ind[1], i)]+[t_ind(ind[2], i)]))
y = np.array(y)

d = order(y)  # correct order

dtrj_1 = h5.File(f'2nd_clustering/cluster_0/distance_dtraj.h5', 'r')
dtrj_2 = h5.File(
    f'2nd_clustering/cluster_1/distance_dtraj_transformed.h5', 'r')
dtrj_3 = h5.File(
    f'2nd_clustering/cluster_2/distance_dtraj_transformed.h5', 'r')

with h5.File(f'2nd_clustering/origin_dtraj.h5', 'a') as f:
    for i, trj in enumerate(dtrj_1):
        if (trj in dtrj_2.keys()) == True and (trj in dtrj_3.keys()) == True:
            data = np.concatenate(
                [np.array(dtrj_1[trj])] + [np.array(dtrj_2[trj])] + [np.array(dtrj_3[trj])])
        elif (trj in dtrj_2.keys()) == True and (trj in dtrj_3.keys()) == False:
            data = np.concatenate(
                [np.array(dtrj_1[trj])] + [np.array(dtrj_2[trj])])
        elif (trj in dtrj_2.keys()) == False and (trj in dtrj_3.keys()) == True:
            data = np.concatenate(
                [np.array(dtrj_1[trj])] + [np.array(dtrj_3[trj])])
        elif (trj in dtrj_2.keys()) == False and (trj in dtrj_3.keys()) == False:
            data = np.concatenate([np.array(dtrj_1[trj])])
        f.create_dataset(trj, data=data, dtype=int)

full_dtraj = h5.File(f'2nd_clustering/origin_dtraj.h5', 'r')
with h5.File(f'2nd_clustering/final_dtraj.h5', 'a') as l:
    for m, n in enumerate(full_dtraj):
        data = np.array(full_dtraj[n])[d[m]]  # correct order
        l.create_dataset(n, data=data, dtype=int)

print('Correction indexes of traj in dtraj file and reorganize indice of each data point matching to original F@H trajectories are all done')
