import numpy as np
import pickle
import os
import h5py as h5
from multiprocessing import Pool
import itertools
import gc


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


def traj_rw(info):
    import pickle
    import h5py as h5
    import os
    import numpy as np
    cluster = pickle.load(open(info[1], 'rb'))[info[3]]
    data = h5.File(info[2], 'r')[f'{info[0]:04d}']
    traj = []
    if info[0] in cluster[:, 0]:
        ls_ind = cluster[:, 0] == info[0]
        pos = np.where(ls_ind)[0]
        for i, vl in enumerate(pos):
            fr = data[cluster[vl][1]]
            traj.append(fr)
    return traj
path = '/home/tnguyen/FAH_MSM_4/'

if os.path.isdir(path + 'tica_data/cluster_traj') == False:
	os.mkdir(path + 'tica_data/cluster_traj')
if os.path.isdir(path + 'tica_data/cluster_traj/1st_cluster') == False:
    os.mkdir(path + 'tica_data/cluster_traj/1st_cluster')

cluster_numbers = 3
# path of cluster file
n_cluster = path + f'clustering/1st_clustering_{cluster_numbers}_mstates/distance_cluster_indexes.pickle'
n_data = path + 'tica_data/chi1_2_stride_1.h5'  # path of tica data file
n = 0  # number of cluster choosen for rewritting

if os.path.isdir(path + f'tica_data/cluster_traj/1st_cluster/cluster_{n}') == False:
    os.mkdir(path + f'tica_data/cluster_traj/1st_cluster/cluster_{n}')
os.chdir(path + f'tica_data/cluster_traj/1st_cluster/cluster_{n}')
dtraj = h5.File(n_data, 'r')
list_ind = range(len(dtraj))  # len of the trajectories of data
x = []
for ind in list_ind:
    # combination of iterable information for each trajectory needed to be written in 'n' cluster
    x += [[ind, n_cluster, n_data, n]]

# chuncked the interable information into 5 for each chunck
size = int(10)
for i, c in enumerate(chunked_iterable(x, size=size)):
    with Pool() as p:
        # expect to have a map of 10 each time for calculation
        traj_list = p.map(traj_rw, c)
        with h5.File(f'cl_{n}_trj_{i*size}_to_{i*size+len(c)-1}.h5', 'w') as f:
            for j, traj_dt in enumerate(traj_list):
                f.create_dataset(f'{j+i*size:04d}',  data=traj_dt, dtype='float')
        del traj_list
        gc.collect()
print(f'Re-writting tica data for cluster {n} is done')
