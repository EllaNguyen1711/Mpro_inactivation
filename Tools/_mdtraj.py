import sys
import numpy as np
from multiprocessing import Pool
import MDAnalysis as mda
import itertools
import os


def chunked_iterable(iterable, size):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


def writeInput_Figngerprint_Cal(name):
    # name: name of the file PDB or path to that PDB file
    with open(name[:-4]+'.txtPOVME.in', 'w') as fn:

        info = ['GridSpacing 0.5',  # GridSpacing
                'PointsInclusionSphere 58.426  69.482  28.203 6',  # PIS1 coords of S in A-Cys145
                'PointsInclusionSphere 92.127  87.957  47.868 6',  # PIS2 coords of S in B-Cys145
                'DistanceCutoff 1.09',
                'OutputFilenamePrefix '+'./data/'+name[:-4]+'_',
                'CompressOutput false',
                'PDBFileName '+name,
                'ConvexHullExclusion true',
                'NumProcessors 1', 
		'UseDiskNotMemory false',
		'SaveIndividualPocketVolumes false',
		'SavePocketVolumesTrajectory false',
		'OutputEqualNumPointsPerFrame false',
		'SaveTabbedVolumeFile false',
		'SaveVolumetricDensityMap false',
		'SaveFingerprint true'
                ]

        for key in info:
            fn.write(key+'\n')


def Fingerprint_Cal(name):
    writeInput_Figngerprint_Cal(name)
    input_file = name[: -4]+'.txtPOVME.in'
    command= 'python /home/tnguyen/FAH_MSM_5/POVME_2_0_1/POVME2.py '+input_file
    os.system(command)
    output_file = name[: -4]+'_output.txt'
    os.remove(input_file)
    # os.remove(name)


def traj_to_pdbs_1(top, traj):
    u= mda.Universe(top, traj)
    for i in range(len(u.trajectory)):
        u.select_atoms("protein").write(
            'frame_%s.pdb' % i, frames=u.trajectory[[i]])

def write_selpdb(info): #info = [top, traj, ind]
    u= mda.Universe(info[0], info[1])
    u.select_atoms("protein").write(
            'frame_%s_%s.pdb' % (info[2][0], info[2][1]), frames=u.trajectory[info[2][0]:info[2][1]])

def multi_processing_traj_to_pdbs(top, traj, stride):
	u= mda.Universe(top, traj)
	info = []
	range_x = range(len(u.trajectory))[::20]
	for i, r in enumerate(range_x):
		if i < len(range_x)-1:
			info.append([top, traj, [range_x[i], range_x[i+1]]])
		else:
			info.append([top, traj, [range_x[i], len(u.trajectory)]])
	pool = Pool()
	pool.map(write_selpdb, info)
	pool.close()

class MD_load():
    '''This code will be used for working with a meta-data of MD simulations'''

    def __init__(self, top, trajs):
        self.top = top
        self.trajs = trajs

    def traj_to_pdbs(ind, prefix):  # Should align all trajs before conversion
        u = mda.Universe(self.top, self.trajs[ind])
        for i in range(len(u.trajectory)):
            u.select_atoms("protein").write(
                prefix+'frame_%s.pdb' % i, frames=u.trajectory[[i]])

    def traj_to_pdb(ind, prefix):  # Should align all trajs before conversion
        if ind == 'all':
            for n, traj in enumerate(self.trajs):
                u = mda.Universe(self.top, traj)
                u.select_atoms("protein").write(
                    prefix+'traj_%s.pdb' % n, frames='all')
        else:
            for n in ind:
                u = mda.Universe(self.top, self.trajs[n])
                u.select_atoms("protein").write(
                    prefix+'traj_%s.pdb' % n, frames='all')
