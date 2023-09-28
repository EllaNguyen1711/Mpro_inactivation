import mdtraj as md
import numpy as np
from multiprocessing import Pool
from MDAnalysis.analysis import align
import glob
import os
import sys
sys.path.insert(0, '/home/tnguyen/Tools')
from _Distance import *
os.chdir('/home/tnguyen')

def align_trajs(pdb, traj, name):
	ref = mda.Universe(pdb)
	mobile = mda.Universe(pdb, traj)
	align.AlignTraj(mobile, ref, select='all', filename = name).run() 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--reference', type = str, default = None, help='Location of reference pdb file')
parser.add_argument('--trajectory', type = str, default= None, help='Location of trajectory .xtc, .dcd, .pdb files')
parser.add_argument('--prefix', type = str, default = None, help='Path to Directory storing output files')
parser.add_argument('--project', type = str, default = None, help='Path to project of trajectories')
args = parser.parse_args()

if not os.path.isfile(args.reference):
  raise Exception('Reference missing or is not a file!')

if args.trajectory == None:
	path_to_projects = '/home/tnguyen/fah_cut_data'
	projects_names = ['PROJ14234', 'PROJ14542', 'PROJ14543']
	if args.project == projects_names[0]:
		trajs = glob.glob(path_to_projects+'/'+args.project+'/*/*.xtc')
	else:
		trajs = glob.glob(path_to_projects+'/'+args.project+'/*.xtc')
	ref = md.load(args.reference)
	if os.path.isdir(args.prefix + '/' + args.project) == False:
		os.mkdir(args.prefix + '/' + args.project)
	info = []
	for i, path in enumerate(trajs):
		info.append([args.reference, path, args.prefix+ '/' + args.project+ f'/CLONE{i}.xtc'])
	for l, c in enumerate(chunked_iterable(info, size=10)):
		with Pool() as p:
			dis_run = p.starmap(align_trajs, c)
			
else: 
	trajs = args.trajectory
	info = []
	for i, path in enumerate(trajs):
		info.append([args.reference, path, args.prefix+ '/' + f'traj{i}.xtc'])
	for l, c in enumerate(chunked_iterable(info, size=10)):
		with Pool() as p:
			dis_run = p.starmap(align_trajs, c)
print ('Alignment of trajectories is done!')
