import MDAnalysis as mda
import mdtraj as md
from MDAnalysis.analysis import align
import os
import numpy as np

def fix_traj(ref, traj, indices):
    sim = mda.Universe(ref, traj)
    ref_top = mda.Universe(ref)
    ls = range(np.min(indices), np.max(indices)+1)
    for i, vl in enumerate(ls):
        ag = sim.select_atoms('protein')
        ag.write(f'fr_{vl}.xtc', frames = sim.trajectory[vl:vl+1])
    frames = [f'fr_{vl}.xtc' for vl in ls]
    for j in frames:
        fix_frame(ref, j)
    fix_traj_ls = [f'{frame[:int(len(frame)-4)]}.pdb' for frame in frames]
    fix_traj = md.load(fix_traj_ls)
    fix_traj.superpose(md.load(ref))
    fix_traj.save_xtc('fixed_traj.xtc')
    ags = sim.select_atoms('protein')
    ags.write('before.xtc', frames = sim.trajectory[:np.min(indices)])
    ags.write('after.xtc', frames = sim.trajectory[np.max(indices)+1:])
    full_traj = mda.Universe(ref, ['before.xtc', 'fixed_traj.xtc', 'after.xtc'])
    protein = full_traj.select_atoms('protein')
    protein.write(traj[:int(len(traj)-4)] + '_fixed.xtc', frames = full_traj.trajectory[:])
    os.remove('before.xtc')
    os.remove('fixed_traj.xtc')
    os.remove('after.xtc')
    for file in frames:
        os.remove(file)
        os.remove(f'{file[:int(len(file)-4)]}.pdb')

def fix_frame(ref, frame):
    ref_top = mda.Universe(ref)
    pdb = mda.Universe(ref, frame)
    mono_A = pdb.select_atoms('segid A')
    mono_B = pdb.select_atoms('segid B')
    align.alignto(mono_A, ref_top, select = 'segid A')
    mono_A.atoms.write('A.pdb')
    align.alignto(mono_B, ref_top, select = 'segid B')
    mono_B.atoms.write('B.pdb')
    n_A = mda.Universe('A.pdb')
    n_B = mda.Universe('B.pdb')
    fix_fr = mda.Merge(n_A.atoms, n_B.atoms)
    fix_fr.atoms.write(f'{frame[:int(len(frame)-4)]}.pdb')
    os.remove('A.pdb')
    os.remove('B.pdb')

