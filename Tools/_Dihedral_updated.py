import pyemma
import numpy as np
import MDAnalysis as mda
import pyemma.coordinates as coor
import h5py as h5
import scipy.integrate as integrate
import os
from scipy.stats import *
from scipy.stats import entropy
import scipy.integrate as integrate

class Dihedral_dim:

    def __init__(self, pdb, trajs, n_residues = 'all'): #n_residues could be a list of residues ordered from 1
        self.pdb = pdb
        self.trajs = trajs
        self.n_residues = n_residues

    def feature(self):
        feat = coor.featurizer(self.pdb)
        feat.add_backbone_torsions(cossin=False)
        feat.add_sidechain_torsions(which=['chi1', 'chi2'], cossin=False)
        full = feat.active_features
        return full

    def get_mapping_dihedral_angles(self):
        full = self.feature()
        dt = []
        for d in full:
            dt.append(d.describe())
        dt = np.concatenate(dt)
        A = []  # mapping ind for dihedral angles protomer A
        B = []  # mapping ind for dihedral angles protomer B
        E = []  # list of sorted dihedral angles in protomer A
        H = []  # list of sorted dihedral angles in protomer B
        for n in range(1, 307):
            for i, vl in enumerate(dt):
                if int(vl[10:]) == int(n):
                    if int(vl[4:6]) == 0:
                        A.append(i)
                        E.append(vl)
                    elif int(vl[4:6]) == 1:
                        B.append(i)
                        H.append(vl)
        A = np.array(A)
        B = np.array(B)
        E = np.array(E)
        H = np.array(H)
        return A, E, B, H

    def get_selected_ordering_(self):
        lst_resid = self.n_residues
        mapping_ = self.get_mapping_dihedral_angles()
        lst = np.concatenate([mapping_[1], mapping_[3]])
        x = []
        for n, resid in enumerate(lst_resid):
            if resid < 307:
                for i, vl in enumerate(lst):
                    if int(vl[4:6]) == 0 and int(vl[10:]) == int(resid):
                        x.append(i)
            elif resid > 306:
                for i, vl in enumerate(lst):
                    if int(vl[4:6]) == 1 and int(vl[10:]) == int(resid-306):
                        x.append(i)
        x = np.array(x)
        return x

    def get_full_ordered_list_dihedral_angles(self):
        return np.concatenate([self.get_mapping_dihedral_angles()[1],
                               self.get_mapping_dihedral_angles()[3]])

    def get_dihedral_atoms_indexes(self):
        F = np.concatenate([self.get_mapping_dihedral_angles()[0],
                            self.get_mapping_dihedral_angles()[2]])
        full = self.feature()
        indexes = np.concatenate(
            [full[0].angle_indexes, full[1].angle_indexes])
        f_ind = indexes[F]
        return f_ind

    def get_angles(self, save_FN = None):
        if self.n_residues == 'all':
            f = coor.featurizer(self.pdb)
            f.add_dihedrals(self.get_dihedral_atoms_indexes())
            data = coor.source(self.trajs, features=f)
        elif self.n_residues != 'all':
            f = coor.featurizer(self.pdb)
            f.add_dihedrals(self.get_dihedral_atoms_indexes()
                            [self.get_selected_ordering_()])
            data = coor.source(self.trajs, features=f)
        if save_FN != None:
            data.write_to_hdf5(save_FN)
        return data

    def write_data(self, name):
        self.get_angles().write_to_hdf5(name)

    def get_data(self):
        return self.get_angles().get_output()

    def write_resdues_h5(self, path):
        top = mda.Universe(self.pdb)
        dt = np.concatenate(self.get_data())
        chain = ['A', 'B']
        if self.n_residues == 'all':
            ls_angles = self.get_full_ordered_list_dihedral_angles()
            for j, k in enumerate(chain):
                with h5.File(f'{path}/protomer_{k}.h5', 'a') as f:
                    for n, num in enumerate(range(1, int(len(top.residues)/2)+1)):
                        resid = []
                        for i, vl in enumerate(ls_angles):
                            if int(vl[4:6]) == j and int(vl[10:]) == int(n+1):
                                resid.append(dt[:, i])
                        resid = np.array(resid)
                        f.create_dataset(
                            f'resid{n+1}', data=resid, dtype=float)
        elif self.n_residues != 'all':
            ls_angles = self.get_full_ordered_list_dihedral_angles()[
                self.get_selected_ordering_()]
            for r, res_ind in enumerate(self.n_residues):
                with h5.File(f'{path}/selected_residues.h5', 'a') as f:
                    if res_ind > 306:
                        result = []
                        for i, vl in enumerate(ls_angles):
                            if int(vl[4:6]) == 1 and int(vl[10:]) == int(res_ind-306):
                                result.append(dt[:, i])
                        result = np.array(result)
                        f.create_dataset(
                            f'resid{res_ind}', data=result, dtype=float)
                    elif res_ind < 307:
                        result = []
                        for i, vl in enumerate(ls_angles):
                            if int(vl[4:6]) == 0 and int(vl[10:]) == int(res_ind):
                                result.append(dt[:, i])
                        result = np.array(result)
                        f.create_dataset(
                            f'resid{res_ind}', data=result, dtype=float)

class Dihedral_mono:

    def __init__(self, pdb, trajs, n_residues = 'all'):
        self.pdb = pdb
        self.trajs = trajs
        self.n_residues = n_residues
        ref = mda.Universe(pdb)
        self.residues = len(ref.residues)

    def feature(self):
        feat = coor.featurizer(self.pdb)
        feat.add_backbone_torsions(cossin=False)
        feat.add_sidechain_torsions(which=['chi1', 'chi2'], cossin=False)
        full = feat.active_features
        return full

    def get_mapping_dihedral_angles(self):
        full = self.feature()
        dt = []
        for d in full:
            dt.append(d.describe())
        dt = np.concatenate(dt)
        A = []  # mapping ind for dihedral angles protomer A
        B = []  # mapping ind for dihedral angles protomer B
        E = []  # list of sorted dihedral angles in protomer A
        H = []  # list of sorted dihedral angles in protomer B
        for n in range(1, self.residues+1):
            for i, vl in enumerate(dt):
                if int(vl[10:]) == int(n):
                    if int(vl[4:6]) == 0:
                        A.append(i)
                        E.append(vl)
        A = np.array(A)
        E = np.array(E)
        return A, E

    def get_selected_ordering_(self):
        lst_resid = self.n_residues
        mapping_ = self.get_mapping_dihedral_angles()
        lst = mapping_[1]
        x = []
        for n, resid in enumerate(lst_resid):
            if resid < self.residues+1:
                for i, vl in enumerate(lst):
                    if int(vl[4:6]) == 0 and int(vl[10:]) == int(resid):
                        x.append(i)
        x = np.array(x)
        return x

    def get_full_ordered_list_dihedral_angles(self):
        return self.get_mapping_dihedral_angles()[1]

    def get_dihedral_atoms_indexes(self):
        F = self.get_mapping_dihedral_angles()[0]
        full = self.feature()
        indexes = np.concatenate(
            [full[0].angle_indexes, full[1].angle_indexes])
        f_ind = indexes[F]
        return f_ind

    def get_angles(self):
        if self.n_residues == 'all':
            f = coor.featurizer(self.pdb)
            f.add_dihedrals(self.get_dihedral_atoms_indexes())
            data = coor.source(self.trajs, features=f)
        elif self.n_residues != 'all':
            f = coor.featurizer(self.pdb)
            f.add_dihedrals(self.get_dihedral_atoms_indexes()
                            [self.get_selected_ordering_()])
            data = coor.source(self.trajs, features=f)
        return data

    def write_data(self, name):
        self.get_angles().write_to_hdf5(name)

    def get_data(self):
        return self.get_angles().get_output()

    def write_resdues_h5(self, path):
        top = mda.Universe(self.pdb)
        dt = np.concatenate(self.get_data())
        if self.n_residues == 'all':
            ls_angles = self.get_full_ordered_list_dihedral_angles()
            with h5.File(f'{path}/dihedral.h5', 'a') as f:
                for n, num in enumerate(range(1, len(top.residues)+1)):
                    resid = []
                    for i, vl in enumerate(ls_angles):
                        if int(vl[10:]) == int(n+1):
                            resid.append(dt[:, i])
                    resid = np.array(resid)
                    f.create_dataset(
                        f'resid{n+1}', data=resid, dtype=float)

def js(P, Q):
    _P = P
    _Q = Q
    _M = 0.5 * (_P + _Q)
    return entropy(_M) - 0.5*(entropy(_P)+entropy(_Q))

def norm_const(f):
    N, error = integrate.quad(lambda x: f(x), -np.pi, np.pi)
    return N

def JSD_gaussian(x1, x2, weight1, weight2): #For data of torsion angles determined within -pi to pi
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    if weight1 != None and weight2 != None:
        new_w1 = np.concatenate([weight1]*3)
        new_w2 = np.concatenate([weight2]*3)
    else:
        new_w1 = None
        new_w2 = None
    kde_x1 = gaussian_kde(X1, weights=new_w1, bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights=new_w2, bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 10000)
    N1 = norm_const(kde_x1)
    N2 = norm_const(kde_x2)
    print (N1, N2)
    z1 = kde_x1(points)/N1
    z2 = kde_x2(points)/N2
    jsd = js(z1, z2)
    return jsd

def JSD_histogram(x1, x2, w1, w2): #For a faster calculation!
    N1, _ = np.histogram(x1, density=True, bins = 1000, weights=w1)
    N2, _ = np.histogram(x2, density=True, bins = 1000, weights=w2)
    jsd = js(N1, N2)
    print ('Value of JSD is: ', jsd)
    return jsd
