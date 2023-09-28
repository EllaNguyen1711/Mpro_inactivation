import pyemma
import pyemma.coordinates as coor
import h5py as h5
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
from scipy.stats import entropy
from numpy.linalg import norm


class Dihedral:

    def __init__(self, pdb, trajs, n_residues = 'all'):
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
            if resid < 306:
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
                    elif res_ind < 306:
                        result = []
                        for i, vl in enumerate(ls_angles):
                            if int(vl[4:6]) == 0 and int(vl[10:]) == int(res_ind):
                                result.append(dt[:, i])
                        result = np.array(result)
                        f.create_dataset(
                            f'resid{res_ind}', data=result, dtype=float)


def js(P, Q):
    _P = P / norm(P, ord=1)
    _Q = Q / norm(Q, ord=1)
    _M = 0.5 * (_P + _Q)
    return entropy(_M) - 0.5*(entropy(_P)+entropy(_Q))


def JSD(x1, x2, weight1, weight2):
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    kde_x1 = gaussian_kde(X1, weights=normalize(
        np.concatenate([weight1]*3)), bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights=normalize(
        np.concatenate([weight2]*3)), bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 10000)
    z1 = kde_x1(points)*3
    z2 = kde_x2(points)*3
    jsd = js(z1, z2)
    return jsd


def weight_samples(ls_clusters, f_probs, cluster_indexes, path):
    cl_probs = f_probs[ls_clusters]
    x = []
    for j, vl in enumerate(ls_clusters):
        cl_ind = cluster_indexes[vl]
        traj = path + '/cluster_%s.xtc' % vl
        u = mda.Universe(traj)
        x.append([cl_probs[j]/len(cl_ind)]*len(u.trajectory))
    x = np.concatenate(x)
    return x


def weight_full_traj(dtrajs, cluster_indexes, f_probs):
    probs_fr = []
    for i, vl1 in enumerate(f_probs):
        probs_fr.append(vl1/len(cluster_indexes[i]))
    traj_probs = []
    for j, vl2 in enumerate(np.array(dtrajs)):
        traj = []
        for l, vl3 in enumerate(dtrajs[vl2]):
            traj.append(probs_fr[int(vl3)])
        traj_probs.append(np.array(traj))
    traj_probs = np.array(traj_probs)
    return traj_probs


def KDE_plot_(data, st1, st2, weight1, weight2, residue, angle, path_to_f):
    list_torsion = ['phi', 'psi', 'chi1', 'chi2']
    x1 = data[st1][angle]
    x2 = data[st2][angle]
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    kde_x1 = gaussian_kde(X1, weights=normalize(
        np.concatenate([weight1]*3)), bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights=normalize(
        np.concatenate([weight2]*3)), bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 10000)  # defined probability space X
    z1 = kde_x1(points)*3
    z2 = kde_x2(points)*3
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    ax.plot(points, z1, label='State%s' % st1, color='deepskyblue')
    ax.plot(points, z2, label='State%s' % st2, color='orange')
    ax.legend()
    ax.set_xlabel(r'$\%s$' % list_torsion[angle])
    ax.set_ylabel('Probability distribitions')
    ax.set_title(r'Residue%s_PDF of angle $\%s$ over states %s and %s' % (residue, list_torsion[angle],
                                                                          st1, st2))
    plt.savefig(
        f'{path_to_f}/Residue%s_PDF-of-angle-%s-over-states-%s-and-%s.png' % (residue, list_torsion[angle], st1, st2))


def KDE_plot(data, st1, st2, weight1, weight2, residue, angle, peaks, path_to_f):
    list_torsion = ['phi', 'psi', 'chi1', 'chi2']
    x1 = data[st1][angle]
    x2 = data[st2][angle]
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    kde_x1 = gaussian_kde(X1, weights=normalize(
        np.concatenate([weight1]*3)), bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights=normalize(
        np.concatenate([weight2]*3)), bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 10000)  # defined probability space X
    z1 = kde_x1(points)*3
    z2 = kde_x2(points)*3
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(points, z1, label='State%s' % st1, color='deepskyblue')
    ax.plot(points, z2, label='State%s' % st2, color='orange')
    if st1 == 0:
        ax.vlines(x=x1[peaks], ymin=0, ymax=1.4,
                  color='black', linestyle='--', linewidth=2.0)
    if st1 == 3:
        ax.vlines(x=x2[peaks], ymin=0, ymax=1.4,
                  color='black', linestyle='--', linewidth=2.0)
    ax.legend()
    ax.set_xlabel(r'$\%s$' % list_torsion[angle])
    ax.set_ylabel('Probability distribitions')
    ax.set_title(r'Residue%s_PDF of angle $\%s$ over states %s and %s' % (residue, list_torsion[angle],
                                                                          st1, st2))
    fig.tight_layout()
    plt.savefig(
        f'{path_to_f}/Residue%s_PDF-of-angle-%s-over-states-%s-and-%s.png' % (residue, list_torsion[angle], st1, st2))


def KDE_plot_residues(data, resid1, resid2, weight, angle, path_to_f):
    list_torsion = ['phi', 'psi', 'chi1', 'chi2']
    x1 = data[f'resid{resid1}'][angle]
    x2 = data[f'resid{resid2}'][angle]
    X1 = np.concatenate([x1-2*np.pi, x1, x1+2*np.pi])
    X2 = np.concatenate([x2-2*np.pi, x2, x2+2*np.pi])
    kde_x1 = gaussian_kde(X1, weights= np.concatenate([weight]*3), bw_method=.05)
    kde_x2 = gaussian_kde(X2, weights= np.concatenate([weight]*3), bw_method=.05)
    points = np.linspace(-np.pi, np.pi, 1000)  # defined probability space X
    z1 = kde_x1(points)*3
    z2 = kde_x2(points)*3
    fig, ax = plt.subplots(figsize=(8, 6))
    fig.tight_layout()
    fig.subplots_adjust(top=0.9)
    ax.plot(points, z1, label='A%s' % (resid1), color='deepskyblue')
    ax.plot(points, z2, label='B%s' % (resid2-306), color='orange')
    ax.legend()
    ax.set_xlabel(r'$\%s$' % list_torsion[angle])
    ax.set_ylabel('Probability distribitions')
    ax.set_title(r'A%s$\%s$ and B%s$\%s$' % (resid1, list_torsion[angle], resid2-306, list_torsion[angle]))
    plt.savefig(
        f'{path_to_f}/Residue%s_PDF_of_%s.png' % (resid1, list_torsion[angle]))

