# ---------------------------------------------------------------------------
# Copyright (C) 2026 Huanchen Zhai. All Rights Reserved.
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ---------------------------------------------------------------------------


import zipfile
import itertools
import numpy as np, sys, tempfile
from pyscf import lib
from pyscf.tools.fcidump import to_scf
import scipy.linalg

fid = 5

lib.param.TMPDIR = '/tmp'

tmpdir = tempfile.mkdtemp()

with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
    zip_ref.extract('FEMOCO.FCIDUMP', tmpdir)

mf = to_scf(tmpdir + '/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

gmf = mf.to_ghf()

n_sites = mf.mol.nao
h1e = mf.get_hcore()
ovlp = mf.get_ovlp()

gmf.get_hcore = lambda *args: scipy.linalg.block_diag(h1e, h1e)
gmf.get_ovlp = lambda *args: scipy.linalg.block_diag(ovlp, ovlp)

gmf = gmf.newton()

xstrs = open('04-fm-gccsd-sorted.out', 'r').readlines()
hf_xstrs = open('01-fm-ghf-all.out', 'r').readlines()

e_hfs = {l.split()[3]: float(l.split()[4]) for l in hf_xstrs}

fe_idxs = [
    [2, 3, 4, 5, 6], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25],
    [26, 27, 28, 29, 30], [44, 45, 46, 47, 48], [49, 50, 51, 52, 53],
    [54, 55, 56, 57, 58]
]
mo_idxs = [[65, 69, 70, 71, 72]]

all_idxs = [np.array(x) for x in fe_idxs + mo_idxs]

for ix in range(16):
    print('---- IX = ', ix)
    xdstr = xstrs[ix].split()[3]
    ref_e = e_hfs[xdstr]
    iiix, iiiy = int(xstrs[ix].split()[0]), int(xstrs[ix].split()[2])

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    gdm = scipy.linalg.block_diag(dm0[0], dm0[1])

    np.random.seed(2025)
    gdm += (np.random.random(gdm.shape) - 0.5) * 0.1

    gmf.verbose = 4
    gmf.max_cycle = 800
    gmf.conv_tol = 1E-13
    gmf.ah_conv_tol = 1E-16
    gmf.kernel(dm0=gdm)

    if abs(gmf.e_tot - ref_e) >= 1E-8:
        print('!!WARNING', gmf.e_tot, ref_e)
    
    gdm = gmf.make_rdm1()

    vecs = []

    for idxs in all_idxs:
        pauli = 0.5 * lib.PauliMatrices
        xov = ovlp[idxs, :][:, idxs]
        gidxs = np.array([x + p for p in [0, len(ovlp)] for x in idxs])
        dms = np.einsum('ij,xst->xsitj', xov, pauli).reshape(3, len(xov) * 2, len(xov) * 2).conj()
        dms = (dms + dms.transpose(0,2,1).conj()) / 2
        xgdm = gdm[gidxs, :][:, gidxs]
        svec = np.einsum('ij,xij->x', xgdm, dms)
        assert np.linalg.norm(np.imag(svec)) < 1E-10
        print(idxs, svec.real, np.linalg.norm(svec.real))
        vecs.append(svec.real)

    def rotate(vs, target):
        from scipy.spatial.transform import Rotation as R
        p1_norm = vs[0] / np.linalg.norm(vs[0])
        p1_target_norm = target / np.linalg.norm(target)
        rotation_axis = np.cross(p1_norm, p1_target_norm)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        rotation_angle = np.arccos(np.dot(p1_norm, p1_target_norm))
        rotation_vector = rotation_axis * rotation_angle
        rotation = R.from_rotvec(rotation_vector)
        new_vs = [rotation.apply(x) for x in vs]
        return np.array(new_vs)

    print(" ".join(["%.3f" % np.linalg.norm(v) for v in vecs]))

    vecs = rotate(vecs, [-0.9, -0.2, 0.5])
    vecs[-1] *= 2.5

    temp = open('05-fm-template.vesta', 'r').readlines()
    istart = temp.index("VECTR\n")
    iend = temp.index("SPLAN\n")

    with open('05-fm-ghf-spin-%d.vesta' % ix, 'w') as f:
        f.writelines(temp[:istart])
        f.write("VECTR\n")
        for i in range(len(vecs)):
            f.write("%4d %10.5f %10.5f %10.5f 0\n" % (i + 1, *(vecs[i] / 2)))
            f.write("%5d %4d %4d %4d %4d\n" % (i + 1, 0, 0, 0, 0))
            f.write(" 0 0 0 0 0\n")
        f.write(" 0 0 0 0 0\n")
        f.write("VECTT\n")
        for i in range(len(vecs)):
            f.write("%4d  0.250   0 255   0 1\n" % (i + 1))
        f.write(" 0 0 0 0 0\n")
        f.writelines(temp[iend:])
