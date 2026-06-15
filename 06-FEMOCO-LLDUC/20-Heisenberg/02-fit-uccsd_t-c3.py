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


import numpy as np, zipfile
from pyscf.tools.fcidump import to_scf
from pyscf import ao2mo
import tempfile

with tempfile.TemporaryDirectory() as tmpdir:
    with zipfile.ZipFile('../FEMOCO.FCIDUMP.zip', 'r') as zip_ref:
        zip_ref.extract('FEMOCO.FCIDUMP', tmpdir)
    mf = to_scf(tmpdir + '/FEMOCO.FCIDUMP', molpro_orbsym=True, mf=None)

intor_symmetric = mf.mol.intor_symmetric
mf.mol.intor_symmetric = lambda intor, **kwargs: np.eye(mf.mol.nao) \
    if intor == 'int1e_ovlp' else intor_symmetric(intor, **kwargs)

mf.mol.groupname = 'c1'
mf.mol.symmetry = True
mf.mol.verbose = 4

umf = mf.to_uhf().newton()

femo_idxs = [
    [2, 3, 4, 5, 6], [16, 17, 18, 19, 20], [21, 22, 23, 24, 25],
    [26, 27, 28, 29, 30], [44, 45, 46, 47, 48], [49, 50, 51, 52, 53],
    [54, 55, 56, 57, 58], [65, 69, 70, 71, 72]
]

conn = [
    [[0, 1], [0, 2], [0, 3]],
    [[1, 2], [1, 3], [2, 3]],
    [[4, 5], [4, 6], [5, 6]],
    [[4, 7], [5, 7], [6, 7]],
    [[1, 4], [2, 5], [3, 6]],
    [[2, 4], [3, 4], [1, 5], [3, 5], [1, 6], [2, 6]],
]

xa = np.zeros((35, 7))
xb = np.zeros((35,))
n_sites = umf.mol.nao

for iiid in range(35):
    l = open("../01-UCCSD/07-fm-uccsd-min.out", "r").readlines()[iiid].split()
    xdstr, e_ccsd_t = l[3], float(l[5])
    print(iiid, '::', xdstr, e_ccsd_t)

    dm0 = np.zeros((2, n_sites, n_sites))
    occ = np.array([{'2': [1.0, 1.0], 'a': [1.0, 0.0], 'b': [0.0, 1.0], '0': [0.0, 0.0]}[x] for x in xdstr])
    dm0[0, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 0]
    dm0[1, np.mgrid[:n_sites], np.mgrid[:n_sites]] = occ[:, 1]

    umf.verbose = 4
    umf.max_cycle = 200
    umf.conv_tol = 1E-12
    umf.ah_conv_tol = 1E-16
    umf.kernel(dm0=dm0)

    print('EHF[%d] =' % iiid, umf.e_tot)
    dm1a, dm1b = umf.make_rdm1()
    dm2aa, dm2ab, dm2bb = umf.make_rdm2()
    g2e = ao2mo.restore(1, mf._eri, mf.mol.nao)

    pdm_e = np.einsum('ij,ij->', dm1a + dm1b, mf.get_hcore()) \
        + 0.5 * np.einsum('ijkl,ijkl->', dm2aa + dm2bb + dm2ab + dm2ab.transpose((2, 3, 0, 1)), g2e) \
        + mf.mol.energy_nuc()
    print('EPDM[%d] =' % iiid, pdm_e)
    assert abs(pdm_e - umf.e_tot) < 1E-8

    for ii, cxx in enumerate(conn):
        q = 0
        for ix, jx in cxx:
            p = 0
            for i in femo_idxs[ix]:
                for j in femo_idxs[jx]:
                    p += 0.25 * dm2aa[i, i, j, j]
                    p -= 0.25 * dm2ab[i, i, j, j]
                    p -= 0.25 * dm2ab[j, j, i, i]
                    p += 0.25 * dm2bb[i, i, j, j]
                    p -= 0.50 * dm2ab[i, j, j, i]
                    p -= 0.50 * dm2ab[j, i, i, j]
            if jx != 7:
                q += 4.00 if p > 0 else -6.00
            else:
                q += 1.00 if p > 0 else -1.50
        xa[iiid, ii] = q
    xa[iiid, -1] = 1

    xb[iiid] = e_ccsd_t

x, res, rk, s = np.linalg.lstsq(xa, xb, rcond=None)
print(x[-1], x[:-1] * 219474.6)
print(res, 'RMSE = ', res ** 0.5)
print(rk)
print(s)

coeffs = [x[-1]]
for ix, cn in zip(x[:-1], conn):
    for p in range(len(cn)):
        coeffs.append(ix)

np.save('02-fit-uccsd_t-c3-mat.npy', xa)
np.save('02-fit-uccsd_t-c3.npy', coeffs)

conn_nm = ["FeTerm-FeL", "FeL-FeL", "FeR-FeR", "FeR-Mo", "FeL-FeR (axial)", "FeL-FeR (cross)"]

xbpred = xa @ x

for ix, v in enumerate(xbpred.tolist()):
    print('%5d' % ix, '::', v)

print('CONST =', coeffs[0])
jx = 1
for nm, v in zip(conn_nm, conn):
    print('NM = %s' % nm, coeffs[jx] * 219474.6)
    jx += len(v)

print('RMSE =', np.linalg.norm(xbpred - xb))
