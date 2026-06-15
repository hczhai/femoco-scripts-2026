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


import numpy as np

ncas, n_elec, spin = 994, 423, 35
nocca = (n_elec + 35) // 2
noccb = (n_elec - 35) // 2

dm1 = np.load('14-ump2-dm1.npy')
mo = np.load('14-mo-coeff.npy')
print(dm1.shape)
print(mo.shape)
print(nocca, noccb)

dm = np.einsum('xpi,xij,xqj->pq', mo, dm1, mo, optimize=True)

eig, coeff = np.linalg.eigh(dm)
idx = np.argsort(eig)[::-1]
eig = eig[idx]
coeff = coeff[:, idx]

for ii, x in enumerate(eig):
    print(ii, x)

print(len(eig[(eig <= 1.9700) & (eig >= 0.0200)]), 35 + 35 + 47)
print(len(eig[(eig <= 1.9840) & (eig >= 0.0100)]), 35 + 54 + 91)
print(len(eig[(eig <= 1.9960) & (eig >= 0.0050)]), 35 + 91 + 158)
print(len(eig[(eig <= 1.9985) & (eig >= 0.0025)]), 35 + 121 + 247)

np.save('17-uno-occ.npy', eig)
np.save('17-uno-coeff.npy', coeff)

pfid = 5

xcmo = np.load('%02d-sel-core-mo-coeff.npy' % pfid)
xamo = np.load('%02d-sel-mo-coeff.npy' % pfid)
xvmo = np.load('%02d-sel-virt-mo-coeff.npy' % pfid)

print(coeff.shape)
print(xcmo.shape)
print(xamo.shape)
print(xvmo.shape)

ump2_coeff = np.zeros((xcmo.shape[0], xcmo.shape[0]))
ump2_coeff[:, :xcmo.shape[1]] = xcmo
ump2_coeff[:, xcmo.shape[1]:xcmo.shape[1] + xamo.shape[1]] = xamo @ coeff
ump2_coeff[:, xcmo.shape[1] + xamo.shape[1]:] = xvmo

np.save('17-mo-coeff.npy', ump2_coeff)

ovlp = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-full-overlap.npy' % (18, 10))
orig_mo_coeff = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-orig-mo-coeff.npy' % (18, 11))

mo_ovlp = np.zeros((xamo.shape[1], ))

for ii in range(xamo.shape[1]):
    v = np.einsum('a,ab,bj->j', ump2_coeff[:, ii + xcmo.shape[1]], ovlp, orig_mo_coeff[:, 517:517 + 76], optimize=True)
    print(ii, np.linalg.norm(v))
    mo_ovlp[ii] = np.linalg.norm(v)

np.save('17-mo-ovlp.npy', mo_ovlp)
