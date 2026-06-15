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

eig = np.load('17-uno-occ.npy')
coeff = np.load('17-uno-coeff.npy')

print(len(eig[(eig <= 1.9700) & (eig >= 0.0200)]), 35 + 35 + 47, sum(eig[(eig <= 1.9700) & (eig >= 0.0200)]), 35 + 35 * 2)
print(len(eig[(eig <= 1.9840) & (eig >= 0.0100)]), 35 + 54 + 91, sum(eig[(eig <= 1.9840) & (eig >= 0.0100)]), 35 + 54 * 2)
print(len(eig[(eig <= 1.9960) & (eig >= 0.0050)]), 35 + 91 + 159, sum(eig[(eig <= 1.9960) & (eig >= 0.0050)]), 35 + 91 * 2)
print(len(eig[(eig <= 1.9985) & (eig >= 0.0025)]), 35 + 121 + 248, sum(eig[(eig <= 1.9985) & (eig >= 0.0025)]), 35 + 121 * 2)

pfid = 12

ecore = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/ecore.npy' % pfid)
h1e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/h1e.npy' % pfid)
g2e = np.load('/mnt/ceph/users/hzhai/%02d-fm-active/g2e.npy' % pfid)

mo1 = coeff[:, (eig <= 1.9700) & (eig >= 0.0200)]
mo2 = coeff[:, (eig <= 1.9840) & (eig >= 0.0100)]
mo3 = coeff[:, (eig <= 1.9960) & (eig >= 0.0050)]
mo4 = coeff[:, (eig <= 1.9985) & (eig >= 0.0025)]

moc1 = coeff[:, (eig > 1.9700)]
moc2 = coeff[:, (eig > 1.9840)]
moc3 = coeff[:, (eig > 1.9960)]
moc4 = coeff[:, (eig > 1.9985)]

from pyscf import ao2mo

fid = 18

def unpack(g2e, ml, mr):
    return g2e.reshape((ml, ml, mr, mr))

for im, (mo, moc) in enumerate(zip([mo1, mo2, mo3, mo4], [moc1, moc2, moc3, moc4])):
    print('---', im)
    xg2e = ao2mo.restore(8, ao2mo.general(g2e, (mo, mo, mo, mo)), mo.shape[1])
    np.save('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-g2e.npy' % (fid, im), xg2e)
    xg2e = None
    xh1e = np.einsum('ij,ia,jb->ab', h1e, mo, mo, optimize='optimal')
    upk = unpack(ao2mo.general(g2e, (mo, mo, moc, moc), compact=False), mo.shape[1], moc.shape[1])
    xh1e += 2 * np.einsum('pqii->pq', upk, optimize='optimal')
    upk = None
    xh1e -= np.einsum('piiq->pq', ao2mo.general(g2e, (mo, moc, moc, mo), compact=False)
        .reshape((mo.shape[1], moc.shape[1], moc.shape[1], mo.shape[1])), optimize='optimal')
    xcore = ecore + np.einsum("ij,ij->", 2 * moc @ moc.T, h1e, optimize='optimal')
    cc = ao2mo.restore(1, ao2mo.general(g2e, (moc, moc, moc, moc), compact=False), moc.shape[1])
    xcore += 2 * np.einsum('iijj->', cc, optimize='optimal')
    xcore -= np.einsum('ijji->', cc, optimize='optimal')

    np.save('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-ecore.npy' % (fid, im), xcore)
    np.save('/mnt/ceph/users/hzhai/%02d-fm-active/%02d-h1e.npy' % (fid, im), xh1e)
