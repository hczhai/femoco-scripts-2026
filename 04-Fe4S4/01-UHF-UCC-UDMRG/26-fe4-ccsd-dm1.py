
import numpy as np, sys

fid = 26
pfid = 21
iiid = int(sys.argv[1])

def charge(ixs, jxs, dm1, dm2):
    dm1a, dm1b = dm1
    r = 0
    for i, j in zip(ixs, jxs):
        r += dm1a[i, j] + dm1b[i, j]
    return r

fe_idxs = [[2, 3, 4, 5, 6], [7, 8, 9, 10, 11], [24, 25, 26, 27, 28], [29, 30, 31, 32, 33]]
order = 2

uref_dm1 = np.load('data/%02d-cc-%02d-%02d-pdm1.npy' % (pfid, order, iiid))
uref_dm2 = np.load('data/%02d-cc-%02d-%02d-pdm2.npy' % (pfid, order, iiid))
mo_coeff = np.load("data/%02d-uhf-mo-coeff-%02d.npy" % (pfid, iiid))
dm1 = np.einsum('xpi,xij,xqj->xpq', mo_coeff, uref_dm1, mo_coeff, optimize=True)
dm2aa = np.einsum('pi,qj,ijkl,rk,sl->pqrs', mo_coeff[0], mo_coeff[0], uref_dm2[0], mo_coeff[0], mo_coeff[0], optimize=True)
dm2ab = np.einsum('pi,qj,ijkl,rk,sl->pqrs', mo_coeff[0], mo_coeff[0], uref_dm2[1], mo_coeff[1], mo_coeff[1], optimize=True)
dm2bb = np.einsum('pi,qj,ijkl,rk,sl->pqrs', mo_coeff[1], mo_coeff[1], uref_dm2[2], mo_coeff[1], mo_coeff[1], optimize=True)
dm2 = np.array([dm2aa, dm2ab, dm2bb])

scr = np.zeros((4, 4))
for ii, ix in enumerate(fe_idxs):
    for jj, jx in enumerate(fe_idxs):
        scr[ii, jj] = charge(ix, jx, dm1, dm2)
    print('CHARGE %d ::' % ii, ' '.join(['%12.6f' % x for x in scr[ii]]))
