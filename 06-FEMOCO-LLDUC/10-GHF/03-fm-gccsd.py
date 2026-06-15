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

fid = 3
iiid = int(sys.argv[1])

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

xstrs = open('01-fm-ghf-all.out', 'r').readlines()

def gccsd_t(cc):
    import pytblis
    t1, t2 = cc.t1, cc.t2
    eris = cc.ao2mo(cc.mo_coeff)
    nocc, nvir = t1.shape

    bcei = np.asarray(eris.ovvv).conj().transpose(3,2,1,0)
    majk = np.asarray(eris.ooov).conj().transpose(2,3,0,1)
    bcjk = np.asarray(eris.oovv).conj().transpose(2,3,0,1)
    fvo = eris.fock[nocc:,:nocc]
    mo_e = eris.mo_energy
    eijk = lib.direct_sum('i+j+k->ijk', mo_e[:nocc], mo_e[:nocc], mo_e[:nocc])
    eabc = lib.direct_sum('a+b+c->abc', mo_e[nocc:], mo_e[nocc:], mo_e[nocc:])

    t2T = t2.transpose(2,3,0,1)
    t1T = t1.T
    def get_wv(a, b, c):
        w  = pytblis.contract('ejk,ei->ijk', t2T[a,:], bcei[b,c], alpha=1.0)
        pytblis.contract('im,mjk->ijk', t2T[b,c], majk[:,a], alpha=-1.0, beta=1.0, out=w)
        v = w.copy()
        pytblis.contract('i,jk->ijk', t1T[a], bcjk[b,c], alpha=1.0, beta=1.0, out=v)
        pytblis.contract('i,jk->ijk', fvo[a], t2T [b,c], alpha=1.0, beta=1.0, out=v)
        w = w + w.transpose(2,0,1) + w.transpose(1,2,0)
        return w, v

    et = 0
    for a in range(nvir):
        for b in range(a):
            for c in range(b):
                wabc, vabc = get_wv(a, b, c)
                wcab, vcab = get_wv(c, a, b)
                wbac, vbac = get_wv(b, a, c)

                w = wabc + wcab - wbac
                v = vabc + vcab - vbac
                w /= eijk - eabc[a,b,c]
                et += np.dot(w.ravel(), v.conj().ravel())
    et /= 2
    return et

for ix in range(iiid * 4, (iiid + 1) * 4):
    if ix >= len(xstrs):
        break
    print('---- IX = ', ix)
    xdstr, ref_e = xstrs[ix].split()[3], float(xstrs[ix].split()[4])
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

    from pyscf import cc

    mmc = cc.GCCSD(gmf)
    mmc.diis_space = 36
    mmc.max_cycle = 400
    mmc.kernel()

    e_ccsd_t = mmc.e_tot + gccsd_t(mmc)
    print("%5d :: %5d %s %20.12f %20.12f %d" % (iiix, iiiy, xdstr, mmc.e_tot, e_ccsd_t, int(mmc.converged)))
