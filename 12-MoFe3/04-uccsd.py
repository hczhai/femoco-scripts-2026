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


from pyscf import gto, scf, solvent, lib
import numpy as np
from pyscf.data.elements import chemcore

fid = 4
pfid = 2

fm_basis = {'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH', 'H': 'def2-SVP'}

mol = gto.M(atom="../MoFe3.xyz", basis=fm_basis, charge=-3, spin=3, symmetry=False, verbose=3)
mol.max_memory = 1400000

frozen = chemcore(mol)
print('frozen =', frozen)

mf = scf.sfx2c(scf.UHF(mol))
mf = solvent.ddcosmo.ddcosmo_for_scf(mf)
mf.with_solvent.eps = 4.0
mf.max_cycle = 200
mf.verbose = 4

dm0 = np.load('%02d-uhf-dm1.npy' % pfid)

mf = mf.newton()
mf.max_cycle = 1200
mf.conv_tol = 1E-10
mf.ah_conv_tol = 1E-14
mf.kernel(dm0=dm0)
dm1 = np.array(mf.make_rdm1())

np.save('%02d-mo-coeff.npy' % fid, mf.mo_coeff)
np.save('%02d-mo-occ.npy' % fid, mf.mo_occ)
np.save('%02d-mo-energy.npy' % fid, mf.mo_energy)
np.save('%02d-uhf-dm1.npy' % fid, dm1)

from pyscf.cc.uccsd import UCCSD, _ChemistsERIs
from pyscf import ao2mo

def _make_eris_incore(mycc, mo_coeff=None, ao2mofn=None):
    eris = _ChemistsERIs()
    eris._common_init_(mycc, mo_coeff)

    nocca, noccb = mycc.nocc
    nmoa, nmob = mycc.nmo
    nvira, nvirb = nmoa-nocca, nmob-noccb

    moa = eris.mo_coeff[0]
    mob = eris.mo_coeff[1]
    nmoa = moa.shape[1]
    nmob = mob.shape[1]

    if callable(ao2mofn):
        eri_aa = ao2mofn(moa).reshape([nmoa]*4)
    else:
        eri_aa = ao2mo.restore(1, ao2mo.full(mycc._scf._eri, moa), nmoa)
    eri_aa = eri_aa.reshape(nmoa,nmoa,nmoa,nmoa)

    eris.oooo = eri_aa[:nocca,:nocca,:nocca,:nocca].copy()
    eris.ovoo = eri_aa[:nocca,nocca:,:nocca,:nocca].copy()
    eris.ovov = eri_aa[:nocca,nocca:,:nocca,nocca:].copy()
    eris.oovv = eri_aa[:nocca,:nocca,nocca:,nocca:].copy()
    eris.ovvo = eri_aa[:nocca,nocca:,nocca:,:nocca].copy()
    eris.ovvv = eri_aa[:nocca,nocca:,nocca:,nocca:].copy()
    eris.vvvv = eri_aa[nocca:,nocca:,nocca:,nocca:].copy()
    eri_aa = None

    if callable(ao2mofn):
        eri_bb = ao2mofn(mob).reshape([nmob]*4)
    else:
        eri_bb = ao2mo.restore(1, ao2mo.full(mycc._scf._eri, mob), nmob)
    eri_bb = eri_bb.reshape(nmob,nmob,nmob,nmob)

    eris.OOOO = eri_bb[:noccb,:noccb,:noccb,:noccb].copy()
    eris.OVOO = eri_bb[:noccb,noccb:,:noccb,:noccb].copy()
    eris.OVOV = eri_bb[:noccb,noccb:,:noccb,noccb:].copy()
    eris.OOVV = eri_bb[:noccb,:noccb,noccb:,noccb:].copy()
    eris.OVVO = eri_bb[:noccb,noccb:,noccb:,:noccb].copy()
    eris.OVVV = eri_bb[:noccb,noccb:,noccb:,noccb:].copy()
    eris.VVVV = eri_bb[noccb:,noccb:,noccb:,noccb:].copy()
    eri_bb = None

    if callable(ao2mofn):
        eri_ab = ao2mofn((moa,moa,mob,mob))
    else:
        eri_ab = ao2mo.general(mycc._scf._eri, (moa,moa,mob,mob), compact=False)
    eri_ab = eri_ab.reshape(nmoa,nmoa,nmob,nmob)
    eri_ba = eri_ab.transpose(2,3,0,1)

    eris.ooOO = eri_ab[:nocca,:nocca,:noccb,:noccb].copy()
    eris.ovOO = eri_ab[:nocca,nocca:,:noccb,:noccb].copy()
    eris.ovOV = eri_ab[:nocca,nocca:,:noccb,noccb:].copy()
    eris.ooVV = eri_ab[:nocca,:nocca,noccb:,noccb:].copy()
    eris.ovVO = eri_ab[:nocca,nocca:,noccb:,:noccb].copy()
    eris.ovVV = eri_ab[:nocca,nocca:,noccb:,noccb:].copy()
    eris.vvVV = eri_ab[nocca:,nocca:,noccb:,noccb:].copy()

    #eris.OOoo = eri_ba[:noccb,:noccb,:nocca,:nocca].copy()
    eris.OVoo = eri_ba[:noccb,noccb:,:nocca,:nocca].copy()
    #eris.OVov = eri_ba[:noccb,noccb:,:nocca,nocca:].copy()
    eris.OOvv = eri_ba[:noccb,:noccb,nocca:,nocca:].copy()
    eris.OVvo = eri_ba[:noccb,noccb:,nocca:,:nocca].copy()
    eris.OVvv = eri_ba[:noccb,noccb:,nocca:,nocca:].copy()
    #eris.VVvv = eri_ba[noccb:,noccb:,nocca:,nocca:].copy()
    eri_ba = None
    eri_ab = None

    if not callable(ao2mofn):
        ovvv = eris.ovvv.reshape(nocca*nvira,nvira,nvira)
        eris.ovvv = lib.pack_tril(ovvv).reshape(nocca,nvira,nvira*(nvira+1)//2)
        eris.vvvv = ao2mo.restore(4, eris.vvvv, nvira)

        OVVV = eris.OVVV.reshape(noccb*nvirb,nvirb,nvirb)
        eris.OVVV = lib.pack_tril(OVVV).reshape(noccb,nvirb,nvirb*(nvirb+1)//2)
        eris.VVVV = ao2mo.restore(4, eris.VVVV, nvirb)

        ovVV = eris.ovVV.reshape(nocca*nvira,nvirb,nvirb)
        eris.ovVV = lib.pack_tril(ovVV).reshape(nocca,nvira,nvirb*(nvirb+1)//2)
        vvVV = eris.vvVV.reshape(nvira**2,nvirb**2)
        idxa = np.tril_indices(nvira)
        idxb = np.tril_indices(nvirb)
        eris.vvVV = lib.take_2d(vvVV, idxa[0]*nvira+idxa[1], idxb[0]*nvirb+idxb[1])

        OVvv = eris.OVvv.reshape(noccb*nvirb,nvira,nvira)
        eris.OVvv = lib.pack_tril(OVvv).reshape(noccb,nvirb,nvira*(nvira+1)//2)
    return eris

UCCSD.ao2mo = lambda self, mo_coeff=None: _make_eris_incore(self, mo_coeff)

mmc = UCCSD(mf, frozen=frozen)
mmc.verbose = 5
mmc.diis_space = 12
mmc.max_cycle = 800
mmc.incore_complete = True
mmc.kernel()

dm = mmc.make_rdm1()

np.save('%02d-uccsd-dm1.npy' % fid, dm)
