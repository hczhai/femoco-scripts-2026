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


from pyscf import gto, scf, solvent
import numpy as np

def add_spin_dim(H, spin, non_spin_dim=3):
    H = np.asarray(H)
    if H.ndim == non_spin_dim:
        H = H[None]
    assert H.ndim == (non_spin_dim + 1)
    if H.shape[0] < spin:
        H = np.asarray((H[0],) * spin)
    return H

def transform_rdm1_to_ao(dm_mo_mo, c_ao_mo):
    dm_mo_mo = np.asarray(dm_mo_mo)
    c_ao_mo = np.asarray(c_ao_mo)

    nao = c_ao_mo.shape[-2]
    if c_ao_mo.ndim < dm_mo_mo.ndim:
        c_ao_mo = add_spin_dim(c_ao_mo, dm_mo_mo.shape[0], non_spin_dim=2)
    if c_ao_mo.ndim == 2:
        dm_ao_ao = c_ao_mo @ dm_mo_mo @ c_ao_mo.conj().T
    else:
        spin = c_ao_mo.shape[0]
        dm_mo_mo = add_spin_dim(dm_mo_mo, spin, non_spin_dim=2)
        assert dm_mo_mo.ndim == c_ao_mo.ndim
        dm_ao_ao = np.zeros((spin, nao, nao), dtype=c_ao_mo.dtype)
        for s in range(spin):
            dm_ao_ao[s] = c_ao_mo[s] @ dm_mo_mo[s] @ c_ao_mo[s].conj().T
    return dm_ao_ao

def transform_rdm1_to_mo(rdm1_ao_ao, c_ao_mo, s_ao_ao):
    rdm1_ao_ao = np.asarray(rdm1_ao_ao)
    c_ao_mo = np.asarray(c_ao_mo)
    nao = c_ao_mo.shape[-2]
    if c_ao_mo.ndim < rdm1_ao_ao.ndim:
        c_ao_mo = add_spin_dim(c_ao_mo, rdm1_ao_ao.shape[0], non_spin_dim=2)
    if c_ao_mo.ndim == 2:
        c_inv = c_ao_mo.conj().T.dot(s_ao_ao)
        rdm1_mo_mo = c_inv @ rdm1_ao_ao @ c_inv.conj().T
    else:
        spin = c_ao_mo.shape[0]
        rdm1_ao_ao = add_spin_dim(rdm1_ao_ao, spin, non_spin_dim=2)
        assert rdm1_ao_ao.ndim == c_ao_mo.ndim
        rdm1_mo_mo = np.zeros((spin, nao, nao), dtype=c_ao_mo.dtype)
        for s in range(spin):
            c_inv = c_ao_mo[s].conj().T.dot(s_ao_ao)
            rdm1_mo_mo[s] = c_inv @ rdm1_ao_ao[s] @ c_inv.conj().T
    return rdm1_mo_mo

def get_metal_init_guess(mol, orbs=["3d"], atom_idxs=[0], coupling="+", atomic_spins=[4], nelec_corr=None):
    from pyscf import scf, lo
    from pyscf.lib import logger

    rmf = scf.RHF(mol)
    dm0 = rmf.get_init_guess(key="atom")
    if dm0.ndim == 2:
        dm0 = np.array([dm0, dm0]) / 2.0
    nelec = np.einsum("sij,ji->s", dm0, rmf.get_ovlp())
    logger.info(scf.RHF(mol), "NELEC BEFORE = %.1f %.1f", nelec[0], nelec[1])

    idxs = []
    assert len(atom_idxs) == len(coupling)
    for ia, orx in zip(atom_idxs, orbs):
        if orx.startswith(mol.atom_symbol(ia) + ' '):
            label = "%d %s.*" % (ia, orx)
        else:
            label = "%d %s %s.*" % (ia, mol.atom_symbol(ia), orx)
        idx = mol.search_ao_label(label)
        idxs.append(idx)
        logger.info(scf.RHF(mol), "ORB = %s IDX = %s", label, idx)

    ld = lo.orth_ao(mol, "lowdin", pre_orth_ao="SCF")
    dl0 = transform_rdm1_to_mo(dm0, ld, rmf.get_ovlp())

    dax, dbx = 0, 0
    if nelec_corr == "all":
        dl0[0] += (mol.nelectron / 2.0 - nelec[0]) / len(dl0[0])
        dl0[1] += (mol.nelectron / 2.0 - nelec[1]) / len(dl0[1])
    elif nelec_corr == "metal":
        dax = (mol.nelectron / 2.0 - nelec[0]) / len(atom_idxs)
        dbx = (mol.nelectron / 2.0 - nelec[1]) / len(atom_idxs)

    for i in range(len(atom_idxs)):
        dspin = atomic_spins[i] / 2.0
        idx = idxs[i]
        if coupling[i] == "+":
            dl0[0][idx, idx] += (dax + dspin) / len(idx)
            dl0[1][idx, idx] += (dbx - dspin) / len(idx)
        else:
            dl0[0][idx, idx] += (dax - dspin) / len(idx)
            dl0[1][idx, idx] += (dbx + dspin) / len(idx)

    dm0 = transform_rdm1_to_ao(dl0, ld)
    nelec = np.einsum("sij,ji->s", dm0, rmf.get_ovlp())
    logger.info(scf.RHF(mol), "NELEC AFTER  = %.1f %.1f", nelec[0], nelec[1])
    return dm0

fid = 2

fm_basis = {'Fe': 'TZP-DKH', 'Mo': 'TZP-DKH', 'S': 'TZP-DKH', 'H': 'def2-SVP'}

mol = gto.M(atom="../MoFe3.xyz", basis=fm_basis, charge=-3, spin=3, symmetry=False, verbose=3)
mol.max_memory = 500000

mf = scf.sfx2c(scf.UKS(mol))
mf.xc = 'b3lyp'
mf = solvent.ddcosmo.ddcosmo_for_scf(mf)
mf.with_solvent.eps = 4.0
mf.max_cycle = 200
mf.verbose = 4

dm0 = get_metal_init_guess(mol, orbs=["4d", "3d", "3d", "3d"], atom_idxs=[0, 1, 2, 3],
    coupling="-+++", atomic_spins=[5, 5, 5, 5], nelec_corr=None)

mf.kernel(dm0=dm0)
dm1 = np.array(mf.make_rdm1())

mf = mf.newton()
mf.max_cycle = 1200
mf.conv_tol = 1E-10
mf.ah_conv_tol = 1E-14
mf.kernel(dm0=dm1)
dm1 = np.array(mf.make_rdm1())

np.save('%02d-mo-coeff.npy' % fid, mf.mo_coeff)
np.save('%02d-mo-occ.npy' % fid, mf.mo_occ)
np.save('%02d-mo-energy.npy' % fid, mf.mo_energy)
np.save('%02d-uhf-dm1.npy' % fid, dm1)
