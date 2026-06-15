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
from pyblock2.driver.core import DMRGDriver, SymmetryTypes

local_spin = 4
n_sites = 8
total_spin = (local_spin * 7 + 1) % 2

conn = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3], [4, 5], [4, 6], [5, 6], [4, 7], [5, 7], [6, 7],
    [1, 4], [2, 5], [3, 6], [2, 4], [3, 4], [1, 5], [3, 5], [1, 6], [2, 6]]

jfit_ccsd_t = np.load('./02-fit-uccsd_t-c3.npy')

topology = [(i, j, k) for (i, j), k in zip(conn, jfit_ccsd_t[1:])]
ecore = jfit_ccsd_t[0] + 22140.0

driver = DMRGDriver(scratch='/tmp', symm_type=SymmetryTypes.SZ | SymmetryTypes.SAny, n_threads=64)

driver.set_symmetry_groups("U1")
Q = driver.bw.SX

# [Part A] Set states and matrix representation of operators in local Hilbert space
site_basis, site_ops = [], []

for k in range(n_sites):
    xspin = local_spin if k < n_sites - 1 else 1
    basis = [(Q(tm), 1) for tm in range(-xspin, xspin + 1, 2)] # [01]
    ident = np.zeros((len(basis), len(basis)))
    sp = np.zeros((len(basis), len(basis)))
    sm = np.zeros((len(basis), len(basis)))
    sz = np.zeros((len(basis), len(basis)))
    for im, tm in enumerate(range(-xspin, xspin + 1, 2)):
        ident[im, im] = 1.0
    for im, tm in list(enumerate(range(-xspin, xspin + 1, 2)))[:-1]:
        sp[im + 1, im] = np.sqrt((xspin - tm) * (xspin + tm + 2) // 4)
    for im, tm in list(enumerate(range(-xspin, xspin + 1, 2)))[1:]:
        sm[im - 1, im] = np.sqrt((xspin + tm) * (xspin - tm + 2) // 4)
    for im, tm in enumerate(range(-xspin, xspin + 1, 2)):
        sz[im, im] = tm / 2.0
    ops = {"": ident, "P": sp, "M": sm, "Z": sz}
    site_basis.append(basis)
    site_ops.append(ops)

driver.initialize_system(n_sites=n_sites, vacuum=Q(0), target=Q(total_spin), hamil_init=False)
driver.ghamil = driver.get_custom_hamiltonian(site_basis, site_ops)

b = driver.expr_builder()
for (i, j, v) in topology:
    assert i < j
    b.add_term("PM", [i, j], 0.5 * v)
    b.add_term("MP", [i, j], 0.5 * v)
    b.add_term("ZZ", [i, j], 1.0 * v)
b.add_term("", [], ecore)
heis_mpo = driver.get_mpo(b.finalize(adjust_order=True, fermionic_ops=""), iprint=0)

bond_dims = [200] * 5 + [100] * 5 + [50] * 5
noises = [1E-4] * 4 + [1E-5] * 9 + [0]
dav_thrds = [1E-7] * 5 + [1E-9] * 5 + [1E-10]

all_eners = []
all_ssqs = []
all_zmats = []
all_scorrs = []

global_spins = list(range(total_spin, 4, 2))

b = driver.expr_builder()
for i, j in [(i, j) for i in range(n_sites) for j in range(n_sites)]:
    b.add_term("PM", [i, j], 0.5)
    b.add_term("MP", [i, j], 0.5)
    b.add_term("ZZ", [i, j], 1.0)
sq_mpo = driver.get_mpo(b.finalize(adjust_order=True, fermionic_ops=""), iprint=0)

print("\n === S = %d/2 Heisenberg %d-Site Model === \n" % (local_spin, n_sites))
for spin in global_spins:
    target = driver.target
    target.twos = target.twos_low = spin
    ket = driver.get_random_mps(tag='KET', bond_dim=250, nroots=8, target=Q(spin))
    eners = driver.dmrg(heis_mpo, ket, n_sweeps=20, iprint=2, bond_dims=bond_dims, noises=noises,
        thrds=dav_thrds, dav_max_iter=500, cutoff=0)
    all_eners.append(eners)
    print('2Sz = %d : E =' % spin, ('%9.4f' * len(eners)) % tuple(eners))
    sqs_expts = []
    zmats = []
    for ir in range(ket.nroots):
        iket = driver.split_mps(ket, ir, tag='GSTMP')
        sqs_expts.append(driver.expectation(iket, sq_mpo, iket))
        zmat = np.zeros((n_sites, ))
        for k in range(n_sites):
            xmpo = driver.get_mpo(driver.expr_builder().add_term("Z", [k], 1.0).finalize(adjust_order=False, fermionic_ops=""), iprint=0)
            zmat[k] = driver.expectation(iket, xmpo, iket)
        print('zmat', ir, zmat, np.sum(zmat), np.sum(np.abs(zmat)))
        zmats.append(zmat)
    print('<S^2> =', sqs_expts)
    all_ssqs.append(sqs_expts)
    all_zmats.append(zmats)
    scorr_expts = []
    for ir in range(ket.nroots):
        iket = driver.split_mps(ket, ir, tag='GSTMP')
        scorr = np.zeros((n_sites, n_sites))
        for i, j in [(i, j) for i in range(n_sites) for j in range(n_sites) if i <= j]:
            b = driver.expr_builder()
            b.add_term("PM", [i, j], 0.5)
            b.add_term("MP", [i, j], 0.5)
            b.add_term("ZZ", [i, j], 1.0)
            sc_mpo = driver.get_mpo(b.finalize(adjust_order=True, fermionic_ops=""), iprint=0)
            scorr[i, j] = scorr[j, i] = driver.expectation(iket, sc_mpo, iket)
        scorr_expts.append(scorr)
    all_scorrs.append(scorr_expts)

ener_min = np.min([np.min(x) for x in all_eners])
print('Emin = %20.10f' % ener_min)
print("\nRelative energies:\n")
for spin, eners in zip(global_spins, all_eners):
    eners = np.array(eners) - ener_min
    print('2Sz = %d : E =' % spin, ('%12.6f' * len(eners)) % tuple(eners))

print("\n<S^2> for each state:\n")
for spin, eners, ssqs in zip(global_spins, all_eners, all_ssqs):
    eners = np.array(eners) - ener_min
    print('2Sz = %d : SSQ =' % spin, ('%9.4f' * len(ssqs)) % tuple(ssqs))

print("\nSz for each state:\n")
for spin, zmats in zip(global_spins, all_zmats):
    for ix, zmat in enumerate(zmats):
        print('2Sz = %d : iroot = %d' % (spin, ix), ('%9.4f' * len(zmat)) % tuple(zmat))

print("\nSpin correlation for each state:\n")
for spin, scorr_expts in zip(global_spins, all_scorrs):
    for ix, scorr in enumerate(scorr_expts):
        print('2Sz = %d : iroot = %d' % (spin, ix))
        for i in range(n_sites):
            print('    %d ::' % i, ('%8.3f' * len(scorr[i])) % tuple(scorr[i]))

np.save('24-energies.npy', all_eners)
np.save('24-ssqs.npy', all_ssqs)
np.save('24-spin-corr.npy', all_scorrs)
