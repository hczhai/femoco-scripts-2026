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


import matplotlib.pyplot as plt
import numpy as np, json
import scipy.stats

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 6))
plt.grid(which='major', axis='both', linestyle='--')

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']

fwd_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['0']['forward'][0].items(), key=lambda x: int(x[0]))}
rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['0']['reverse'][0].items(), key=lambda x: int(x[0]))}

bonds = sorted(fwd_dmrg)
dmrg_fwd_dw = [fwd_dmrg[bond]["dw"] for bond in bonds]
dmrg_fwd_energy = [327 + fwd_dmrg[bond]["energy"] for bond in bonds]

bonds = sorted(rev_dmrg)
dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
dmrg_rev_energy = [327 + rev_dmrg[bond]["energy"] for bond in bonds]
reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)

x_reg = np.array([0, 3E-5])
emin, emax = min(dmrg_rev_energy), max(dmrg_rev_energy)

plt.plot(dmrg_fwd_dw, dmrg_fwd_energy, '--', marker='X', markersize=7, mec=grays[2], color=grays[2], mfc='white', mew=1.5, label='SA-DMRG (forward schedule)')
plt.plot(x_reg, reg.intercept + reg.slope * x_reg, '--', linewidth=2.0, color=grays[1])
plt.plot(dmrg_rev_dw, dmrg_rev_energy, ' ', marker='p', markersize=7, mec=colors[4], color=colors[4], mfc='white', mew=2, label='SA-DMRG (reverse schedule)')
plt.text(0.5E-5, emin - 0.0012, "$E(D=\\infty) = %.6f \\pm %.6f$ [$R^2 = %.6f$]" %
    (reg.intercept, abs(reg.intercept - emin) / 5, reg.rvalue ** 2), color=colors[0], fontsize=13)
plt.xlim([0, 3E-5])
plt.ylim([-0.2453, -0.2432])
plt.legend(loc='upper left', frameon=False, fontsize=14, ncols=1, bbox_to_anchor=(0, 1))
plt.xlabel("Largest discarded weight")
plt.ylabel('Total energy + 327 (Hartree)')
plt.tight_layout()
plt.savefig('21-sm-fe4-extra-a.png', dpi=300)
