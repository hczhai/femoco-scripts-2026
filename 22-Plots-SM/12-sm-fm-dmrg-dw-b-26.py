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


import numpy as np, json
import matplotlib.pyplot as plt
import scipy.stats

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']

e_udmrg_18k = 22140 + min([x["energy"] for x in data["udmrg_18k"]["26"]["forward"][0].values() if x["bond"] == 18000])
data_rev = {x['bond']: (x['energy'], x['dw']) for _, x in sorted(data['udmrg_18k']["26"]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
dmrg_rev_bond = [x for x, _ in sorted(data_rev.items())[::-1]]
dmrg_rev_energy = [22140 + y for _, (y, _) in sorted(data_rev.items())[::-1]]
dmrg_rev_dw = [z for _, (_, z) in sorted(data_rev.items())[::-1]]

data_fwd = {x['bond']: (x['energy'], x['dw']) for _, x in sorted(data['udmrg_18k']["26"]['forward'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
dmrg_fwd_bond = [x for x, _ in sorted(data_fwd.items())]
dmrg_fwd_energy = [22140 + y for _, (y, _) in sorted(data_fwd.items())]
dmrg_fwd_dw = [z for _, (_, z) in sorted(data_fwd.items())]

print(dmrg_rev_energy)
print(dmrg_rev_dw)
print(dmrg_rev_bond)

print(dmrg_fwd_energy)
print(dmrg_fwd_dw)
print(dmrg_fwd_bond)

reg = scipy.stats.linregress(dmrg_rev_dw[:3], dmrg_rev_energy[:3])
emin, emax = min(dmrg_rev_energy), max(dmrg_rev_energy)
eex, rex = reg.intercept, abs(reg.intercept - emin) / 5

print(eex)

x_reg = np.array([0, 1.25E-4])

plt.plot(np.array(dmrg_fwd_dw)*10**5, dmrg_fwd_energy, '--', marker='X', markersize=7, mec=grays[2], color=grays[2], mfc='white', mew=1.5, label='UDMRG (forward schedule)')
plt.plot(np.array(dmrg_rev_dw)*10**5, dmrg_rev_energy, '-', marker='p', markersize=7, mec=colors[4], color=colors[4], mfc='white', mew=2, label='UDMRG (reverse schedule)')
plt.plot(np.array(x_reg)*10**5, reg.intercept + reg.slope * x_reg, '--', linewidth=1, color=grays[2], label='linear fit ($R^2 = %.2f$)' % (reg.rvalue ** 2, ), zorder=-5)
plt.plot(np.array(dmrg_rev_dw[:3])*10**5, dmrg_rev_energy[:3], ' ', marker='p', markersize=7, mec='none', color='none', mfc=grays[2], mew=2, alpha=0.34)

plt.axhline(y=reg.intercept, color=grays[1], linewidth=1, linestyle='-', alpha=0.5, zorder=-9)
plt.text(10.9, reg.intercept, "$E_{\\text{extrapolated}} = %.6f$" % (reg.intercept, ),
    ha='right', va='center', color='none', fontsize=14,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8)
plt.text(10.9, reg.intercept, "$E_{\\text{extrapolated}} = %.6f$" % (reg.intercept, ),
    ha='right', va='center', color=grays[1], fontsize=14)

plt.axhline(y=dmrg_fwd_energy[-1], color=grays[1], linewidth=1, linestyle='-', alpha=0.5, zorder=-9)
plt.text(10.9, dmrg_fwd_energy[-1], "$E_{\\text{variational}} = %.6f$" % (dmrg_fwd_energy[-1], ),
    ha='right', va='center', color='none', fontsize=14,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8)
plt.text(10.9, dmrg_fwd_energy[-1], "$E_{\\text{variational}} = %.6f$" % (dmrg_fwd_energy[-1], ),
    ha='right', va='center', color=grays[1], fontsize=14)

plt.legend(loc='upper left', frameon=False, fontsize=14, ncols=1, bbox_to_anchor=(0, 1))
plt.xlim([0, 11])
plt.ylim([-0.4095, -0.380])
plt.xlabel("Largest discarded weight $(\\times 10^{-5})$")
plt.ylabel('Total energy + 22140 (Hartree)')
for dw, ee, bd in zip(dmrg_rev_dw, dmrg_rev_energy, dmrg_rev_bond):
    if bd == 5000:
        plt.text(dw * 10**5 - 0.2, ee + 0.0002, "$D=%d$" % bd, va='center', ha='right', fontsize=10, color=colors[4])
    else:
        plt.text(dw * 10**5 + 0.2, ee + 0.0002, "$D=%d$" % bd, va='top', ha='left', fontsize=10, color=colors[4])

for dw, ee, bd in zip(dmrg_fwd_dw, dmrg_fwd_energy, dmrg_fwd_bond):
    if bd == 18000:
        plt.text(dw * 10**5 - 0.2, ee, "$D=%d$" % bd, va='bottom', ha='right', fontsize=10, color=colors[0])
    else:
        plt.text(dw * 10**5 - 0.2, ee, "$D=%d$" % bd, va='center', ha='right', fontsize=10, color=colors[0])

plt.subplots_adjust(left=0.16, bottom=0.12, right=0.98, top=0.95)
plt.savefig('12-sm-fm-dmrg-dw-b-26.png', dpi=300)
