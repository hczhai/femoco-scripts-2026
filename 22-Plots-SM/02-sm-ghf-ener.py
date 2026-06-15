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

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(9.5, 6))
plt.grid(which='major', axis='y', linestyle='--')

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/07-llduc-ghf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/08-llduc-gccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

bs_map = data['name_x']
colors = data['colors']
inits = [25, 26, 27]

e_ref = min([v['energy'] for _, v in data['composite'].items() if 'energy' in v])

e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_t = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])

gid_map = {x['init']: x['id'] for x in data['ghf']['unique']}
e_ccsd = sorted([x['gccsd'], x['gccsd(t)'], x['init']] for x in data['gccsd'])[:16]
ghf_ids = [gid_map.get(x[2], '?') for x in e_ccsd]

e_gccsd = np.array([x[0] for x in e_ccsd])
e_gccsd_t = np.array([x[1] for x in e_ccsd])

print('Min GHF = %.6f' % min([x['ghf'] for x in data['ghf']['unique']]))
print('Min UHF = %.6f' % min([x['uhf'] for ix in range(35) for x in data['uhf'][str(ix)]]))
print('Diff    = %.2f kcal/mol' % ((
    min([x['uhf'] for ix in range(35) for x in data['uhf'][str(ix)]]) - 
    min([x['ghf'] for x in data['ghf']['unique']])) * ha_to_kcal_mol))

ik = 0
tks = []
for k in inits:
    plt.plot(np.arange(ik, ik + 1), [(e_uccsd[k] - e_ref) * ha_to_kcal_mol], '-', linewidth=2, marker='^', mec=colors[1], markersize=6, color=colors[1], mfc=colors[1], mew=2, label='UCCSD' if k == 25 else None, clip_on=False)
    plt.plot(np.arange(ik, ik + 1), [(e_uccsd_t[k] - e_ref) * ha_to_kcal_mol], '-', linewidth=2, marker='s', mec=colors[2], markersize=5, color=colors[2], mfc=colors[2], mew=2, label='UCCSD(T)' if k == 25 else None, clip_on=False)
    tks.append(ik)
    ik += 2.4

plt.axvline(x=ik - 0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)

xr = list(np.arange(ik, ik + 9)) + list(np.arange(ik + 9, ik + 9 + 7 * 1.5, 1.5))

plt.plot(xr, (e_gccsd - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='^', mec=colors[1], markersize=6, color=colors[1], mfc='white', mew=2, label='GCCSD', clip_on=False)
plt.plot(xr, (e_gccsd_t - e_ref) * ha_to_kcal_mol, '-', linewidth=2, marker='s', mec=colors[2], markersize=5, color=colors[2], mfc='white', mew=2, label='GCCSD(T)', clip_on=False)

tks += xr

plt.ylabel('Estimated total energy error (kcal/mol)')
plt.xticks(tks, [bs_map[str(ix)] for ix in inits] + ['%d' % x for x in ghf_ids])
plt.legend(loc='upper right', frameon=False, fontsize=16, ncols=2, bbox_to_anchor=(0.95, 0.7))
plt.subplots_adjust(left=0.08, right=0.97, top=0.95, bottom=0.075)
plt.savefig('02-sm-ghf-ener.png', dpi=300)
