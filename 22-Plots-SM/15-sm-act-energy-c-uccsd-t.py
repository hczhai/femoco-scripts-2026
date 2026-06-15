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
plt.figure(figsize=(7.5, 6))
plt.grid(which='major', axis='y', linestyle='--')

ha_to_kcal_mol = 627.5094740631
data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/13-large-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))

bs_init = data['bs']
bs_map = data['name_x']
colors = data['colors']
grays = data['grays']
names = [bs_map[str(ix)] for ix in range(35)]
inits = [ix for ix in range(35) if int(names[ix].split('-')[0].split('BS')[1]) >= 7]

e_uhf = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['uhf'][str(ix)]])[0] for ix in inits])
e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in inits])
e_uccsd_t = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[1] for ix in inits])

e_uhf_ls = np.array([[min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in inits] for ia in range(4)])
e_uccsd_ls = np.array([[min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[0] for ix in inits] for ia in range(4)])
e_uccsd_t_ls = np.array([[min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)]])[1] for ix in inits] for ia in range(4)])

ik = 0
tks = []
for k in bs_init[6:]:
    kk = np.array([inits.index(x) for x in k])
    plt.plot(np.arange(ik, ik + len(k)), (e_uccsd_t[kk] - np.min(e_uccsd_t)) * ha_to_kcal_mol, '-.', linewidth=1.25, marker='s', mec=colors[2], markersize=6, color=colors[2], mfc='white', mew=1.5, label='UCCSD(T) (LLDUC)' if k == [18, 20, 25] else None, clip_on=False, alpha=0.7)
    plt.plot(np.arange(ik, ik + len(k)), (e_uccsd_t_ls[2][kk] - np.min(e_uccsd_t_ls[2])) * ha_to_kcal_mol, ':', linewidth=2, marker='s', mec=colors[2], markersize=6, color=colors[2], mfc='white', mew=2, label='UCCSD(T) (285o)' if k == [18, 20, 25] else None, clip_on=False, alpha=1)
    plt.plot(np.arange(ik, ik + len(k)), (e_uccsd_t_ls[3][kk] - np.min(e_uccsd_t_ls[3])) * ha_to_kcal_mol, '-', linewidth=2.5, marker='s', mec=colors[2], markersize=6, color=colors[2], mfc=grays[3], mew=2, label='UCCSD(T) (404o)' if k == [18, 20, 25] else None, clip_on=False, alpha=1)
    tks.append((ik + ik + len(k) - 1) / 2)
    plt.axvline(x=ik + len(k) - 0.5, color='gray', linewidth=0.5, linestyle='--', alpha=0.5)
    ik += len(k)

plt.ylabel('Relative energy (kcal/mol)')
plt.ylim((0, 30))
plt.xticks(tks, ['BS%d' % x for x in range(7, 11)])
plt.xlim((-0.5, 17.5))
plt.legend(loc='upper right', frameon=False, fontsize=16, ncols=1)
plt.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.075)
plt.savefig('15-sm-act-energy-c-uccsd-t.png', dpi=300)
