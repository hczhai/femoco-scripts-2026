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

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(5, 3.5))

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/70-dmrg-arxiv-2506-13386.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

colors = data["blues"]
ener_lmo_5k = data['femoco-lmo']['5000']
ener_emo_5k = data['femoco-emo']['5000']
ener_luo_5k = min([v[0]['udmrg'] for v in data['udmrg_5k'].values()])
e_ref = min([v['energy'] for _, v in data['composite'].items() if 'energy' in v])

print('LMO E = %25.12f' % ener_lmo_5k)
print('EMO E = %25.12f' % ener_emo_5k)
print('LUO E = %25.12f' % ener_luo_5k)
print('Eref  = %25.12f' % e_ref)

plt.bar([1, 2, 3], (np.array([ener_lmo_5k, ener_emo_5k, ener_luo_5k]) - e_ref) * ha_to_kcal_mol,
    width=0.2, color=colors[1], edgecolor=colors[0], linewidth=1.5, label='DMRG ($D$=5k)')
plt.ylabel('Total energy error (kcal/mol)')
plt.xticks([1, 2, 3], ['LMO\nRef.16', 'EMO\nRef.16', 'LUO\nthis work'])
plt.legend(loc='upper right', frameon=False, fontsize=12)
plt.tight_layout()
plt.savefig('10-fig-1f.png', dpi=300)
