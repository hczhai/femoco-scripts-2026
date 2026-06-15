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

plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(5, 5))

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/05-llduc-udmrg-8k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

bs_map = data['name_ijk']
colors = data['lights']
grays = data['grays']
inits = [25, 26, 27, 4]

e_ref = min([v['energy'] for _, v in data['composite'].items() if 'energy' in v])

e_uccsd_t = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])
e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_dmrg_5k = np.array([data['udmrg_5k'][str(ix)][0]['udmrg'] for ix in range(35)])
e_dmrg_8k = np.array([data['udmrg_8k'][str(ix)][0]['udmrg'] if len(data['udmrg_8k'][str(ix)]) != 0 else np.nan for ix in range(35)])

f = lambda x: (x - e_ref) * ha_to_kcal_mol

plt.bar(np.arange(4) - 0.27, f(e_uccsd_t[np.array(inits)]),  color=colors[0], edgecolor=grays[1], linewidth=1.5, width=0.18, label='UCCSD(T)')
plt.bar(np.arange(4) - 0.09, f(e_uccsdt[np.array(inits)]),   color=colors[1], edgecolor=grays[1], linewidth=1.5, width=0.18, label='UCCSDT')
plt.bar(np.arange(4) + 0.09, f(e_dmrg_5k[np.array(inits)]), color=colors[2], edgecolor=grays[1], linewidth=1.5, width=0.18, label='UDMRG ($D$=5k)')
plt.bar(np.arange(4) + 0.27, f(e_dmrg_8k[np.array(inits)]), color=colors[3], edgecolor=grays[1], linewidth=1.5, width=0.18, label='UDMRG ($D$=8k)')
plt.ylim((0, 47))
plt.xlim((-0.54, 3.54))
plt.legend(loc='upper left', frameon=False, fontsize=13, ncols=2, columnspacing=1)
plt.xticks(np.arange(4), [bs_map[str(x)] for x in inits])
plt.ylabel('Estimated total energy error (kcal/mol)')
plt.subplots_adjust(left=0.15, right=0.97, top=0.95, bottom=0.075)
plt.savefig('21-fig-2b.png', dpi=300)
