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
plt.figure(figsize=(6.5, 6))
plt.grid(which='major', axis='y', linestyle='--')

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/01-llduc-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/03-llduc-udmrg-5k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/05-llduc-udmrg-8k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/13-large-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/14-large-uccsd.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/17-large-ump2.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/15-qmmm-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/16-qmmm-uccsd.json", "r", encoding='utf-8'))

bs_map = data['name_ijk']
bs_init = data['bs']
colors = data['colors']

e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_qmmm = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['qmmm_ccsd'][str(ix)]])[0] for ix in range(35)])

ia = 3
e_uhf_03 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in range(35)])
e_uccsd_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x] + [(0, 0)])[0] for ix in range(35)])

ump2_map = {(x['init'], x['proj']): x for ix in range(35) for x in data['large_ump2'][str(ix)]}
e_init_03 = [min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv'], x['init'], x['proj']) for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x] + [(0, 0)])[3:] for ix in range(35)]
e_uhf_03 = np.array([ump2_map[k]['uhf_small'] if k != () else 0.0 for k in e_init_03])
e_ump2_03 = np.array([ump2_map[k]['ump2_small'] if k != () else 0.0 for k in e_init_03])
e_uhf_04 = np.array([ump2_map[k]['uhf_large'] if k != () else 0.0 for k in e_init_03])
e_ump2_04 = np.array([ump2_map[k]['ump2_large'] if k != () else 0.0 for k in e_init_03])
e_uccsd_03_comp = e_uccsd_03 + e_ump2_04 - e_ump2_03

e_uccsd = (e_uccsd - np.min(e_uccsd[~np.isnan(e_uccsd)])) * ha_to_kcal_mol
e_uccsd_qmmm = (e_uccsd_qmmm - np.min(e_uccsd_qmmm[~np.isnan(e_uccsd_qmmm)])) * ha_to_kcal_mol
e_uccsd_03_comp = (e_uccsd_03_comp - np.min(e_uccsd_03_comp[(e_uccsd_03_comp != 0) & ~np.isnan(e_uccsd_03_comp)])) * ha_to_kcal_mol

pp, dx = 0, 0.3

def separate(x):
    y = list(x)
    for i in range(len(x)):
        for j in range(len(x)):
            for k in range(len(x)):
                if y[k] > 5:
                    continue
                orig_i = y[i]
                if y[i] < y[j] and y[j] < y[k] and abs(y[i] - y[k]) < 0.4:
                    y[i], y[j], y[k] = (y[i] + y[j] + y[k]) / 3 - 0.2, (y[i] + y[j] + y[k]) / 3, (y[i] + y[j] + y[k]) / 3 + 0.2
                    if y[k] > 4.84:
                        dd = y[k] - 4.84
                        y[i], y[j], y[k] = y[i] - dd, y[j] - dd, y[k] - dd
                if y[i] > orig_i:
                    y[i] = orig_i
    for i in range(len(x)):
        for j in range(len(x)):
            if y[i] < y[j] and abs(y[i] - y[j]) < 0.2:
                y[i], y[j] = (y[i] + y[j]) / 2 - 0.1, (y[i] + y[j]) / 2 + 0.1
    for k in range(len(x)):
        if y[k] > 4.84:
            dd = y[k] - 4.84
            y[k] = y[k] - dd
    return y

for vpx, (eners, cl) in enumerate(zip(
    [e_uccsd, e_uccsd_qmmm, e_uccsd_03_comp],
    [colors[1], colors[1], colors[1]])):
    xeners = separate(eners)
    for ik, k in enumerate(bs_init):
        for ikk, kk in enumerate(k):
            if eners[kk] > 5:
                continue
            ls = '-'
            plt.plot([pp - dx - 0.45, pp + dx - 0.45], [eners[kk], eners[kk]], ls, linewidth=2.25, color=cl, label='UHF (LLDUC)' if k == [34] else None, clip_on=False)
            plt.text(pp + dx - 0.12 - 0.25, xeners[kk], bs_map[str(kk)], va='center', fontsize=14)
    pp += 2

plt.xticks([x for x in range(0, 12, 2)],
    ['UCCSD\n(%s)' % x for x in ['LLDUC', 'QM/MM', '404o']] +
    ['UCCSD(T)\n(%s)' % x for x in ['LLDUC', 'QM/MM', '404o']]
)

plt.ylabel('Relative energy (kcal/mol)')
plt.ylim((-0.2, 5))
plt.xlim((-1, 5.2))
plt.subplots_adjust(left=0.11, right=0.985, top=0.95, bottom=0.1)
plt.savefig('50-fig-5b.png', dpi=300)
