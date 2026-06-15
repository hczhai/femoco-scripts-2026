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
plt.figure(figsize=(12, 6))
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

bs_init = data['bs']
bs_map = data['name_x']
colors = data['colors']

e_uhf = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['uhf'][str(ix)]])[0] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[0] for ix in range(35)])
e_uccsd_t = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['uccsd'][str(ix)]])[1] for ix in range(35)])
e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_udmrg_5k = np.array([data['udmrg_5k'][str(ix)][0]['udmrg'] for ix in range(35)])
e_udmrg_8k = np.array([data['udmrg_8k'][str(ix)][0]['udmrg'] if len(data['udmrg_8k'][str(ix)]) != 0 else np.nan for ix in range(35)])

e_udmrg_18k = np.array([min([x["energy"] for x in data["udmrg_18k"][str(ix)]["forward"][0].values() if x["bond"] == 18000])
    if str(ix) in data['udmrg_18k'] else np.nan for ix in range(35)])
e_dmrg_ex = np.array([data['udmrg_18k'][str(ix)]['exp-log-fit'] if str(ix) in data['udmrg_18k'] else np.nan for ix in range(35)])

e_fno_uccsdt = np.array([data['fno_ucc'][str(ix)][0]['7']['uccsdt'] if str(ix) in data['fno_ucc'] else np.nan for ix in range(35)])
e_fno_uccsdtq = np.array([data['fno_ucc'][str(ix)][0]['7']['uccsdtq'] if str(ix) in data['fno_ucc'] else np.nan for ix in range(35)])
e_uccsdt_composite = e_uccsdt + e_fno_uccsdtq - e_fno_uccsdt

print('E[dmrg ex     BS8-E] =', e_dmrg_ex[26])
print('E[uccsdtq 54o BS8-E] =', e_uccsdt_composite[26])
print('E[uccsdtq dmrg diff] =', (e_uccsdt_composite[26] - e_dmrg_ex[26]))
print('E[uccsdtq dmrg diff] =', (e_uccsdt_composite[26] - e_dmrg_ex[26]) * ha_to_kcal_mol)
print()
print('E[dmrg ex     BS7-C] =', e_dmrg_ex[25])
print('E[uccsdtq 54o BS7-C] =', e_uccsdt_composite[25])
print('E[uccsdtq dmrg diff] =', (e_uccsdt_composite[25] - e_dmrg_ex[25]))
print('E[uccsdtq dmrg diff] =', (e_uccsdt_composite[25] - e_dmrg_ex[25]) * ha_to_kcal_mol)

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
    return y

f = lambda x: (x - np.min(x[~np.isnan(x)])) * ha_to_kcal_mol

theories = [e_uhf, e_uccsd, e_uccsd_t, e_uccsdt, e_uccsdt_composite,
    e_udmrg_5k, e_udmrg_8k, e_udmrg_18k, e_dmrg_ex]

for vpx, (eners, cl) in enumerate(zip(theories, [*colors[:4], colors[6], *[colors[4]] * 4])):
    eners = f(eners)
    xeners = separate(eners)
    for ik, k in enumerate(bs_init):
        for ikk, kk in enumerate(k):
            if eners[kk] > 5:
                continue
            ls = '-'
            plt.plot([pp - dx - 0.45, pp + dx - 0.45], [eners[kk], eners[kk]], ls, linewidth=2.25, color=cl, label='UHF (LLDUC)' if k == [34] else None, clip_on=False)
            plt.text(pp + dx - 0.12 - 0.25, xeners[kk], bs_map[str(kk)], va='center', fontsize=14)
    pp += 2
    if vpx in [0, 4]:
        plt.axvline(x=pp - 0.87, color='gray', linewidth=1, linestyle='--', alpha=0.5)

plt.xticks([x for x in range(0, len(theories) * 2, 2)],
    ['UHF', 'UCCSD', 'UCCSD(T)', 'UCCSDT', 'UCCSDTQ\n(54o+comp)',
     'UDMRG\n($D=\\mathrm{5k}$)', 'UDMRG\n($D=\\mathrm{8k}$)',
     'UDMRG\n($D=\\mathrm{18k}$)', 'UDMRG\n($D=\\infty$)']
)

plt.ylabel('Relative energy (kcal/mol)')
plt.ylim((-0.2, 5))
plt.xlim((-1, len(theories) * 2 - 0.8))
plt.subplots_adjust(left=0.05, right=0.985, top=0.95, bottom=0.1)
plt.savefig('06-sm-llduc-bs-levels.png', dpi=300)
