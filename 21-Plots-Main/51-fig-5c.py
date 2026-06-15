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
plt.figure(figsize=(8.4, 6))
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
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

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

e_uccsd = (e_uccsd - np.min(e_uccsd[~np.isnan(e_uccsd)])) * ha_to_kcal_mol
e_uccsd_t = (e_uccsd_t - np.min(e_uccsd_t[~np.isnan(e_uccsd_t)])) * ha_to_kcal_mol
e_uccsdt = (e_uccsdt - np.min(e_uccsdt[~np.isnan(e_uccsdt)])) * ha_to_kcal_mol
e_udmrg_5k = (e_udmrg_5k - np.min(e_udmrg_5k[~np.isnan(e_udmrg_5k)])) * ha_to_kcal_mol
e_udmrg_18k = (e_udmrg_18k - np.min(e_udmrg_18k[~np.isnan(e_udmrg_18k)])) * ha_to_kcal_mol
e_dmrg_ex = (e_dmrg_ex - np.min(e_dmrg_ex[~np.isnan(e_dmrg_ex)])) * ha_to_kcal_mol
e_uccsdt_composite = (e_uccsdt_composite - np.min(e_uccsdt_composite[~np.isnan(e_uccsdt_composite)])) * ha_to_kcal_mol

ia = 3
e_uhf_03 = np.array([min([(x['uhf'], x['uhf_conv']) for x in data['large_uhf'][str(ia)][str(ix)]])[0] for ix in range(35)])
e_uccsd_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x] + [(0, 0)])[0] for ix in range(35)])
e_uccsd_t_03 = np.array([min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv']) for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x] + [(0, 0)])[1] for ix in range(35)])

ump2_map = {(x['init'], x['proj']): x for ix in range(35) for x in data['large_ump2'][str(ix)]}
e_init_03 = [min([(x['uccsd'], x['uccsd(t)'], x['uccsd_conv'], x['init'], x['proj']) for x in data['large_uccsd'][str(ia)][str(ix)] if 'uccsd' in x] + [(0, 0)])[3:] for ix in range(35)]
e_uhf_03 = np.array([ump2_map[k]['uhf_small'] if k != () else 0.0 for k in e_init_03])
e_ump2_03 = np.array([ump2_map[k]['ump2_small'] if k != () else 0.0 for k in e_init_03])
e_uhf_04 = np.array([ump2_map[k]['uhf_large'] if k != () else 0.0 for k in e_init_03])
e_ump2_04 = np.array([ump2_map[k]['ump2_large'] if k != () else 0.0 for k in e_init_03])

e_uccsd_03_comp = e_uccsd_03 + e_ump2_04 - e_ump2_03
e_uccsd_t_03_comp = e_uccsd_t_03 + e_ump2_04 - e_ump2_03

e_uccsd_03_comp = (e_uccsd_03_comp - np.min(e_uccsd_03_comp[(e_uccsd_03_comp != 0) & ~np.isnan(e_uccsd_03_comp)])) * ha_to_kcal_mol
e_uccsd_t_03_comp = (e_uccsd_t_03_comp - np.min(e_uccsd_t_03_comp[(e_uccsd_t_03_comp != 0) & ~np.isnan(e_uccsd_t_03_comp)])) * ha_to_kcal_mol

e_cc_comp = {int(k): data['composite'][k]['ucc'] * ha_to_kcal_mol for k in data['composite']}

g_ccsd_new = np.min(e_uccsd[e_uccsd != 0])
g_ccsd_t_new = np.min(e_uccsd_t[e_uccsd_t != 0])
g_ccsdt = np.min(e_uccsdt[e_uccsdt != 0])
g_ccsdtq = e_uccsdt_composite[26] - e_uccsdt_composite[25]
g_cc_comp = e_cc_comp[25] - e_cc_comp[26]
g_dmrg_5k = np.min(e_udmrg_5k[e_udmrg_5k != 0])
g_dmrg_18k = e_udmrg_18k[25] - e_udmrg_18k[26]
g_dmrg_ex = abs(e_dmrg_ex[25] - e_dmrg_ex[26])

g_uccsd_03 = np.min(e_uccsd_03_comp[(e_uccsd_03_comp != 0) & ~np.isnan(e_uccsd_03_comp)])
g_uccsd_t_03 = np.min(e_uccsd_t_03_comp[(e_uccsd_t_03_comp != 0) & ~np.isnan(e_uccsd_t_03_comp)])

plt.plot([0, 2, 4, 6, 8], [g_ccsd_new, g_ccsd_t_new, g_ccsdt, g_ccsdtq, g_cc_comp], '-', linewidth=2, marker='o', mec=colors[1], markersize=8, color=colors[1], mfc='white', mew=2, label='UCC (LLDUC)')
plt.plot([10, 12, 14], [g_dmrg_5k, g_dmrg_18k, g_dmrg_ex], '-', linewidth=2, marker='^', mec=colors[4], markersize=8, color=colors[4], mfc='white', mew=2, label='UDMRG (LLDUC)')
plt.plot([16, 18], [g_uccsd_03, g_uccsd_t_03], '--', linewidth=2, marker='s', mec=colors[2], markersize=8, color=colors[2], mfc='white', mew=2, label='UCC (404o, 277e) $+ \\Delta$UMP2 (994o, 423e)')
plt.axvline(x=9, color='gray', linewidth=1, linestyle='--', alpha=0.5)
plt.axvline(x=15, color='gray', linewidth=1, linestyle='--', alpha=0.5)

plt.xticks([x for x in range(0, 20, 2)],
    ['CCSD', 'CCSD(T)', 'CCSDT', 'CCSDTQ\n(FNO)', 'CC\n(comp)', '$D=5k$', '$D=18k$', '$D\\rightarrow \\infty$', 'CCSD', 'CCSD(T)']
)

plt.ylabel('Spin isomer gap (kcal/mol)')
plt.ylim((-0, 5))
plt.xlim((-0.5, 19.2))
plt.legend(loc='upper left', frameon=False, fontsize=15, ncols=1)
plt.subplots_adjust(left=0.07, right=0.985, top=0.95, bottom=0.1)
plt.savefig('51-fig-5c.png', dpi=300)
