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
import scipy.stats

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/11-llduc-fno-udmrg.json", "r", encoding='utf-8'))

e_udmrg = {int(k): v['exp-log-fit'] for k, v in data['udmrg_18k'].items()}
fno_ccsds = {int(k): np.array([v[0][str(ir)]['uccsd'] for ir in range(1, 8)]) for k, v in data['fno_ucc'].items()}
fno_ccsdts = {int(k): np.array([v[0][str(ir)]['uccsdt'] for ir in range(1, 8)]) for k, v in data['fno_ucc'].items()}
fno_ccsdtqs = {int(k): np.array([v[0][str(ir)]['uccsdtq'] for ir in range(1, 8)]) for k, v in data['fno_ucc'].items()}
fno_noccs = {int(k): np.array(v[0]['occ']) for k, v in data['fno_ucc'].items()}

e_uccsdt = np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
e_uccsd = np.array([min([(x['uccsd'], x['uccsdt_conv'], x['uccsd']) for x in data['uccsdt'][str(ix)]])[2] for ix in range(35)])
e_uhf = np.array([min([(x['uhf'], x['uccsdt_conv'], x['uhf']) for x in data['uccsdt'][str(ix)]])[2] for ix in range(35)])

regs = {k: scipy.stats.linregress(fno_noccs[k][-3:], (e_uccsdt[k] + fno_ccsdtqs[k] - fno_ccsdts[k])[-3:]) for k in fno_noccs}
extraploated_ccsdtq = {k: regs[k].intercept + regs[k].slope * 2.00 for k in fno_noccs}

fe2_mix_corr = [data['fe2s2']['mix'][x] - data['fe2s2']['mix'][y]
    for x, y in zip(['uccsd', 'uccsdt', 'uccsdtq', 'sa_dmrg_12k'],
        ['uhf', 'uccsd', 'uccsdt', 'uccsdtq'])]

femoco_corr = {k: [
    float(e_uccsd[k] - e_uhf[k]),
    float(e_uccsdt[k] - e_uccsd[k]),
    float(extraploated_ccsdtq[k] - e_uccsdt[k])
] for k in fno_noccs}

femoco_fe2_ratio = {k: [x / y for x, y in zip(femoco_corr[k], fe2_mix_corr)] for k in fno_noccs}
pred_post_q = {k: float(fe2_mix_corr[-1] * np.average([x / y for x, y in zip(femoco_corr[k], fe2_mix_corr)]))
     for k in fno_noccs}

composite_ucc = {k: extraploated_ccsdtq[k] + pred_post_q[k] for k in fno_noccs}
conse_composite_final = {k: float(e_udmrg[k] + composite_ucc[k]) / 2 for k in e_udmrg}
conse_composite_error = {k: float(e_udmrg[k] - composite_ucc[k]) / 2 for k in e_udmrg}
conse_composite_min = {k: float(min(e_udmrg[k], composite_ucc[k])) for k in e_udmrg}
composite_final = {k: float(extraploated_ccsdtq[k] + conse_composite_min[k]) / 2 for k in e_udmrg}
composite_error = {k: float(extraploated_ccsdtq[k] - conse_composite_min[k]) / 2 for k in e_udmrg}

print('full post_q', pred_post_q)
print('dmrg extrap', e_udmrg)
print('composite_ucc', composite_ucc)
print('conse composite_final', conse_composite_final)
print('conse composite_min', conse_composite_min)
print('conse composite_error', {k: v * ha_to_kcal_mol for k, v in conse_composite_error.items()})
print('extraploated_ccsdtq', extraploated_ccsdtq)
print('composite_final', composite_final)
print('composite_error', {k: v * ha_to_kcal_mol for k, v in composite_error.items()})

print('%10s E[54o] = %25.12f' % ('25', (e_uccsdt[25] + fno_ccsdtqs[25] - fno_ccsdts[25])[-1]))
print('%10s E[54o] = %25.12f' % ('26', (e_uccsdt[26] + fno_ccsdtqs[26] - fno_ccsdts[26])[-1]))
print('%10s E[fit] = %25.12f' % ('25', regs[25].intercept + regs[25].slope * 2.00))
print('%10s E[fit] = %25.12f' % ('26', regs[26].intercept + regs[26].slope * 2.00))

data = {
    'composite': {k: {
        'ucc': composite_ucc[k],
    } for k in composite_ucc} | {k: {
        'conse-energy-min': conse_composite_min[k],
        'conse-energy': conse_composite_final[k],
        'conse-error': conse_composite_error[k],
        'energy': composite_final[k],
        'error': composite_error[k],
        'ucc': composite_ucc[k],
    } for k in composite_final}
}

with open("26-llduc-composite.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
