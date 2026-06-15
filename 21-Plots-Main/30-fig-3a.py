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
plt.figure(figsize=(9.5, 6))

ha_to_kcal_mol = 627.5094740631

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/02-llduc-uccsdt.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/26-llduc-composite.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']

def exp_log_fit(x, y):
    from scipy.optimize import curve_fit
    M2, E = np.array(x), np.array(y)
    E_model = lambda M2, E_inf, A, kappa: E_inf + A * np.exp(-kappa * M2)

    E_inf_init, A_init, kappa_init = -0.41, 0.28, 0.02
    (E_inf, A, kappa), pcov = curve_fit(E_model, M2, E, xtol=1E-14, gtol=1e-14, ftol=1e-14, 
        p0=[E_inf_init, A_init, kappa_init], maxfev=20000)
    perr = np.sqrt(np.diag(pcov))

    poly_model = lambda xx, E_inf=E_inf, A=A, kappa=kappa: E_model(xx, E_inf, A, kappa)
    return poly_model, None, E_inf

e_udmrg_18k = 22140 + min([x["energy"] for x in data["udmrg_18k"]["26"]["forward"][0].values() if x["bond"] == 18000])
data_x = {x['bond']: (x['energy'], x['dw']) for _, x in sorted(data['udmrg_18k']["26"]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
dmrg_inv_bond = [1 / x for x, _ in sorted(data_x.items())]
dmrg_energy = [22140 + y for _, (y, _) in sorted(data_x.items())]
dmrg_dw = [z for _, (_, z) in sorted(data_x.items())]
fit = exp_log_fit([np.log(x) ** 2 for x, _ in sorted(data_x.items())], dmrg_energy)

nocc_26 = np.array(data['fno_ucc']['26'][0]['occ'])
e_uccsdt = 22140 + np.array([min([(x['uccsdt'], x['uccsdt_conv']) for x in data['uccsdt'][str(ix)]])[0] for ix in range(35)])
ccsdt_act_26 = np.array([data['fno_ucc']['26'][0][str(ir)]['uccsdt'] for ir in range(1, 8)])
ccsdtq_act_26 = np.array([data['fno_ucc']['26'][0][str(ir)]['uccsdtq'] for ir in range(1, 8)])

ccsdt_comp_26 = (e_uccsdt[26] + ccsdtq_act_26 - ccsdt_act_26)[-1]

reg_26 = scipy.stats.linregress(nocc_26[-3:], (e_uccsdt[26] + ccsdtq_act_26 - ccsdt_act_26)[-3:])
intercept_cc_26 = reg_26.intercept + reg_26.slope * 2.00

reg_dw = scipy.stats.linregress(dmrg_dw[-3:], dmrg_energy[-3:])
dmrg_dw_ex = reg_dw.intercept

conse_composite = 22140 + data['composite']['26']['conse-energy']
conse_error = data['composite']['26']['conse-error']

final_composite = 22140 + data['composite']['26']['energy']
final_error = data['composite']['26']['error']

pred_bdim = 393000
print('E[DMRG-18k-variational] = %25.16f' % e_udmrg_18k)
print('E[composite-UCCSDT-54o] = %25.16f' % ccsdt_comp_26)
print('E[DMRG-explog-fitting ] = %25.16f' % fit[2])
print('E[conse composite     ] = %25.16f' % conse_composite)
print('E[final composite     ] = %25.16f' % final_composite)
print('E[composite-UCCSDT-ex ] = %25.16f' % intercept_cc_26)
print('E[DMRG-dw-extrap      ] = %25.16f' % dmrg_dw_ex)
print('E[DMRG-%sk           ] = %25.16f' % (pred_bdim // 1000, fit[0](np.array([[np.log(pred_bdim)**2]])).flatten()[0]))
print('E[DMRG-%sk-error ha  ] = %25.16f' % (pred_bdim // 1000, fit[2] - fit[0](np.array([[np.log(pred_bdim)**2]])).flatten()[0]))
print('E[DMRG-%sk-error     ] = %25.16f' % (pred_bdim // 1000, (fit[2] - fit[0](np.array([[np.log(pred_bdim)**2]])).flatten()[0]) * ha_to_kcal_mol))
print('E[DMRG-18k-var above  ] = %25.16f' % ((e_udmrg_18k - fit[2]) * ha_to_kcal_mol))

x_reg = np.linspace(1E-14, 2.1E-4, 1000)

plt.axhline(y=e_udmrg_18k, color=grays[1], linewidth=1.4, linestyle='-', alpha=0.75, zorder=-9)
plt.text(2.1E-4, e_udmrg_18k, 'rigorous variational upper bound: UDMRG ($D$=18k)  ',
    ha='right', va='center', color=grays[1], fontsize=12,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8)
plt.text(2.12E-4, e_udmrg_18k, '$%.6f$' % e_udmrg_18k,
    ha='left', va='center', color=grays[1], fontsize=12)

plt.axhline(y=ccsdt_comp_26, color=colors[9], linewidth=1.4, linestyle='-', alpha=0.75, zorder=-9)
plt.text(2.1E-4, ccsdt_comp_26, 'high confidence upper bound: composite UCCSDT + $\\Delta$UCCSDTQ(54o)  ',
    ha='right', va='center', color=colors[9], fontsize=12, zorder=-8)
plt.text(2.1E-4, ccsdt_comp_26, 'high confidence upper bound: composite UCCSDT + $\\Delta$UCCSDTQ(54o)  ',
    ha='right', va='center', color='none', fontsize=12,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8.7)
plt.text(2.12E-4, ccsdt_comp_26, '$%.6f$' % ccsdt_comp_26,
    ha='left', va='center', color=colors[9], fontsize=12)

plt.axhline(y=intercept_cc_26, color=colors[9], linewidth=1.4, linestyle='-', alpha=0.75, zorder=-9)
plt.text(2.1E-4, intercept_cc_26, 'medium confidence upper bound: UCCSDTQ from extrapolation  ',
    ha='right', va='center', color=colors[9], fontsize=12, zorder=-8)
plt.text(2.1E-4, intercept_cc_26, 'medium confidence upper bound: UCCSDTQ from extrapolation  ',
    ha='right', va='center', color='none', fontsize=12,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8.7)
plt.text(2.12E-4, intercept_cc_26, '$%.6f$' % intercept_cc_26,
    ha='left', va='center', color=colors[9], fontsize=12)

plt.axhline(y=fit[2], color='#000000', linewidth=1.4, linestyle='-', alpha=0.75, zorder=-8.6)
plt.text(2.08E-4, fit[2] + 0.0001, 'UDMRG $\\exp [-\\kappa (\\ln D)^2]$ extrapolation ',
    ha='right', va='center', color='#000000', fontsize=12, zorder=-7.5)
plt.text(2.08E-4, fit[2] + 0.0001, 'UDMRG $\\exp [-\\kappa (\\ln D)^2]$ extrapolation ',
    ha='right', va='center', color='none', fontsize=12,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8.5)
plt.text(2.12E-4, fit[2] + 0.0001, '$%.6f$' % fit[2],
    ha='left', va='center', color='#000000', fontsize=12)

plt.axhline(y=final_composite, color='#000000', linewidth=1.4, linestyle='-', alpha=0.95, zorder=-8.0)
plt.text(2.08E-4, final_composite + 0.0001, 'final estimated BS8-237 energy',
    ha='right', va='center', color='#000000', fontsize=12, zorder=-7.8, fontweight="bold")
plt.text(2.08E-4, final_composite + 0.0001, 'final estimated BS8-237 energy',
    ha='right', va='center', color='none', fontsize=12, fontweight="bold",
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-7.9)
plt.text(2.11E-4, final_composite + 0.0001, '$\\mathbf{%.6f}$' % final_composite,
    ha='left', va='center', color='#000000', fontsize=12, fontweight="bold")

# plt.axhline(y=conse_composite, color='#000000', linewidth=1.4, linestyle='-', alpha=0.75, zorder=-9)
# plt.text(2.08E-4, conse_composite + 0.0000, 'consensus DMRG and CC extrapolation',
#     ha='right', va='center', color='#000000', fontsize=12, zorder=-7.6)
# plt.text(2.08E-4, conse_composite + 0.0000, 'consensus DMRG and CC extrapolation',
#     ha='right', va='center', color='none', fontsize=12,
#     bbox=dict(facecolor="white", edgecolor="none"), zorder=-8.7)
# plt.text(2.12E-4, conse_composite + 0.0000, '$%.6f$' % conse_composite,
#     ha='left', va='center', color='#000000', fontsize=12)

ucc_composite = 22140 + data['composite']['26']['ucc']
plt.axhline(y=ucc_composite, color=colors[6], linewidth=1.4, linestyle='-', alpha=0.75, zorder=-8.6)
plt.text(2.08E-4, ucc_composite - 0.0001, 'composite UCC extrapolation',
    ha='right', va='center', color=colors[6], fontsize=12, zorder=-8)
plt.text(2.08E-4, ucc_composite - 0.0001, 'composite UCC extrapolation',
    ha='right', va='center', color='none', fontsize=12,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8.55)
plt.text(2.12E-4, ucc_composite - 0.0001, '$%.6f$' % ucc_composite,
    ha='left', va='center', color=colors[6], fontsize=12)

plt.plot(1 / np.array([18000]), [e_udmrg_18k], linewidth=0, marker='X', markersize=7, mec=grays[1], color=grays[1], mfc='white', mew=1.75, label='UDMRG (forward schedule)', zorder=5)
plt.plot(dmrg_inv_bond, dmrg_energy, linewidth=0, marker='p', markersize=7, mec=colors[4], color=colors[4], mfc='white', mew=2, label='UDMRG (reverse schedule)', zorder=5)

plt.plot(x_reg, fit[0]((np.log(1 / x_reg) ** 2).reshape(-1, 1)).flatten(), '-',  linewidth=2.0,  color=grays[0], label='$E=E_0+A\\ \\exp [-\\kappa (\\ln D)^2]$ fit', zorder=4)

plt.axhspan(final_composite - final_error, final_composite + final_error,
    color=grays[2], hatch='////', linewidth=0, alpha=0.15, zorder=-7.5,
    label='uncertainty: $\\pm %.2f$ kcal/mol' % (final_error * ha_to_kcal_mol))

plt.legend(loc='upper left', frameon=False, fontsize=14, ncols=1)
plt.xlim((0.0, 2.1E-4))
plt.ylabel('Total energy + 22140 (Hartree)')
plt.xlabel('$1/D$')
plt.xticks([0, 1/18000, 1/12000, 1/9000, 1/7000, 1/6000, 1/5000], ['0', '1/18k', '1/12k', '1/9k', '1/7k', '1/6k', '1/5k'])
plt.subplots_adjust(left=0.14, right=0.88, top=0.97)
plt.savefig('30-fig-3a.png', dpi=300)
