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
plt.figure(figsize=(9.5, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/22-fe4s4-udmrg.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['0']['reverse'][0].items(), key=lambda x: int(x[0]))}

bonds = sorted(rev_dmrg)
dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
dmrg_rev_energy = [rev_dmrg[bond]["energy"] for bond in bonds]
reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)
dmrg_sa_extra_0 = 327 + reg.intercept

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['2']['reverse'][0].items(), key=lambda x: int(x[0]))}

bonds = sorted(rev_dmrg)
dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
dmrg_rev_energy = [rev_dmrg[bond]["energy"] for bond in bonds]
reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)
dmrg_sa_extra_2 = 327 + reg.intercept

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-udmrg']['0']['reverse'][0].items(), key=lambda x: int(x[0]))}
eners_0 = 327 + np.array([x['energy'] for _, x in rev_dmrg.items()])

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-udmrg']['1']['reverse'][0].items(), key=lambda x: int(x[0]))}
eners_1 = 327 + np.array([x['energy'] for _, x in rev_dmrg.items()])

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-udmrg']['2']['reverse'][0].items(), key=lambda x: int(x[0]))}
eners_2 = 327 + np.array([x['energy'] for _, x in rev_dmrg.items()])

dmrg_12k_0 = 327 + [x['energy'] for _, x in sorted(data['fe4s4-udmrg']['0']['forward'][0].items(),
    key=lambda x: int(x[0])) if x['bond'] == 12000][-1]
dmrg_12k_1 = 327 + [x['energy'] for _, x in sorted(data['fe4s4-udmrg']['1']['forward'][0].items(),
    key=lambda x: int(x[0])) if x['bond'] == 12000][-1]
dmrg_12k_2 = 327 + [x['energy'] for _, x in sorted(data['fe4s4-udmrg']['2']['forward'][0].items(),
    key=lambda x: int(x[0])) if x['bond'] == 12000][-1]

bdims = np.array([x['bond'] for _, x in rev_dmrg.items()])

xt = 2

def exp_log_fit(x, y):
    from scipy.optimize import curve_fit
    M = np.array(x)
    E = np.array(y)

    def E_model(M, E_inf, A, kappa):
        return E_inf + A * np.exp(-kappa * M)

    E_inf_init = -0.25
    A_init = 0.5
    kappa_init = 0.01
    p0 = [E_inf_init, A_init, kappa_init]

    popt, pcov = curve_fit(E_model, M, E, p0=p0, maxfev=10000)
    E_inf, A, kappa = popt
    perr = np.sqrt(np.diag(pcov))

    print(f"E_inf  = {E_inf:.8f} ± {perr[0]:.2e}")
    print(f"A      = {A:.8f} ± {perr[1]:.2e}")
    print(f"kappa = {kappa:.8f} ± {perr[2]:.2e}")

    poly_model = lambda xx, E_inf=E_inf, A=A, kappa=kappa: E_model(xx, E_inf, A, kappa)
    return poly_model, None, E_inf

model_0, coefficients, intercept_0 = exp_log_fit(np.log(bdims) ** xt, eners_0)
model_1, coefficients, intercept_1 = exp_log_fit(np.log(bdims) ** xt, eners_1)
model_2, coefficients, intercept_2 = exp_log_fit(np.log(bdims) ** xt, eners_2)

x_reg_x = np.linspace(1E-14, 3.5E-4, 600)

plt.plot(1 / np.array([12000]), [dmrg_12k_0], linewidth=0, marker='X', markersize=7, mec=colors[0], color=colors[0], mfc='white', mew=1.4, label='UDMRG (forward, BS1)', zorder=5)
plt.plot(1 / np.array(bdims), eners_0, linewidth=0, marker='p', markersize=7, mec=colors[0], color=colors[0], mfc='white', mew=2, label='UDMRG (reverse, BS1)', zorder=5)
plt.plot(x_reg_x, model_0((np.log(1 / x_reg_x) ** xt).reshape(-1, 1)).flatten(), '-',  linewidth=2,  color=grays[2], label='$\\exp [-\\kappa (\\ln D)^2]$ fit ($E_{\\mathrm{BS1}}=%.6f$)' %  (intercept_0, ), zorder=4)

plt.plot(1 / np.array([12000]), [dmrg_12k_1], linewidth=0, marker='X', markersize=7, mec=colors[1], color=colors[1], mfc='white', mew=1.4, label='UDMRG (forward, BS2)', zorder=5)
plt.plot(1 / np.array(bdims), eners_1, linewidth=0, marker='p', markersize=7, mec=colors[1], color=colors[1], mfc='white', mew=2, label='UDMRG (reverse, BS2)', zorder=5)
plt.plot(x_reg_x, model_1((np.log(1 / x_reg_x) ** xt).reshape(-1, 1)).flatten(), '-',  linewidth=2,  color=grays[2], label='$\\exp [-\\kappa (\\ln D)^2]$ fit ($E_{\\mathrm{BS2}}=%.6f$)' %  (intercept_1, ), zorder=4)

plt.plot(1 / np.array([12000]), [dmrg_12k_2], linewidth=0, marker='X', markersize=7, mec=colors[4], color=colors[4], mfc='white', mew=1.4, label='UDMRG (forward, BS3)', zorder=5)
plt.plot(1 / np.array(bdims), eners_2, linewidth=0, marker='p', markersize=7, mec=colors[4], color=colors[4], mfc='white', mew=2, label='UDMRG (reverse, BS3)', zorder=5)
plt.plot(x_reg_x, model_2((np.log(1 / x_reg_x) ** xt).reshape(-1, 1)).flatten(), '-',  linewidth=2,  color=grays[2], label='$\\exp [-\\kappa (\\ln D)^2]$ fit ($E_{\\mathrm{BS3}}=%.6f$)' %  (intercept_2, ), zorder=4)

plt.axhline(y=dmrg_sa_extra_0, color=grays[1], linewidth=1, linestyle='-', alpha=0.5, zorder=-9)
plt.text(3.47E-4, dmrg_sa_extra_0, "$E_{\\text{SA-DMRG extrapolated [BS1]}} = %.6f$" % (dmrg_sa_extra_0, ),
    ha='right', va='center', color='none', fontsize=14,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8)
plt.text(3.47E-4, dmrg_sa_extra_0, "$E_{\\text{SA-DMRG extrapolated [BS1]}} = %.6f$" % (dmrg_sa_extra_0, ),
    ha='right', va='center', color=grays[1], fontsize=14)

plt.axhline(y=dmrg_sa_extra_2, color=grays[1], linewidth=1, linestyle='-', alpha=0.5, zorder=-9)
plt.text(3.47E-4, dmrg_sa_extra_2, "$E_{\\text{SA-DMRG extrapolated [BS3]}} = %.6f$" % (dmrg_sa_extra_2, ),
    ha='right', va='center', color='none', fontsize=14,
    bbox=dict(facecolor="white", edgecolor="none"), zorder=-8)
plt.text(3.47E-4, dmrg_sa_extra_2, "$E_{\\text{SA-DMRG extrapolated [BS3]}} = %.6f$" % (dmrg_sa_extra_2, ),
    ha='right', va='center', color=grays[1], fontsize=14)

plt.legend(loc='lower right', frameon=False, fontsize=10, ncols=2)
plt.xlim((0.0, 3.5E-4))
plt.ylim((-0.253, -0.237))
plt.ylabel('Total energy + 327 (Hartree)')
plt.xlabel('$1/D$')
plt.xticks([0, 1/12000, 1/8000, 1/6000, 1/5000, 1/4000, 1/3000], ['0', '  1/12k', '1/8k', '1/6k', '1/5k', '1/4k', '1/3k'])
plt.subplots_adjust(left=0.142, right=0.97, top=0.97)
plt.savefig('23-sm-fe4-udmrg-fit.png', dpi=300)
