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

plt.rcParams.update({'font.size': 16})
plt.figure(figsize=(8.5, 6))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

bs_map = data['name_x']
colors = data['colors']
grays = data['grays']

def exp_log_fit(x, y):
    from scipy.optimize import curve_fit
    M = np.array(x)
    E = np.array(y)

    def E_model(M, E_inf, A, kappa):
        return E_inf + A * np.exp(-kappa * M)

    E_inf_init = -0.41
    A_init = 0.25
    kappa_init = 0.08
    p0 = [E_inf_init, A_init, kappa_init]

    popt, pcov = curve_fit(E_model, M, E, p0=p0, maxfev=10000)
    E_inf, A, kappa = popt
    perr = np.sqrt(np.diag(pcov))

    print(f"E_inf  = {E_inf:.8f} ± {perr[0]:.2e}")
    print(f"A      = {A:.8f} ± {perr[1]:.2e}")
    print(f"kappa = {kappa:.8f} ± {perr[2]:.2e}")

    poly_model = lambda xx, E_inf=E_inf, A=A, kappa=kappa: E_model(xx, E_inf, A, kappa)
    return poly_model, None, E_inf

e_udmrg_18k_25 = 22140 + min([x["energy"] for x in data["udmrg_18k"]["25"]["forward"][0].values() if x["bond"] == 18000])
data_25 = {x['bond']: (x['energy'], x['dw']) for _, x in sorted(data['udmrg_18k']["25"]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
dmrg_log_bond_25 = [np.log(x) ** 2 for x, _ in sorted(data_25.items())]
dmrg_inv_bond_25 = [1 / x for x, _ in sorted(data_25.items())]
dmrg_energy_25 = [22140 + y for _, (y, _) in sorted(data_25.items())]
print(bs_map['25'] + ":")
fits_25 = exp_log_fit(dmrg_log_bond_25, dmrg_energy_25)

e_udmrg_18k_26 = 22140 + min([x["energy"] for x in data["udmrg_18k"]["26"]["forward"][0].values() if x["bond"] == 18000])
data_26 = {x['bond']: (x['energy'], x['dw']) for _, x in sorted(data['udmrg_18k']["26"]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
dmrg_log_bond_26 = [np.log(x) ** 2 for x, _ in sorted(data_26.items())]
dmrg_inv_bond_26 = [1 / x for x, _ in sorted(data_26.items())]
dmrg_energy_26 = [22140 + y for _, (y, _) in sorted(data_26.items())]
print(bs_map['26'] + ":")
fits_26 = exp_log_fit(dmrg_log_bond_26, dmrg_energy_26)

plt.plot(1 / np.array([18000]), [e_udmrg_18k_25], linewidth=0, marker='X', markersize=7,
    mec=colors[4], color=colors[4], mfc='white', mew=1.5, label='UDMRG (forward, %s)' % bs_map['25'], zorder=5)
plt.plot(dmrg_inv_bond_25, dmrg_energy_25, linewidth=0, marker='p', markersize=7,
    mec=colors[4], color=colors[4], mfc='white', mew=1.5, label='UDMRG (reverse, %s)' % bs_map['25'], zorder=5)

plt.plot(1 / np.array([18000]), [e_udmrg_18k_26], linewidth=0, marker='X', markersize=7,
    mec=colors[6], color=colors[6], mfc='white', mew=1.5, label='UDMRG (forward, %s)' % bs_map['26'], zorder=5)
plt.plot(dmrg_inv_bond_26, dmrg_energy_26, linewidth=0, marker='p', markersize=7,
    mec=colors[6], color=colors[6], mfc='white', mew=1.5, label='UDMRG (reverse, %s)' % bs_map['26'], zorder=5)

x_reg = np.linspace(1E-14, 1/4500, 100)

plt.plot(x_reg, fits_25[0]((np.log(1 / x_reg) ** 2).reshape(-1, 1)).flatten(), '--',  linewidth=1.5,  color=grays[0],
    label='$\\exp [-\\kappa (\\ln D)^2]$ fit ($E_{\\mathrm{%s}}=%.6f$)' %  (bs_map['25'], fits_25[2], ), zorder=4)

plt.plot(x_reg, fits_26[0]((np.log(1 / x_reg) ** 2).reshape(-1, 1)).flatten(), '-',  linewidth=2.0,  color=grays[0],
    label='$\\exp [-\\kappa (\\ln D)^2]$ fit ($E_{\\mathrm{%s}}=%.6f$)' %  (bs_map['26'], fits_26[2], ), zorder=4)

plt.legend(loc='lower right', frameon=False, fontsize=12, ncols=1)
plt.xlim((0, 1 / 4500))
plt.ylabel('Total energy + 22140 (Hartree)')
plt.xlabel('$1/D$')
plt.xticks([0, 1/18000, 1/12000, 1/9000, 1/7000, 1/6000, 1/5000], ['0', '1/18k', '1/12k', '1/9k', '1/7k', '1/6k', '1/5k'])
plt.subplots_adjust(left=0.16, right=0.96, top=0.97)
plt.savefig('11-sm-fm-dmrg-energy.png', dpi=300)
