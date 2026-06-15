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

plt.rcParams.update({'font.size': 12})
plt.figure(figsize=(9.5, 4))

data = json.load(open("../20-Data/90-colors.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

colors = data['colors']
grays = data['grays']

rev_dmrg_26 = {x['bond']: x for _, x in
    sorted(data['udmrg_18k'][str(26)]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
bonds_26 = [x[0] for x in sorted(rev_dmrg_26.items())]
flops_dmrg_26 = [x[1]['flop'] for x in sorted(rev_dmrg_26.items())]

pred_bdim = 393000
reg_drmg = scipy.stats.linregress(np.log(bonds_26), np.log(flops_dmrg_26))
flop_pred = np.exp(reg_drmg.intercept + reg_drmg.slope * np.log(pred_bdim))

noccs_26 = [58 - data['fno_ucc']['26'][0][str(ir)]['nfrozen'] for ir in range(1, 8)]
flops_cc_26 = data['fno_ucc']['26'][0]['flops']
reg_cc = scipy.stats.linregress(np.log(noccs_26[-2:]), np.log(flops_cc_26[-2:]))
flop_58 = np.exp(reg_cc.intercept + reg_cc.slope * np.log(58))

time_rev = sum(flops_dmrg_26) * 4 / 1.102E18
time_pred = flop_pred * 2 / 1.102E18
time_54o = flops_cc_26[-1] * data['fno_ucc']['26'][0]['7']['uccsdtq_niter'] / 1.102E18
time_76o = flop_58 * 100 / 1.102E18

print('UDMRG reverse =', time_rev)
print('UDMRG D=%sk  =' % (pred_bdim // 1000), time_pred)
print('UCC 54o       =', time_54o)
print('UCC 76o       =', time_76o)

costs = [time_54o, time_rev, time_76o, time_pred]

plt.barh(0, costs[0],  color=colors[6], edgecolor=grays[1], linewidth=0, height=0.2, alpha=0.75)
plt.barh(0, costs[0],  color="none", edgecolor=grays[1], linewidth=1.5, height=0.2)
plt.barh(1, costs[1],  color=colors[4], edgecolor=grays[1], linewidth=0, height=0.2, alpha=0.75)
plt.barh(1, costs[1],  color="none", edgecolor=grays[1], linewidth=1.5, height=0.2)
plt.barh(2, costs[2],  color=colors[6], edgecolor="white", hatch='////', linewidth=0, height=0.2, alpha=0.65)
plt.barh(2, costs[2],  color="none", edgecolor=grays[2], linewidth=1.5, height=0.2)
plt.barh(3, costs[3],  color=colors[4], edgecolor="white", hatch='////', linewidth=0, height=0.2, alpha=0.65)
plt.barh(3, costs[3],  color="none", edgecolor=grays[2], linewidth=1.5, height=0.2)

plt.gca().invert_yaxis()

plt.text(costs[0]*1.1, 0, ' %.2f seconds' % time_54o, ha='left', va='center')
plt.text(costs[1]*1.1, 1, ' %.1f seconds' % time_rev, ha='left', va='center')
plt.text(costs[2]*1.1, 2, ' %.1f seconds' % time_76o, ha='left', va='center')
plt.text(costs[3]*1.1, 3, ' %.1f hours' % (time_pred / 3600), ha='left', va='center')

plt.xlim((1E-1, 6E6))
plt.xscale('log')
plt.yticks(np.arange(4), [
    'UCCSTQ (54o, 69e)\n(extrapolation + composite estimate ~1 kcal/mol error)',
    'UDMRG ($D=18\\mathrm{k}$) (76o, 113e)\n(extrapolation to ~1 kcal/mol error)',
    'UCCSTQ (76o, 113e)\n(composite estimate ~1 kcal/mol error)',
    'UDMRG ($D=%s\\mathrm{k}$) (76o, 113e)\n(~1 kcal/mol extrapolation distance)' % (pred_bdim // 1000)])
plt.xlabel('Estimated wall time (in seconds) using Frontier CPU + GPU FLOPs          ')
plt.subplots_adjust(left=0.51, right=0.95, top=0.95, bottom=0.15)
plt.savefig('04-fig-4.png', dpi=300)
