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

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

bs_map = data['name_x']

sci_notation = lambda x: "$%.2f\\times 10^{%.0f}$" % tuple([float(x) for x in ("%.2E" % x).split('E+')])
ftime = lambda x: ("%5.1f" if x >= 10 else "%5.2f") % x

rev_dmrg_26 = {x['bond']: x for _, x in
    sorted(data['udmrg_18k'][str(26)]['reverse'][0].items(), key=lambda x: int(x[0])) if x['bond'] >= 5000}
bonds_26 = [x[0] for x in sorted(rev_dmrg_26.items())]
times_26 = [x[1]['time'] / 3600 * x[1]["nnodes"] * 96 for x in sorted(rev_dmrg_26.items())]
flops_26 = [x[1]['flop'] for x in sorted(rev_dmrg_26.items())]

pred_bdim = 393000
reg = scipy.stats.linregress(np.log(bonds_26), np.log(times_26))
cpu_hours_pred = np.exp(reg.intercept + reg.slope * np.log(pred_bdim))

reg = scipy.stats.linregress(np.log(bonds_26), np.log(flops_26))
flop_pred = np.exp(reg.intercept + reg.slope * np.log(pred_bdim))

for ii, ix in enumerate([25, 26]):
    bond = 18000
    fwd_dmrg = {x['bond']: x for _, x in
        sorted(data['udmrg_18k'][str(ix)]['forward'][0].items(), key=lambda x: int(x[0]))}
    fwd_nsweeps = [x['bond'] for _, x in data['udmrg_18k'][str(ix)]['forward'][0].items()].count(bond)

    print('%10s & %10s & %7d & %4d & %4d & %20s & %5s hours & %9d \\\\' % (
        bs_map[str(ix)], "forward", bond, fwd_dmrg[bond]["nnodes"], fwd_nsweeps,
        sci_notation(fwd_dmrg[bond]["flop"]), ftime(fwd_dmrg[bond]["time"] / 3600),
        round(fwd_dmrg[bond]["time"] / 3600 * fwd_dmrg[bond]["nnodes"] * 96),
    ))

    total_cpu_hours = fwd_dmrg[bond]["time"] / 3600 * fwd_dmrg[bond]["nnodes"] * 96 * fwd_nsweeps
    for bond in range(16000, 4000, -1000):
        rev_dmrg = {x['bond']: x for _, x in
            sorted(data['udmrg_18k'][str(ix)]['reverse'][0].items(), key=lambda x: int(x[0]))}
        rev_nsweeps = [x['bond'] for _, x in data['udmrg_18k'][str(ix)]['reverse'][0].items()].count(bond)
        print('%10s & %10s & %7d & %4d & %4d & %20s & %5s hours & %9d \\\\' % (
            '', "reverse", bond, rev_dmrg[bond]["nnodes"], rev_nsweeps,
            sci_notation(rev_dmrg[bond]["flop"]), ftime(rev_dmrg[bond]["time"] / 3600),
            round(rev_dmrg[bond]["time"] / 3600 * rev_dmrg[bond]["nnodes"] * 96),
        ))
        total_cpu_hours += rev_dmrg[bond]["time"] / 3600 * rev_dmrg[bond]["nnodes"] * 96 * rev_nsweeps
    print('%10s & & & & & \\multicolumn{2}{r}{\\textbf{total CPU hours:}} & \\textbf{%7s} \\\\' % (
            '', round(total_cpu_hours)
        ))
    print("     \\hline")
print('     \\textbf{estimated} & & %s & & & %s & & %s \\\\' % (
    pred_bdim, sci_notation(flop_pred), sci_notation(cpu_hours_pred)
))
