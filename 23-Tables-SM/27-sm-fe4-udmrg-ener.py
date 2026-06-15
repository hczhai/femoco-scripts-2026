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


import json
import scipy.stats

data = json.load(open("../20-Data/22-fe4s4-udmrg.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))

sci_notation = lambda x: "$%.2f\\times 10^{-%.0f}$" % tuple([float(x) for x in ("%.2E" % x).split('E-')])

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['0']['reverse'][0].items(), key=lambda x: int(x[0]))}
bonds = sorted(rev_dmrg)
dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
dmrg_rev_energy = [rev_dmrg[bond]["energy"] for bond in bonds]
reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)
e_ref_0 = reg.intercept

rev_dmrg = {x['bond']: x for _, x in
    sorted(data['fe4s4-sadmrg-20k']['2']['reverse'][0].items(), key=lambda x: int(x[0]))}
bonds = sorted(rev_dmrg)
dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
dmrg_rev_energy = [rev_dmrg[bond]["energy"] for bond in bonds]
reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)
e_ref_2 = reg.intercept

for ix in range(3):
    print('BS%d' % (ix + 1))
    fwd_dmrg = {x['bond']: x for _, x in
        sorted(data['fe4s4-udmrg'][str(ix)]['forward'][0].items(), key=lambda x: int(x[0]))}
    bond = 12000
    print('& %16s & %8d & %15.6f & %20s & %+12.6f \\\\' % (
        "forward", bond, fwd_dmrg[bond]["energy"], sci_notation(fwd_dmrg[bond]["dw"]),
        fwd_dmrg[bond]["energy"] - (e_ref_0 if ix == 0 else e_ref_2)
    ))
    rev_dmrg = {x['bond']: x for _, x in
        sorted(data['fe4s4-udmrg'][str(ix)]['reverse'][0].items(), key=lambda x: int(x[0]))}
    for bond in sorted(rev_dmrg)[::-1]:
        print('& %16s & %8d & %15.6f & %20s & %12s \\\\' % (
            "reverse", bond, rev_dmrg[bond]["energy"], sci_notation(rev_dmrg[bond]["dw"]),
            ''
        ))
    print('& %16s & %8s & %15.6f & %20s & %+12.6f \\\\' % (
            "$\\exp[-\\kappa(\\log D)^2]$ fit", "$\\infty$", data['fe4s4-udmrg'][str(ix)]['exp-log-fit'], '',
            data['fe4s4-udmrg'][str(ix)]['exp-log-fit'] - (e_ref_0 if ix == 0 else e_ref_2)))
    if ix != 2:
        print("    \\hline")
