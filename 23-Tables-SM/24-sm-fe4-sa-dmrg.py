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

data = json.load(open("../20-Data/18-fe4s4-sa-dmrg.json", "r", encoding='utf-8'))

sci_notation = lambda x: "$%.2f\\times 10^{-%.0f}$" % tuple([float(x) for x in ("%.2E" % x).split('E-')])

for ix in range(3):

    fwd_dmrg = {x['bond']: x for _, x in
        sorted(data['fe4s4-sadmrg-20k']['%s' % ix]['forward'][0].items(), key=lambda x: int(x[0]))}
    if ix != 1:
        rev_dmrg = {x['bond']: x for _, x in
            sorted(data['fe4s4-sadmrg-20k']['%s' % ix]['reverse'][0].items(), key=lambda x: int(x[0]))}

        bonds = sorted(rev_dmrg)
        dmrg_rev_dw = [rev_dmrg[bond]["dw"] for bond in bonds]
        dmrg_rev_energy = [rev_dmrg[bond]["energy"] for bond in bonds]
        reg = scipy.stats.linregress(dmrg_rev_dw, dmrg_rev_energy)

    bond = 20000
    print('%16s & %16s & %8d & %15.6f & %20s \\\\' % (
        "\\multirow{%d}{*}{BS%d}" % (1 if ix == 1 else 7, ix + 1), "forward",
        bond, fwd_dmrg[bond]["energy"], sci_notation(fwd_dmrg[bond]["dw"]),
    ))
    if ix != 1:
        print("    \\cline{2-5}")
        print("    & \\multirow{5}{*}{reverse}", end='')
        for bond in range(18000, 13000, -1000):
            print('    %s %s & %8d & %15.6f & %20s \\\\' % (
                "&" if bond != 18000 else "", "", bond, rev_dmrg[bond]["energy"], sci_notation(rev_dmrg[bond]["dw"]),
            ))
        print("    \\cline{2-5}")
        print('    & %16s & %8s & %15.6f & %20s \\\\' % (
            "extrapolated", "$\\infty$", reg.intercept, '',
        ))
    if ix != 2:
        print("    \\hline")
