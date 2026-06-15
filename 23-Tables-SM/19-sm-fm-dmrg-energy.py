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

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

bs_map = data['name_x']

sci_notation = lambda x: "$%.2f\\times 10^{-%.0f}$" % tuple([float(x) for x in ("%.2E" % x).split('E-')])

for ii, ix in enumerate([25, 26]):
    bond = 18000
    fwd_dmrg = {x['bond']: x for _, x in
        sorted(data['udmrg_18k'][str(ix)]['forward'][0].items(), key=lambda x: int(x[0]))}
    print('%10s & %16s & %8d & %15.6f & %20s \\\\' % (
        bs_map[str(ix)], "forward", bond,
        fwd_dmrg[bond]["energy"], sci_notation(fwd_dmrg[bond]["dw"]),
    ))
    for bond in range(16000, 4000, -1000):
        rev_dmrg = {x['bond']: x for _, x in
            sorted(data['udmrg_18k'][str(ix)]['reverse'][0].items(), key=lambda x: int(x[0]))}
        print('%10s & %16s & %8d & %15.6f & %20s \\\\' % (
            '', "reverse", bond,
            rev_dmrg[bond]["energy"], sci_notation(rev_dmrg[bond]["dw"]),
        ))
    print('%10s & %16s & %8s & %15.6f & %20s \\\\' % (
        '', "$\\exp[-\\kappa(\\log D)^2]$ fit", "$\\infty$",
        data['udmrg_18k'][str(ix)]['exp-log-fit'], '',
    ))
    if ii != 1:
        print("     \\hline")
