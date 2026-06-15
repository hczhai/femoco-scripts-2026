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

data = {'fe4s4-uhf': {ix: [] for ix in range(3)}}

for ix in range(3):
    with open("../04-Fe4S4/01-UHF-UCC-UDMRG/00-fe4-init-%d.out" % (ix, ), "r") as f:
        uhf_conv, uhf_ssq = True, None
        for l in f.readlines():
            if '::' in l and ' -1 ' in l:
                ll = l.strip().split()
                data['fe4s4-uhf'][ix].append({
                    'init': ll[2],
                    'uhf': float(ll[3]),
                    'uhf_conv': uhf_conv,
                    'uhf_ssq': uhf_ssq,
                })
                uhf_conv = True
            elif l.startswith('SCF not converged'):
                uhf_conv = False
            elif '<S^2> =' in l:
                uhf_ssq = [float(l.split()[-4]), float(l.split()[-1])]

with open("19-fe4s4-uhf.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Fe4S4: Number of unique UHF solutions: ', sum([len(data['fe4s4-uhf'][ix]) for ix in range(3)]))
