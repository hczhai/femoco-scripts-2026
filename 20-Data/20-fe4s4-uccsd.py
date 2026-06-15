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

data = {'fe4s4-uccsd': {ix: [] for ix in range(3)}}

for ix in range(3):
    with open("../04-Fe4S4/01-UHF-UCC-UDMRG/02-fe4-ccsd-%d.out" % (ix, ), "r") as f:
        uhf, uhf_conv = 0.0, True
        for l in f.readlines():
            if '::' in l and int(l.strip().split()[-1]) <= 1:
                ll = l.strip().split()
                data['fe4s4-uccsd'][ix].append({
                    'init': ll[2],
                    'uccsd': float(ll[3]),
                    'uccsd(t)': float(ll[4]),
                    'uccsd_conv': int(ll[5]) == 1,
                    'uhf': uhf,
                    'uhf_conv': uhf_conv,
                })
                uhf, uhf_conv = 0.0, True
            elif l.startswith('SCF not converged'):
                uhf_conv = False
            elif 'SCF energy =' in l:
                uhf = float(l.split('SCF energy =')[1].split()[0])

with open("20-fe4s4-uccsd.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Fe4S4: Number of UCCSD calculations: ', sum([len(data['fe4s4-uccsd'][ix]) for ix in range(3)]))
