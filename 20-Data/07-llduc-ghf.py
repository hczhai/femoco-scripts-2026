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

data = {'ghf': {ix: [] for ix in range(35)}}

for ix in range(35):
    with open("../06-FEMOCO-LLDUC/10-GHF/00-fm-ghf-%d.out" % (ix, ), "r") as f:
        ghf_conv = True
        for l in f.readlines():
            if '::' in l and ' -1 ' in l:
                ll = l.strip().split()
                data['ghf'][ix].append({
                    'init_id': int(ll[1]),
                    'init': ll[2],
                    'ghf': float(ll[3]),
                    'ghf_conv': ghf_conv,
                })
                ghf_conv = True
            elif l.startswith('SCF not converged'):
                ghf_conv = False

all_ghfs = sorted([x for ix in range(35) for x in data['ghf'][ix]], key=lambda x: (x['ghf'], x['init_id'], x['init']))
data['ghf']['unique'] = []
for x in all_ghfs:
    if len(data['ghf']['unique']) == 0 or abs(x['ghf'] - data['ghf']['unique'][-1]['ghf']) > 1E-8:
        x['id'] = len(data['ghf']['unique'])
        data['ghf']['unique'].append(x)

with open("07-llduc-ghf.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of unique GHF solutions: ', len(data['ghf']['unique']))
