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

data = {'gccsd': []}

for ix in range(72):
    with open("../06-FEMOCO-LLDUC/10-GHF/03-fm-gccsd-%d.out" % (ix, ), "r") as f:
        ghf, ghf_conv = 0.0, True
        for l in f.readlines():
            if '::' in l and int(l.strip().split()[-1]) <= 1:
                ll = l.strip().split()
                data['gccsd'].append({
                    'init': ll[3],
                    'gccsd': float(ll[4]),
                    'gccsd(t)': float(ll[5]),
                    'gccsd_conv': int(ll[6]) == 1,
                    'ghf': ghf,
                    'ghf_conv': ghf_conv,
                })
                ghf, ghf_conv = 0.0, True
            elif l.startswith('SCF not converged'):
                ghf_conv = False
            elif 'SCF energy =' in l:
                ghf = float(l.split('SCF energy =')[1].split()[0])

with open("08-llduc-gccsd.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of GCCSD calculations: ', len(data['gccsd']))
