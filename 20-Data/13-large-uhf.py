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

data = {'large_uhf': {ia: {ix: [] for ix in range(35)} | {-1: []} for ia in range(4)}}

for ia in range(4):
    for ix in range(35):
        with open("../07-FEMOCO-Large/01-UHF/%d-uhf-spin-%d.out" % (ia + 30, ix), "r") as f:
            uhf_conv, uhf_ssq = True, None
            for l in f.readlines():
                if '::' in l and ' -1 ' in l:
                    ll = l.strip().split()
                    ixx = int(l.strip().split('BS')[1].split()[0])
                    data['large_uhf'][ia][ixx].append({
                        'init': ll[3],
                        'proj': l.strip().split('BS')[0].split()[-1].strip(),
                        'uhf': float(ll[4]),
                        'uhf_conv': uhf_conv,
                        'uhf_ssq': uhf_ssq,
                    })
                    uhf_conv = True
                elif l.startswith('SCF not converged'):
                    uhf_conv = False
                elif '<S^2> =' in l:
                    uhf_ssq = [float(l.split()[-4]), float(l.split()[-1])]
    print('Model %d: Number of unique UHF solutions: ' % ia, sum([len(data['large_uhf'][ia][ix]) for ix in range(-1, 35)]))

with open("13-large-uhf.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
