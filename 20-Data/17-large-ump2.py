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

data = {'large_ump2': {ix: [] for ix in range(35)} | {-1: []}}

for ix in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
    bs_map = {(l.split()[3], l.split()[-3]): int(l.split('BS')[1].strip())
        for l in open("../07-FEMOCO-Large/01-UHF/33-uhf-spin-%d.out" % (ix, ), "r").readlines() if '::' in l}
    with open("../07-FEMOCO-Large/06-UMP2-A4/74-ump2-composite-%d.out" % (ix, ), "r") as f:
        uhf_conv = True
        for l in f.readlines():
            if '::' in l:
                ll = l.strip().split()
                iix = bs_map[(ll[2], ll[-1])]
                data['large_ump2'][iix].append({
                    'init': ll[2],
                    'proj': ll[-1].strip(),
                    'uhf_small': float(ll[3]),
                    'ump2_small': float(ll[4]),
                    'uhf_large': float(ll[5]),
                    'ump2_large': float(ll[6]),
                    'uhf_conv': uhf_conv,
                })
                uhf_conv = True
            elif l.startswith('SCF not converged'):
                uhf_conv = False


for ix in [6, 7, 8, 9, 10]:
    with open("../07-FEMOCO-Large/06-UMP2-A4/76-ump2-composite-%d.out" % (ix, ), "r") as f:
        uhf_conv = True
        for l in f.readlines():
            if '::' in l:
                ll = l.strip().split()
                assert (ll[2], ll[-1]) not in [(x['init'], x['proj']) for x in data['large_ump2'][ix]]
                data['large_ump2'][ix].append({
                    'init': ll[2],
                    'proj': ll[-1].strip(),
                    'uhf_small': float(ll[3]),
                    'ump2_small': float(ll[4]),
                    'uhf_large': float(ll[5]),
                    'ump2_large': float(ll[6]),
                    'uhf_conv': uhf_conv,
                })
                uhf_conv = True
            elif l.startswith('SCF not converged'):
                uhf_conv = False

with open("17-large-ump2.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UMP2 calculations: ', sum([len(data['large_ump2'][ix]) for ix in range(35)]))
