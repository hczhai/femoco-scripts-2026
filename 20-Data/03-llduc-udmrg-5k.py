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


import json, numpy as np

data = {'udmrg_5k': {ix: [] for ix in range(35)}}

for ix in range(35):
    if ix in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
        fid, re = 0, ''
    elif ix in [1, 5, 12, 13, 15, 23, 30, 33, 34]:
        fid, re = 20, '-re'
    elif ix in [32]:
        fid, re = 20, '-re2'
    else:
        fid, re = 20, ''
    d = {}
    with open("../06-FEMOCO-LLDUC/04-UDMRG-5k/%02d-fm-dmrg-uccsd%s-%d.out" % (fid, re, ix), "r") as f:
        for l in f.readlines():
            if l.startswith('xdstr ='):
                d['init'] = l.split()[2]
            elif l.startswith('STATE 0 DMRG energy ='):
                d['udmrg'] = float(l.split()[5])
        assert 'init' in d and 'udmrg' in d
        data['udmrg_5k'][ix].append(d)

for ix in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
    dm1 = np.load('../06-FEMOCO-LLDUC/04-UDMRG-5k/data/02-%02d-1pdm.npy' % ix)
    mo = np.load('../06-FEMOCO-LLDUC/04-UDMRG-5k/data/00-mo-coeff-%02d.npy' % ix)
    sfdm = mo[0] @ dm1[0] @ mo[0].T + mo[1] @ dm1[1] @ mo[1].T
    data['udmrg_5k'][ix][0]['nat-occ'] = np.linalg.eigh(sfdm)[0][::-1].tolist()

with open("03-llduc-udmrg-5k.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)

print('Number of UDMRG (D=5k) calculations: ', sum([len(data['udmrg_5k'][ix]) for ix in range(35)]))
