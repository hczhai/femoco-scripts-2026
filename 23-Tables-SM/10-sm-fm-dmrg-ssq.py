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

data = json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

ix = 26
ds = sorted([int(x) for x in data['udmrg_18k'][str(ix)]['reverse_ssq'][0]])

e_dmrg = {x['bond']: x['energy'] for _, x in
    sorted(data['udmrg_18k'][str(ix)]['reverse'][0].items(), key=lambda x: int(x[0]))}

for k in range(len(ds) // 2):
    print('%7d & %12.6f & %8.2f & %7d & %12.6f & %8.2f \\\\' % (
        ds[k], 22140 + e_dmrg[ds[k]], data['udmrg_18k'][str(ix)]['reverse_ssq'][0][str(ds[k])][0],
        ds[k + 3], 22140 + e_dmrg[ds[k + 3]], data['udmrg_18k'][str(ix)]['reverse_ssq'][0][str(ds[k + 3])][0],
    ))
