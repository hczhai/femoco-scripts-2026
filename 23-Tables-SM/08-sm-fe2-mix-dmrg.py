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

data = json.load(open("../20-Data/10-fe2s2.json", "r", encoding='utf-8'))

ref = data['fe2s2']['mix']['udmrg']

e_udmrg = {x['bond']: x['energy'] for _, x in sorted(ref['reverse'].items(), key=lambda x: int(x[0]))}

for d in [100, 200, 400, 600, 800, 1000, 1200, 1600, 2000, 2400]:
    print('%6s & %12.6f & %10.6f & %10.6f & %10.6f & %10.6f \\\\' % (
        d, e_udmrg[d], ref['reverse_ssq'][str(d)][0], ref['reverse_ssq'][str(d)][1] - 1,
        ref['reverse_sz'][str(d)][0], ref['reverse_sz'][str(d)][1]
    ))
