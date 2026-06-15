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
data |= json.load(open("../20-Data/06-llduc-fno-ucc.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/11-llduc-fno-udmrg.json", "r", encoding='utf-8'))

bs_map = data['name_x']

for ii, ix in enumerate([25, 26, 27]):
    for ir in range(1, 8):
        print('%10s & %10s & %8s & %8s & %8s & %8s & %8s \\\\' % (
            bs_map[str(ix)] if ir == 1 else '',
            '(%2do, %2de)' % (26 + ir * 4, 13 + ir * 8),
            '%8.2f' % data['fno_ucc'][str(ix)][0][str(ir)]['uccsd_ssq'][0],
            '%8.2f' % data['fno_ucc'][str(ix)][0][str(ir)]['uccsdt_ssq'][0],
            '%8.2f' % data['fno_ucc'][str(ix)][0][str(ir)]['uccsdtq_ssq'][0] if ir == 1 else '',
            '%8.2f' % data['fno_udmrg'][str(ix)][0][str(ir)]['udmrg_5000_ssq'][0] if ir <= 6 else '',
            '%8.2f' % data['fno_udmrg'][str(ix)][0][str(ir)]['udmrg_8000_ssq'][0] if ir <= 6 else '',
        ))
    if ii != 2:
        print("     \\hline")
