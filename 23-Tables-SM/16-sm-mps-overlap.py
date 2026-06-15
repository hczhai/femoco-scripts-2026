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
import scipy.stats

data = json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/04-llduc-udmrg-18k.json", "r", encoding='utf-8'))

bs_map = data['name_x']
dets = [[x for x in data["udmrg_18k"][str(ix)]["dets"][1] if abs(x) > 0.05] for ix in [25, 26]]

for ii, ix in enumerate([25, 26]):
    for ir in range(len(dets[ii]) // 5):
        print('%10s & %10.6f & %10.6f & %10.6f & %10.6f & %10.6f \\\\' % (
            bs_map[str(ix)] if ir == 0 else '', *dets[ii][ir * 5:(ir + 1) * 5]))
    if ii != 1:
        print("     \\hline")
