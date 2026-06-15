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

data = json.load(open("../20-Data/00-llduc-uhf.json", "r", encoding='utf-8'))
data |= json.load(open("../20-Data/80-bs-states.json", "r", encoding='utf-8'))

names = [data['name_x'][str(ix)] for ix in range(35)]
cnts = [len(data['uhf'][str(ix)]) for ix in range(35)]

pairs = sorted(zip(names, cnts), key=lambda x: (int(x[0].split('-')[0][2:]), x[0].split('-')[1:]))

nc = (len(pairs) + 4 - 1) // 4
for i in range(0, nc):
    print('%10s & %5d & %10s & %5d & %10s & %5d & %10s & %5s \\\\' % (
        *pairs[i], *pairs[i + nc], *pairs[i + nc * 2], *(pairs[i + nc * 3]
        if i + nc * 3 < len(pairs) else ['\\textbf{total}', '\\textbf{%d}' % sum(cnts)])))
