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


import numpy as np

conn = [
    [[0, 1], [0, 2], [0, 3]],
    [[1, 2], [1, 3], [2, 3]],
    [[4, 5], [4, 6], [5, 6]],
    [[4, 7], [5, 7], [6, 7]],
    [[1, 4], [2, 5], [3, 6]],
    [[2, 4], [3, 4], [1, 5], [3, 5], [1, 6], [2, 6]],
]

idx = np.array([1, 4, 3, 2, 5, 7, 6, 8])

names = ['left terminal', 'left cubane ', 'right cubane', 'right (Mo) terminal', 'central (axial)', 'central (cross)']
js = ['LT', 'LC', 'RC', 'RT', 'CA', 'CC']

for nm, j, cn in zip(names, js, conn):
    print('%25s & $J_{\\mathrm{%s}}$ & %5s & %5s & %5s & %5s & %5s & %5s \\\\' % (
        nm, j, *sorted(["(%d,%d)" % tuple(sorted((idx[x], idx[y]))) for x, y in cn]),
        *(['', '', ''] if j != 'CC' else [])
    ))
