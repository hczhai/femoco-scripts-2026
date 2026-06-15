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


import json, itertools

bs_init = [
    [34], # BS1
    [15], # BS2
    [0, 1, 5], # BS3
    [24, 29, 31], # BS4
    [22, 23, 28, 30, 32, 33], # BS5
    [12, 13, 14], # BS6
    [18, 20, 25], # BS7
    [16, 17, 19, 21, 26, 27], # BS8
    [2, 7, 11], # BS9
    [3, 4, 6, 8, 9, 10] # BS10
]
bs_flat = [x for x in bs_init for x in x]
assert sorted(bs_flat) == list(range(35))

name_x = {}
for ix, xs in enumerate(bs_init):
    for jx, x in enumerate(xs):
        name_x[x] = 'BS%d-%s' % (ix + 1, 'ABCDEF'[jx]) if len(xs) != 1 else 'BS%d' % (ix + 1, )

name_ijk = {}
for ix, ixx in enumerate(itertools.combinations(range(7), 3)):
    ijk = sorted([[0, 1, 4, 3, 2, 5, 7, 6][x + 1] for x in ixx])
    name_ijk[ix] = name_x[ix].split('-')[0] + '-' + ''.join('%s' % x for x in ijk)

# from Table S1 in SI of https://onlinelibrary.wiley.com/doi/full/10.1002/qua.25627
name_ref = {}
for ix, xs in enumerate(bs_init):
    if ix + 1 == 8: # special rule for BS8
        xmap = ['236', '246', '237', '347', '245', '345']
    else:
        xmap = sorted([name_ijk[x].split('-')[1] for x in xs])
    for x in xs:
        name_ref[x] = 'BS%d-%d' % (ix + 1, xmap.index(name_ijk[x].split('-')[1]) + 1) if len(xs) != 1 else 'BS%d' % (ix + 1, )

data = {
    'bs': bs_init,
    'name_x': name_x,
    'name_ijk': name_ijk,
    'name_ref': name_ref,
}

with open("80-bs-states.json", "w", encoding='utf-8') as f:
    json.dump(data, f, ensure_ascii=False, indent=4)
