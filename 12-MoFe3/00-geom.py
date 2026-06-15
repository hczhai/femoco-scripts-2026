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
from pyscf.data.nist import BOHR

print(BOHR)

coord = np.array([
    [ 0.0000,  0.0000,  2.3619],
    [-2.9316,  0.0000, -1.8670],
    [ 1.4658,  2.5388, -1.8670],
    [ 1.4658, -2.5388, -1.8670],
    [ 4.0008,  0.0000,  0.4302],
    [-2.0004,  3.4648,  0.4302],
    [-2.0004, -3.4648,  0.4302],
    [ 0.0000,  0.0000, -5.0061],
    [-3.3942,  0.0000,  5.8276],
    [ 1.6971,  2.9395,  5.8276],
    [ 1.6971, -2.9395,  5.8276],
    [-6.6850,  0.0000, -3.7913],
    [ 3.3425,  5.7894, -3.7913],
    [ 3.3425, -5.7894, -3.7913],
    [-5.6926,  0.0000,  4.7848],
    [ 2.8463,  4.9300,  4.7848],
    [ 2.8463, -4.9300,  4.7848],
    [-8.5351,  0.0000, -2.0745],
    [ 4.2676,  7.3917, -2.0745],
    [ 4.2676, -7.3917, -2.0745],
])

coord *= BOHR

elems = ['Mo'] + ['Fe'] * 3 + ['S'] * 10 + ['H'] * 6

with open('00-mofe3.xyz', 'w') as f:
    f.write('%d\n\n' % len(coord))
    for em, x in zip(elems, coord):
        f.write('%4s%12.6f%12.6f%12.6f\n' % (em, *x))

print(np.linalg.norm(coord[0] - coord[1]))