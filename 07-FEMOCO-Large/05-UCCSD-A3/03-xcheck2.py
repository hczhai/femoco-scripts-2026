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


import os

for j in range(18):
    for i in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
        if open('03-uccsd-act-%d-%d.out' % (i, j), 'r').readlines()[-1].startswith('::'):
            continue
        elif os.path.exists('03-xrusty-uccsd-act-%d-%d.out' % (i, j)):
            if open('03-xrusty-uccsd-act-%d-%d.out' % (i, j), 'r').readlines()[-1].startswith('::'):
                continue
        print('i = %d j = %d' % (i, j))
