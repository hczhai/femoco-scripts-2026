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


x = [[] for _ in range(35)]

for iiid in range(35):
    for l in open('05-fm-ccsd-%d.out' % iiid, 'r').readlines():
        if l.startswith('--- '):
            xtt = l.split()[-2]
        elif '::' in l and int(l.strip().split()[-1]) <= 1:
            ll = l.strip().split()
            x[iiid].append([float(ll[4]), float(ll[5]), l.rstrip() + ' ' + xtt])
    x[iiid].sort()

for iiid in range(35):
    for f, f2, l in x[iiid][:1]:
        print(l)
