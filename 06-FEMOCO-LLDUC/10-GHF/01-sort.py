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


x = []

for iiid in range(35):
    for l in open('00-fm-ghf-%d.out' % iiid, 'r').readlines():
        if '::' in l and '-1' in l:
            ll = l.strip().split()
            x.append([float(ll[3]), l.strip(), iiid])

x.sort()

x_uniq = []
for xx in x:
    if len(x_uniq) == 0 or abs(xx[0] - x_uniq[-1][0]) > 1E-8:
        x_uniq.append(xx)

for f, l, iiid in x_uniq:
    print("%5d" % iiid, l)
