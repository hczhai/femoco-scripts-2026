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


eners = [[] for _ in range(36)]

for ix in range(35):
    for l in open('02-fm-init-%d.out' % ix, 'r').readlines():
        if '::' in l and '-1 ' in l:
            ll = l.strip().split()
            eners[int(ll[0])].append([float(ll[4]), l.rstrip() + " DIA"])
    for l in open('03-fm-init-pj-%d.out' % ix, 'r').readlines():
        if '::' in l and '-1 ' in l:
            ll = l.strip().split()
            eners[int(ll[0])].append([float(ll[4]), l.rstrip() + " OPT"])

xcnt = 0

for ix in range(36):
    exx = sorted(eners[ix])
    exv = []
    for xx in exx:
        if len(exv) == 0 or abs(xx[0] - exv[-1][0]) > 1E-8:
            exv.append(xx)
    for xx in exv[:24]:
        print(xx[1])
    xcnt += len(exv)

# print(xcnt)
