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


fid = 30

eners = [[] for _ in range(35)]

for ix in range(35):
    with open('%02d-uhf-spin-%d.out' % (fid, ix), 'r') as f:
        for l in f.readlines():
            if '::' in l and int(l.split()[-1]) != -1:
                jx = int(l.split()[2])
                ex = float(l.split()[4])
                iix = int(l.split()[-1])
                eners[iix].append((jx, ex, l.strip()))

xcnt = 0

for ix in range(35):
    exx = sorted(eners[ix], key=lambda x: x[1])
    exv = []
    for xx in exx:
        if len(exv) == 0 or abs(xx[1] - exv[-1][1]) > 1E-8:
            exv.append(xx)
    for xx in exv[:24]:
        print("%5d :: %s" % (ix, xx[2].split(' :: ')[1].split('BS')[0]))
    xcnt += len(exv)

print(xcnt)