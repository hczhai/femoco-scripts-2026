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


mp = {}

for l in open('08-uccsd-min-24-ref.out', 'r').readlines():
    xd = l.split()[3]
    mp[xd] = l.strip()

for l in open('08-uccsd-min-24-add.out', 'r').readlines():
    xd = l.split()[3]
    assert xd not in mp
    mp[xd] = l.strip()

eners = [[] for _ in range(35)]

for l in open('../28-FM-ACT/33-uhf-min-24.out', 'r').readlines():
    iiix = int(l.split()[0])
    if iiix not in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
        continue
    ll = mp[l.split()[3]]
    eners[iiix].append((float(ll.split()[5]), ll.split(' :: ')[1], int(ll.split()[0])))

for ix in range(35):
    exx = sorted(eners[ix], key=lambda x: x[0])
    for xx in exx[:1]:
        print('%5d :: %s INIT %d' % (ix, xx[1], xx[2]))
