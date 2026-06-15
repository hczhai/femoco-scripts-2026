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
    for l in open('04-fm-ccsd-%d.out' % iiid, 'r').readlines():
        if '::' in l and int(l.strip().split()[-1]) <= 1:
            ll = l.strip().split()
            x[iiid].append([float(ll[3]), float(ll[4]), l.strip()])
    # x[iiid].sort()

ccsd_new = []
ccsd_t_new = []

for iiid in range(35):
    for f, f2, l in x[iiid]:
        ccsd_new.append(f - -22140)
        ccsd_t_new.append(f2 - -22140)
        print("%5d" % iiid, l)

# print('e_ccsd_new =', ccsd_new)
# print('e_ccsd_t_new =', ccsd_t_new)
