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

hf_xstrs = open('01-fm-ghf-all.out', 'r').readlines()
e_hfs = {l.split()[3]: (il, float(l.split()[4])) for il, l in enumerate(hf_xstrs)}

for ix in range(72):
    for l in open('03-fm-gccsd-%d.out' % ix, 'r').readlines():
        if '::' in l and int(l.strip().split()[-1]) == 1:
            ll = l.strip().split()
            x.append([float(ll[4]), float(ll[5]), l.strip(), *e_hfs[ll[3]]])

for f, f2, l, idx, e_hf in sorted(x):
    print('%5d' % int(l.split()[0]), '::', l.split('::')[1], '%5d' % idx, e_hf)
