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

fid = 3

eners = [[] for _ in range(35)]

for ix in [2, 3, 4, 6, 7, 8, 9, 10, 11, 16, 17, 18, 19, 20, 21, 25, 26, 27]:
    for px in range(0, 24):
        found = False
        if os.path.exists('%02d-uccsd-act-%d-%d.out' % (fid, ix, px)):
            with open('%02d-uccsd-act-%d-%d.out' % (fid, ix, px), 'r') as f:
                for l in f.readlines():
                    if '::' in l:
                        jx = int(l.split()[1])
                        ex = float(l.split()[-4])
                        eners[ix].append((jx, ex, l.strip() + " %s" % l.split()[-1]))
                        found = True
        if os.path.exists('%02d-xrusty-uccsd-act-%d-%d.out' % (fid, ix, px)):
            with open('%02d-xrusty-uccsd-act-%d-%d.out' % (fid, ix, px), 'r') as f:
                for l in f.readlines():
                    if '::' in l:
                        jx = int(l.split()[1])
                        ex = float(l.split()[-4])
                        eners[ix].append((jx, ex, l.strip() + " %s" % l.split()[-1]))
                        found = True
        if not found:
            e_hf = 0.0
            e_ccsd = None
            with open('%02d-uccsd-act-%d-%d.out' % (fid, ix, px), 'r') as f:
                for l in f.readlines():
                    if 'converged SCF energy =' in l:
                        e_hf = float(l.split('converged SCF energy =')[1].split()[0])
                    elif 'E_corr(UCCSD) =' in l:
                        ex = float(l.split('E_corr(UCCSD) =')[1].split()[0])
                        e_ccsd = ((-1, ex + e_hf, "%.12f " % (ex + e_hf) + l.strip()))
                        found = True
            if found:
                eners[ix].append(e_ccsd)
        if not found:
            print(ix, px)

for ix in range(35):
    exx = sorted(eners[ix], key=lambda x: x[1])[:1]
    exv = []
    for xx in exx:
        if len(exv) == 0 or abs(xx[1] - exv[-1][1]) > 1E-8:
            exv.append(xx)
    for xx in exv:
        print('%5d %s' % (ix, xx[2]))
