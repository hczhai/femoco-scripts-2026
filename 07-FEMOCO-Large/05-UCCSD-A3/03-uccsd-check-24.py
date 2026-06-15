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

xdstrs = open('03-uhf-min-24.out', 'r').readlines()

xd_map = {(ix, jx): xstr for l in xdstrs for ll in [l.split()] for ix, jx, xstr in [[int(ll[0]), int(ll[2]), ll[3]]]}

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
                    if l.startswith('--- '):
                        iix, jjx, ttx = int(l.split()[1]), int(l.split()[3]), l.split()[4]
                    elif 'converged SCF energy =' in l:
                        e_hf = float(l.split('converged SCF energy =')[1].split()[0])
                    elif 'E_corr(UCCSD) =' in l:
                        ex = float(l.split('E_corr(UCCSD) =')[1].split()[0])
                        ll = ":: %5d %s %20.12f %20.12f %20.12f %d %s" % (jjx, xd_map[(iix, jjx)], e_hf, ex + e_hf, 0.0, -2, ttx)
                        e_ccsd = (jjx, ex + e_hf, ll.strip() + " %s" % ll.split()[-1])
                        found = True
            if found:
                eners[ix].append(e_ccsd)
        if not found:
            print(ix, px)

for ix in range(35):
    exx = sorted(eners[ix], key=lambda x: x[1])
    exv = []
    for xx in exx:
        if len(exv) == 0 or abs(xx[1] - exv[-1][1]) > 1E-8:
            exv.append(xx)
    for xx in exv[:24]:
        print('%5d %s' % (ix, xx[2]))
