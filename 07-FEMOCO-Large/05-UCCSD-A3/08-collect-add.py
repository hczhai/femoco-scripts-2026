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

fid = 7

for ix in range(20):
    if os.path.exists('%02d-uccsd-act-add-%d.out' % (fid, ix)):
        e_hf = 0.0
        e_ccsd = None
        with open('%02d-uccsd-act-add-%d.out' % (fid, ix), 'r') as f:
            for l in f.readlines():
                if l.startswith('--- '):
                    iix, jjx, ttx = int(l.split()[1]), int(l.split()[3]), l.split()[4]
                elif l.startswith('xdstr ='):
                    xdstr = l.split()[2]
                elif 'converged SCF energy =' in l:
                    e_hf = float(l.split('converged SCF energy =')[1].split()[0])
                elif 'E_corr(UCCSD) =' in l:
                    ex = float(l.split('E_corr(UCCSD) =')[1].split()[0])
                    ll = ":: %5d %s %20.12f %20.12f %20.12f %d %s" % (jjx, xdstr, e_hf, ex + e_hf, 0.0, -3, ttx)
                    e_ccsd = ll.strip() + " %s" % ll.split()[-1]
        print('%5d %s' % (iix, e_ccsd))
