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


for iiid in [25, 26, 27]:
    data = []
    for ix in range(24):
        with open('00-fm-ccsdt-%d-%d.out' % (iiid, ix), 'r') as f:
            ll = f.readlines()
            if '::' in ll[-1]:
                data.append((float(ll[-1].split()[-2]), ix, ll[-1].strip(), 1))
            else:
                e_hf = 0
                e_corr = 0
                for l in ll:
                    if l.startswith('EHF  ='):
                        e_hf = float(l.split()[2])
                    elif 'E_corr(UCCSDT) =' in l:
                        e_corr = float(l.split('E_corr(UCCSDT) =')[1].split()[0])
                data.append((e_hf + e_corr, ix, ll[-1].strip(), -1))
    data.sort()
    for p in data[:1]:
        print(iiid, *p[1:])