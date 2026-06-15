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

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'font.size': 11})
plt.figure(figsize=(9.5, 6))

mds = np.load('04-mds.npy')
meners = (np.load('04-meners.npy') - -0.3359688458) * 1000
mssqs = np.load('04-mssqs.npy')
mzs = np.load('04-mzs.npy') / 8

idx = np.argmin(np.abs(meners - 1))

fig, axs = plt.subplots(3, 1, sharex=True)

axs[0].plot(mds, meners, color='#000')
axs[1].plot(mds, mssqs - (1.5 * 2.5), color='#000')
axs[2].plot(mds, np.sum(np.abs(mzs), axis=1), color='#000')
axs[0].vlines(mds[idx], 0, 2, linestyle='--', color='#AAA')
axs[1].vlines(mds[idx], 0, 2.5, linestyle='--', color='#AAA')
axs[2].vlines(mds[idx], 0, 1.75, linestyle='--', color='#AAA')
axs[0].set_ylim((0, 2))
axs[1].set_ylim((0, 2.5))
axs[2].set_ylim((0, 1.75))

axs[0].set_ylabel('$E - E_{\\mathrm{exact}}$ (mHa)')
axs[1].set_ylabel('$|\\langle S^2 \\rangle - S(S+1)|$')
axs[2].set_ylabel('$\\sum_k |\\langle S^z_k \\rangle| / L$')
axs[2].set_xlabel('MPS bond dimension $D$')

plt.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.11, wspace=0.4, hspace=0.0)

plt.savefig('05-nsa-plot.png', dpi=600)