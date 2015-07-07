#!/usr/bin/env python
#
# Copyright 2015 Free Software Foundation, Inc.
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
#

from helper_functions import *
import matplotlib.pyplot as plt


def bhattacharyya_bounds(block_size, design_snr):
    '''
    Harish Vangala, Emanuele Viterbo, Yi Hong: 'A Comparative Study of Polar Code Constructions for the AWGN Channel', 2015
    In this paper it is called Bhattacharyya bounds channel construction and is abbreviated PCC-0
    Best design SNR for block_size = 2048, R = 0.5, is 0dB.
    Compare with Arikan: 'Channel Polarization: A Method for Constructing Capacity-Achieving Codes for Symmetric Binary-Input Memoryless Channels.
    Proposition 5. inequalities turn into equalities for BEC channel. Otherwise they represent an upper bound.
    :return Z-parameters in natural bit-order. Choose according to desired rate.
    '''
    # minimum design snr = -1.5917 corresponds to BER = 0.5
    block_power = power_of_2_int(block_size)
    s = 10 ** (design_snr / 10)  # 'initial z parameter'.
    print(block_size, block_power, design_snr, s)
    z_params = np.zeros(block_size, dtype=float)
    z_params[0] = np.exp(-s)

    for j in range(block_power):
        u = 2 ** j
        for t in range(u):
            z_val = z_params[t]
            z_params[t] = (2 * z_val) - (z_val ** 2)
            z_params[u / 2 + t] = z_val ** 2
    # algorithm operates in bit reversed order. Rewind to natural bit-order.
    return z_params[bit_reverse_vector(np.arange(block_size), block_power)]


def main():
    print 'channel construction Bhattacharyya bounds by Arikan'
    n = 10
    m = 2 ** n
    k = m // 2
    eta = 0.3
    p = 0.1
    design_snr = -1.0

    z_params = bhattacharyya_bounds(m, design_snr)
    plt.plot(z_params)
    plt.show()



if __name__ == '__main__':
    main()

