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

'''
[0] Erdal Arikan: 'Channel Polarization: A Method for Constructing Capacity-Achieving Codes for Symmetric Binary-Input Memoryless Channels', 2009
foundational paper for polar codes.
'''


import numpy as np
from channel_construction_bec import calculate_bec_channel_capacities, calculate_bec_channel_z_parameters
from channel_construction_bec import bhattacharyya_bounds, design_snr_to_bec_eta
from channel_construction_bsc import tal_vardy_tpm_algorithm
import matplotlib.pyplot as plt


def bsc_channel(p):
    '''
    binary symmetric channel (BSC)
    output alphabet Y = {0, 1} and
    W(0|0) = W(1|1) and W(1|0) = W(0|1)

    this function returns a prob's vector for a BSC
    p denotes an erroneous transistion
    '''
    if not (p >= 0.0 and p <= 1.0):
        print "given p is out of range!"
        return np.array([], dtype=float)

    # 0 -> 0, 0 -> 1, 1 -> 0, 1 -> 1
    W = np.array([[1 - p, p], [p, 1 - p]], dtype=float)
    return W


def bec_channel(eta):
    '''
    binary erasure channel (BEC)
    for each y e Y
    W(y|0) * W(y|1) = 0 or W(y|0) = W(y|1)
    transistions are 1 -> 1 or 0 -> 0 or {0, 1} -> ? (erased symbol)
    '''

    # looks like BSC but should be interpreted differently.
    W = np.array((1 - eta, eta, 1 - eta), dtype=float)
    return W


def mutual_information(w):
    '''
    calculate mutual information I(W)
    I(W) = sum over y e Y ( sum over x e X ( ... ) )
    .5 W(y|x) log frac { W(y|x) }{ .5 W(y|0) + .5 W(y|1) }
    '''
    ydim, xdim = np.shape(w)
    i = 0.0
    for y in range(ydim):
        for x in range(xdim):
            v = w[y][x] * np.log2(w[y][x] / (0.5 * w[y][0] + 0.5 * w[y][1]))
            i += v
    i /= 2.0
    return i


def bhattacharyya_parameter(w):
    '''bhattacharyya parameter is a measure of similarity between two prob. distributions'''
    # sum over all y e Y for sqrt( W(y|0) * W(y|1) )
    dim = np.shape(w)
    ydim = dim[0]
    xdim = dim[1]
    z = 0.0
    for y in range(ydim):
        z += np.sqrt(w[y][0] * w[y][1])
    # need all
    return z


def get_frozen_bit_indices_from_capacities(chan_caps, nfrozen):
    indexes = np.array([], dtype=int)
    while indexes.size < nfrozen:
        index = np.argmin(chan_caps)
        indexes = np.append(indexes, index)
        chan_caps[index] = 1.0
    return np.sort(indexes)


def get_frozen_bit_indices_from_z_parameters(z_params, nfrozen):
    indexes = np.array([], dtype=int)
    while indexes.size < nfrozen:
        index = np.argmax(z_params)
        indexes = np.append(indexes, index)
        z_params[index] = 0.0
    return np.sort(indexes)


def get_bec_frozen_indices(nblock, kfrozen, eta):
    bec_caps = calculate_bec_channel_capacities(eta, nblock)
    positions = get_frozen_bit_indices_from_capacities(bec_caps, kfrozen)
    return positions


def frozen_bit_positions(block_size, info_size, design_snr=0.0):
    if not design_snr > -1.5917:
        print('bad value for design_nsr, must be > > -1.5917! default=0.0')
        design_snr = 0.0
    eta = design_snr_to_bec_eta(design_snr)
    return get_bec_frozen_indices(block_size, block_size - info_size, eta)


def main():
    print 'channel construction Bhattacharyya bounds by Arikan'
    n = 8
    m = 2 ** n
    k = m // 2
    eta = 0.3
    p = 0.1
    design_snr = 2
    mu = 16

    ztv = tal_vardy_tpm_algorithm(m, design_snr, mu)
    scaling_factor = 1. / np.max(ztv)
    ztv *= scaling_factor
    plt.plot(ztv)
    # plt.plot(ztv * (1.4 * (10 ** 19)))

    z_params = calculate_bec_channel_z_parameters(eta, m)
    plt.plot(z_params)
    plt.show()

    ztv_indices = get_frozen_bit_indices_from_z_parameters(ztv, k)
    upper_indices = get_frozen_bit_indices_from_z_parameters(z_params, k)

    first = np.zeros(m)
    first[ztv_indices] = 1.0
    second = np.zeros(m)
    second[upper_indices] = 1.0
    plt.plot(first)
    plt.plot(second)
    for i in range(len(first)):
        if not first[i] == second[i]:
            plt.axvline(i, color='r')
    plt.show()


if __name__ == '__main__':
    main()


