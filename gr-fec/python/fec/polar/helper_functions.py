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

import numpy as np
import time, sys
import copy


def power_of_2_int(num):
    return int(np.log2(num))


def is_power_of_two(num):
    if type(num) != int:
        return False  # make sure we only compute integers.
    return num != 0 and ((num & (num - 1)) == 0)


def bit_reverse(value, n):
    # is this really missing in NumPy???
    seq = np.int(value)
    rev = np.int(0)
    rmask = np.int(1)
    lmask = np.int(2 ** (n - 1))
    for i in range(n // 2):
        shiftval = n - 1 - (i * 2)
        rshift = np.left_shift(np.bitwise_and(seq, rmask), shiftval)
        lshift = np.right_shift(np.bitwise_and(seq, lmask), shiftval)
        rev = np.bitwise_or(rev, rshift)
        rev = np.bitwise_or(rev, lshift)
        rmask = np.left_shift(rmask, 1)
        lmask = np.right_shift(lmask, 1)
    if not n % 2 == 0:
        rev = np.bitwise_or(rev, np.bitwise_and(seq, rmask))
    return rev


def bit_reverse_vector(vec, n):
    return np.array([bit_reverse(e, n) for e in vec], dtype=vec.dtype)


def unpack_byte(byte, nactive):
    if np.amin(byte) < 0 or np.amax(byte) > 255:
        return None
    if not byte.dtype == np.uint8:
        byte = byte.astype(np.uint8)
    if nactive == 0:
        return np.array([], dtype=np.uint8)
    return np.unpackbits(byte)[-nactive:]


def pack_byte(bits):
    if len(bits) == 0:
        return 0
    if np.amin(bits) < 0 or np.amax(bits) > 1:  # only '1' and '0' in bits array allowed!
        return None
    bits = np.concatenate((np.zeros(8 - len(bits), dtype=np.uint8), bits))
    res = np.packbits(bits)[0]
    return res


def show_progress_bar(ndone, ntotal):
    nchars = 50

    fract = (1. * ndone / ntotal)
    percentage = 100. * fract
    ndone_chars = int(nchars * fract)
    nundone_chars = nchars - ndone_chars
    sys.stdout.write('\r[{0}{1}] {2:5.2f}% ({3} / {4})'.format('=' * ndone_chars, ' ' * nundone_chars, percentage, ndone, ntotal))



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
    z = 0.0
    for y in range(ydim):
        z += np.sqrt(w[0, y] * w[1, y])
    # need all
    return z


def shuffle_vector(vec, bound):
    shuffle_pos = bit_reverse_vector(np.arange(bound), 4)
    return vec[shuffle_pos]


def interleave_lanes(l0, l1, bound):
    i_pos = np.reshape(np.arange(bound, dtype=int), (2, -1)).T.flatten()
    l0 = np.reshape(l0, (2, -1))
    l1 = np.reshape(l1, (2, -1))
    l0r = np.append(l0[0], l1[0])
    l1r = np.append(l0[1], l1[1])
    l0r = l0r[i_pos]
    l1r = l1r[i_pos]
    return l0r, l1r


def interleave_stages(lanes, bound):
    num_lanes = int(np.shape(lanes)[0])
    n = int(np.log2(num_lanes))

    parts = 2 ** np.arange(n)
    stages = parts[::-1]
    print(stages)
    for stage, part in zip(stages, parts):
        print(stage, part)
        for p in range(part):
            for i in range(stage):
                pos = 2 * i * part + p
                print(stage, part, p, i, pos)
                lanes[pos], lanes[pos + part] = interleave_lanes(lanes[pos], lanes[pos + part], bound)
        print(lanes)

    return lanes


def bit_reverse_vector_elements(vec):
    m = vec.size
    bound = 16
    num_lanes = m // bound
    lanes = np.zeros((num_lanes, bound), dtype=int)
    for i in range(0, num_lanes):
        p = i * bound
        part = vec[p: p + bound]
        lanes[i] = part

    print('reved lanes')
    print(lanes)

    # SHUFFLE!
    for i in range(num_lanes):
        lanes[i] = shuffle_vector(lanes[i], bound)
    print('\nshuffled lanes')
    print(lanes)

    interleave_stages(lanes, bound)
    print('\ninterleaved lanes 1')
    print(lanes)

    num_lanes = int(np.shape(lanes)[0])
    n = int(np.log2(num_lanes))
    lane_rev = bit_reverse_vector(np.arange(num_lanes, dtype=int), n)
    lanes = lanes[lane_rev]
    print('\nreversed lanes')
    print(lanes)
    return lanes.flatten()




def main():
    print 'helper functions'

    for i in range(8):
        print(i, 'is power of 2: ', is_power_of_two(i))
    n = 6
    m = 2 ** n
    k = m // 2
    eta = 0.3

    pos = np.arange(m)
    rev_pos = bit_reverse_vector(pos, n)
    print(pos)
    print(rev_pos)

    reved = bit_reverse_vector_elements(rev_pos)

    print(np.all(pos == reved))

    # bound = 16
    # shuffle_pos = bit_reverse_vector(np.arange(bound), 4)
    # for i in range(0, m, bound):
    #     print("\nlook at this part")
    #     part = pos[i: i + bound]
    #     rev = bit_reverse_vector(part, n)
    #     sorted_rev = np.sort(rev)
    #     print(part)
    #     print(rev)
    #     print(sorted_rev)
    #     sorted_part = rev[shuffle_pos]
    #     print(sorted_part)



if __name__ == '__main__':
    main()
