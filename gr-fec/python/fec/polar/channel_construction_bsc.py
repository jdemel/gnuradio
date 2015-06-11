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
from helper_functions import *
import matplotlib.pyplot as plt


# def bit_reverse(value, n):
#     # is this really missing in NumPy???
#     bits = np.zeros(n, type(value))
#     for index in range(n):
#         mask = 1
#         mask = np.left_shift(mask, index)
#         bit = np.bitwise_and(value, mask)
#         bit = np.right_shift(bit, index)
#         bits[index] = bit
#     bits = bits[::-1]
#     result = 0
#     for index, bit in enumerate(bits):
#         bit = np.left_shift(bit, index)
#         result += bit
#     return result


def get_Bn(n):
    # this is a bit reversal matrix.
    lw = int(np.log2(n))  # number of used bits
    indexes = [bit_reverse(i, lw) for i in range(n)]
    Bn = np.zeros((n, n), type(n))
    for i, index in enumerate(indexes):
        Bn[i][index] = 1
    return Bn


def get_Fn(n):
    # this matrix defines the actual channel combining.
    if n == 1:
        return np.array([1, ])
    F2 = np.array([[1, 0], [1, 1]], np.int)
    nump = int(np.log2(n)) - 1  # number of Kronecker products to calculate
    Fn = F2
    for i in range(nump):
        Fn = np.kron(Fn, F2)
    return Fn

def get_Gn(n):
    # this matrix is called generator matrix
    if not is_power_of_two(n):
        print "invalid input"
        return None
    if n == 1:
        return np.array([1, ])
    Bn = get_Bn(n)
    Fn = get_Fn(n)
    Gn = np.dot(Bn, Fn)
    return Gn


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


def w_transition_prob(y, u, p):
    return 1 - p if y == u else p


def w_split_prob(y, u, G, p):
    ''' Calculate channel splitting probabilities. '''
    N = len(y)  # number of combined channels
    df = N - len(u)  # degrees of freedom.
    prob = 0.0
    for uf in range(2 ** df):
        utemp = unpack_byte(np.array([uf, ]), df)
        ub = np.concatenate((u, utemp))
        # print "y=", y, "\tu=", ub
        x = np.dot(ub, G) % 2
        # print "x=", x
        w = 1.0
        for i in range(N):
            w *= w_transition_prob(y[i], x[i], p)
        # print "for u1=", uf, "prob=", w
        prob += w
    divider = 1.0 / (2 ** (N - 1))
    return divider * prob


def wn_split_channel(N, p):
    G = get_Gn(N)
    y = np.zeros((N, ), dtype=np.uint8)
    n = 1
    u = np.zeros((n + 1, ), dtype=np.uint8)

    pn = w_split_prob(y, u, G, p)
    # print "pn=", pn
    z_params = []
    c_params = []
    for w in range(N):
        nentries = (2 ** N) * (2 ** w)
        print "for w=", w, " nentries=", nentries
        w_probs = np.zeros((nentries, 2), dtype=float)
        for y in range(2 ** N):
            yb = unpack_byte(np.array([y, ]), N)
            # print "\n\nyb", yb
            for u in range(2 ** (w + 1)):
                ub = unpack_byte(np.array([u, ]), w + 1)
                # print "ub", ub
                wp = w_split_prob(yb, ub, G, p)
                ufixed = ub[0:-1]
                yindex = y * (2 ** w) + pack_byte(ufixed)
                uindex = ub[-1]
                # print "y=", y, "ub", u, " prob=", wp, "index=[", yindex, ", ", uindex, "]"
                w_probs[yindex][uindex] = wp
        print "\n", np.sum(w_probs, axis=0)
        z = bhattacharyya_parameter(w_probs)
        z_params.append(z)
        c = mutual_information(w_probs)
        c_params.append(c)
        print "W=", w, "Z=", z, "capacity=", c

    return z_params, c_params


def main():
    print 'channel construction BSC main'
    n = 2 ** 6
    k = n // 2
    eta = 0.3
    p = 0.1
    z, c = wn_split_channel(4, p)
    print(z)
    print(c)


if __name__ == '__main__':
    main()
