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
Based on 2 papers:
[1] Ido Tal, Alexander Vardy: 'How To Construct Polar Codes', 2013
for an in-depth description of a widely used algorithm for channel construction.

[2] Harish Vangala, Emanuele Viterbo, Yi Hong: 'A Comparative Study of Polar Code Constructions for the AWGN Channel', 2015
for an overview of different approaches
'''


import numpy as np
from scipy.optimize import fsolve
from scipy.special import erfc
from helper_functions import *
import matplotlib.pyplot as plt


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
    return p[y == u]
    # return 1 - p if y == u else p

# @profile
def calculate_joint_transition_probability(N, y, x, transition_probs):
    single_ws = np.empty(N)
    single_ws[y == x] = transition_probs[True]
    single_ws[y != x] = transition_probs[False]
    return np.prod(single_ws)


# @profile
def w_split_prob(y, u, G, transition_probs):
    ''' Calculate channel splitting probabilities. '''
    N = len(y)  # number of combined channels
    df = N - len(u)  # degrees of freedom.
    prob = 0.0
    for uf in range(2 ** df):
        utemp = unpack_byte(np.array([uf, ]), df)
        ub = np.concatenate((u, utemp))
        x = np.dot(ub, G) % 2
        true_num = np.sum(y == x)
        false_num = N - true_num
        w = (transition_probs[True] ** true_num) * (transition_probs[False] ** false_num)
        # w = calculate_joint_transition_probability(N, y, x, transition_probs)
        prob += w
    divider = 1.0 / (2 ** (N - 1))
    return divider * prob

# @profile
def wn_split_channel(N, p):
    G = get_Gn(N)
    y = np.zeros((N, ), dtype=np.uint8)
    n = 1
    u = np.zeros((n + 1, ), dtype=np.uint8)
    transition_probs = np.array([p, 1 - p], dtype=float)

    z_params = []
    c_params = []
    for w in range(N):
        nentries = (2 ** N) * (2 ** w)
        print "for w=", w, " nentries=", nentries
        w_probs = np.zeros((nentries, 2), dtype=float)
        for y in range(2 ** N):
            yb = unpack_byte(np.array([y, ]), N)
            for u in range(2 ** (w + 1)):
                ub = unpack_byte(np.array([u, ]), w + 1)
                wp = w_split_prob(yb, ub, G, transition_probs)
                ufixed = ub[0:-1]
                yindex = y * (2 ** w) + pack_byte(ufixed)
                uindex = ub[-1]
                w_probs[yindex][uindex] = wp

        z = bhattacharyya_parameter(w_probs)
        z_params.append(z)
        c = mutual_information(w_probs)
        c_params.append(c)
        print "W=", w, "Z=", z, "capacity=", c

    return z_params, c_params


def calculate_z_param(x):
    # variables etc taken from paper bei Ido Tal et al.
    # name there is f(x)
    # x is the cross over probability of a BSC.
    return 2 * np.sqrt(x * (1 - x))


def calculate_capacity(x):
    # in paper it is called g(x)
    return -1. * x * np.log(x) - (1 - x) * np.log(1 - x)


def solver_equation(val, s):
    cw_lambda = codeword_lambda_callable(s)
    ic_lambda = instantanious_capacity_callable()
    return lambda y: ic_lambda(cw_lambda(y)) - val


def solve_capacity(a, s):
    eq = solver_equation(a, s)
    res = fsolve(eq, 1)
    return np.abs(res[0])  # only positive values needed.


def codeword_lambda_callable(s):
    return lambda y: np.exp(-2 * y * np.sqrt(2 * s))


def codeword_lambda(y, s):
    return codeword_lambda_callable(s)(y)


def instantanious_capacity_callable():
    return lambda x : 1 - np.log2(1 + x) + (x * np.log2(x) / (1 + x))


def instantanious_capacity(x):
    return instantanious_capacity_callable()(x)


def q_function(x):
    # Q(x) = (1 / sqrt(2 * pi) ) * integral (x to inf) exp(- x ^ 2 / 2) dx
    return .5 * erfc(x / np.sqrt(2))


def discretize_awgn(mu, design_snr):
    '''
    needed for Binary-AWGN channels.
    in [1] described in Section VI
    in [2] described as a function of the same name.
    in both cases reduce infinite output alphabet to a finite output alphabet of a given channel.
    idea:
    1. instantaneous capacity C(x) in interval [0, 1]
    2. split into mu intervals.
    3. find corresponding output alphabet values y of likelihood ratio function lambda(y) inserted into C(x)
    4. Calculate probability for each value given that a '0' or '1' is was transmitted.
    '''
    s = 10 ** (design_snr / 10)
    a = np.zeros(mu + 1, dtype=float)
    a[-1] = np.inf
    for i in range(1, mu):
        a[i] = solve_capacity(1. * i / mu, s)

    factor = np.sqrt(2 * s)
    tpm = np.zeros((2, mu))
    for j in range(mu):
        tpm[0][j] = q_function(factor + a[j]) - q_function(factor + a[j + 1])
        tpm[1][j] = q_function(-1. * factor + a[j]) - q_function(-1. * factor + a[j + 1])
    return tpm


def calculate_delta_I(a, b, at, bt):
    c = lambda a, b: -1. * (a + b) * np.log2((a + b) / 2) + a * np.log2(a) + b * np.log2(b)
    return c(a, b) + c(at, bt) - c(a + at, b + bt)


def quantize_to_size(tpm, mu):
    L = np.shape(tpm)[1]
    delta_i_vec = np.zeros(L - 1)
    for i in range(L - 1):
        delta_i_vec[i] = calculate_delta_I(tpm[0, i], tpm[1, i], tpm[0, i + 1], tpm[1, i + 1])

    for i in range(L - mu):
        d = np.argmin(delta_i_vec)
        ap = tpm[0, d] + tpm[0, d + 1]
        bp = tpm[1, d] + tpm[1, d + 1]
        if d > 0:
            delta_i_vec[d - 1] = calculate_delta_I(tpm[0, d - 1], tpm[1, d - 1], ap, bp)
        if d < delta_i_vec.size - 1:
            delta_i_vec[d + 1] = calculate_delta_I(ap, bp, tpm[0, d + 1], tpm[1, d + 1])
        delta_i_vec = np.delete(delta_i_vec, d)
        tpm = np.delete(tpm, d, axis=1)

        tpm[0, d] = ap
        tpm[1, d] = bp
    return tpm


def tal_vardy_tpm_algorithm(block_size, design_snr, mu):
    block_power = power_of_2_int(block_size)
    channels = np.zeros((block_size, 2, mu))
    channels[0] = discretize_awgn(mu, design_snr)

    for j in range(0, block_power):
        u = 2 ** j
        for t in range(u):
            print('u=', u, ', t=', t)
            ch1 = upper_convolve(channels[t], mu)
            ch2 = lower_convolve(channels[t], mu)
            channels[t] = quantize_to_size(ch1, mu)
            channels[u + t] = quantize_to_size(ch2, mu)


    z = np.zeros(block_size)
    for i in range(block_size):
        z[i] = np.sum(channels[i][1])
    return z



def merge_lr_based(q):
    lrs = q[0] / q[1]
    vals, indices, inv_indices = np.unique(lrs, return_index=True, return_inverse=True)
    unq_cnt = np.bincount(inv_indices)
    # compare [1] (20). Ordering of representatives according to LRs.
    temp = q[:, indices]  # so far it does no harm. Does it do, what it is supposed to do?
    temp *= unq_cnt
    return temp


def upper_convolve(tpm, mu):
    q = np.zeros((2, mu ** 2))
    idx = -1
    for i in range(mu):
        idx += 1
        q[0, idx] = (tpm[0, i] ** 2 + tpm[1, i] ** 2) / 2
        q[1, idx] = tpm[0, i] * tpm[1, i]
        for j in range(i + 1, mu):
            idx += 1
            q[0, idx] = tpm[0, i] * tpm[0, j] + tpm[1, i] * tpm[1, j]
            q[1, idx] = tpm[0, i] * tpm[1, j] + tpm[1, i] * tpm[0, j]
            if q[0, idx] < q[1, idx]:
                q[0, idx], q[1, idx] = swap_values(q[0, idx], q[1, idx])

    idx += 1
    q = np.delete(q, np.arange(idx, np.shape(q)[1]), axis=1)
    q = merge_lr_based(q)
    return q


def swap_values(first, second):
    return second, first


def lower_convolve(tpm, mu):
    q = np.zeros((2, mu * (mu + 1)))
    idx = -1
    for i in range(0, mu):
        idx += 1
        q[0, idx] = tpm[0, i] ** 2 / 2
        q[1, idx] = tpm[1, i] ** 2 / 2
        if q[0, idx] < q[1, idx]:
            q[0, idx], q[1, idx] = swap_values(q[0, idx], q[1, idx])
        idx += 1
        q[1, idx] = q[0, idx] = tpm[0, i] * tpm[1, i]

        for j in range(i + 1, mu):
            idx += 1
            q[0, idx] = tpm[0, i] * tpm[0, j]
            q[1, idx] = tpm[1, i] * tpm[1, j]
            if q[0, idx] < q[1, idx]:
                q[0, idx], q[1, idx] = swap_values(q[0, idx], q[1, idx])
            idx += 1
            q[0, idx] = tpm[0, i] * tpm[1, j]
            q[1, idx] = tpm[1, i] * tpm[0, j]
            if q[0, idx] < q[1, idx]:
                q[0, idx], q[1, idx] = swap_values(q[0, idx], q[1, idx])
    q = merge_lr_based(q)
    return q


def main():
    print 'channel construction BSC main'
    n = 6
    m = 2 ** n
    k = m // 2
    eta = 0.3
    p = 0.1

    design_snr = 0.0
    mu = 16

    z_params = tal_vardy_tpm_algorithm(m, design_snr, mu)
    plt.plot(z_params)
    plt.show()



if __name__ == '__main__':
    main()
