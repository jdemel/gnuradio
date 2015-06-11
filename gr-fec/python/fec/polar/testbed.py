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
from encoder import PolarEncoder
from decoder import PolarDecoder

import matplotlib.pyplot as plt


def test_enc_dec_chain():
    ntests = 1000
    n = 32
    k = 16
    frozenbits = np.zeros(n - k)
    frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24), dtype=int)
    # frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    for i in range(ntests):
        bits = np.random.randint(2, size=k)
        encoder = PolarEncoder(n, k, frozenbitposition, frozenbits)
        decoder = PolarDecoder(n, k, frozenbitposition, frozenbits)
        encoded = encoder.encode(bits)
        rx = decoder.decode(encoded)
        if not is_equal(bits, rx):
            raise ValueError('Test #', i, 'failed, input and output differ', bits, '!=', rx)
            return


def is_equal(first, second):
    if not (first == second).all():
        result = first == second
        for i in range(len(result)):
            print '{0:4}: {1:2} == {2:1} = {3}'.format(i, first[i], second[i], result[i])
        return False
    return True


def exact_value(la, lb):
    return np.log((np.exp(la + lb) + 1) / (np.exp(la) + np.exp(lb)))


def approx_value(la, lb):
    return np.sign(la) * np.sign(lb) * np.minimum(np.abs(la), np.abs(lb))


def main():
    # n = 16
    # k = 8
    # frozenbits = np.zeros(n - k)
    # frozenbitposition8 = np.array((0, 1, 2, 4), dtype=int)
    # frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    # print frozenbitposition

    test_enc_dec_chain()


if __name__ == '__main__':
    main()