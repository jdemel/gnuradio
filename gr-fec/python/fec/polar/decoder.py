__author__ = 'johannes'

import numpy as np
from common import PolarCommon

# for dev
from encoder import PolarEncoder
from matplotlib import pyplot as plt


class PolarDecoder(PolarCommon):
    def __init__(self, n, k, frozen_bit_position, frozenbits=None, reverse=True):
        PolarCommon.__init__(self, n, k, frozen_bit_position, frozenbits, reverse)

        self.error_probability = 0.1
        bsc_error_probs = ((1 - self.error_probability) / self.error_probability, self.error_probability / (1 - self.error_probability))
        self.bsc_llrs = np.log(bsc_error_probs)

    def _bit_llr(self, bit):
        if type(bit) == np.ndarray:
            bit = bit[0]
        llr = self.bsc_llrs[bit]
        print "bit_llr: bit=", bit, "llr=", llr
        return llr

    def _llr_odd_func(self, la, lb):
        return np.sign(la) * np.sign(lb) * np.minimum(np.abs(la), np.abs(lb))

    def _llr_even_func(self, la, lb, f):
        exp_val = float(1 - 2 * f)
        return (la ** exp_val) * lb

    def _llr(self, i_in, y, u_hat):
        # i_in is one based for now. Just like in the paper.
        if np.size(y) == 1:
            return self._bit_llr(y)
        i = (i_in + 1) // 2
        N = np.size(y) // 2
        u_odd = u_hat[np.arange(0, 2 * i - 2, 2)]
        u_even = u_hat[np.arange(1, 2 * i - 2, 2)]
        la = self._llr(i, y[0:N], (u_odd + u_even) % 2)
        lb = self._llr(i, y[N:], u_even)
        if i_in % 2 == 0:
            eo = "even"
        else:
            eo = "odd "
        print "llr", eo, ":", i_in, y, u_hat, la, lb
        if i_in % 2 == 0:
            return self._llr_even_func(la, lb, u_hat[-1])
        else:
            return self._llr_odd_func(la, lb)

    def _decode_info_bit(self, dv, data, pos):
        print "decode_info_bit:", pos + 1
        llr = self._llr(pos + 1, data, dv)
        print "decode_info_bit with LLR=", llr
        if llr >= 0:
            return 0
        else:
            return 1
        return 4

    def _decode_bit(self, dv, data, pos):
        f_index = np.where(self.frozen_bit_position == pos)[0]
        if not np.shape(f_index)[0] == 0:
            bit = self.frozenbits[f_index[0]]
        else:
            bit = self._decode_info_bit(dv, data, pos)
        return bit

    def _decode_sc(self, data):
        print data
        print self.frozenbits
        print self.frozen_bit_position
        res = np.array([], dtype=int)
        for i in range(len(data)):
            print "bitnum=", i, "y=", data, "u_hat=", res
            bit = self._decode_bit(res, data, i)
            res = np.append(res, bit)
        return res

    def _extract_info_bits(self, data):
        return data[self.info_bit_position]

    def decode(self, data, is_packed=False):
        if len(data) is not self.N:
            raise ValueError("len(data)={0} is not equal to n={1}!".format(len(data), self.N))
        if is_packed:
            data = np.unpackbits(data)
        data = self._decode_sc(data)
        data = self._extract_info_bits(data)
        if is_packed:
            data = np.packbits(data)
        return data


def traverse_graph(graph, bit_reversed):
    num = 1
    for r in bit_reversed:
        graph, num = handle_node(graph, r, 0, num)
    return graph


def handle_node(graph, r, e, num):
    if not graph[r][e] == 0:
        # print "found recursion break", r, e, num
        return graph, num
    p = np.shape(graph)[1] - 1
    nrow = np.shape(graph)[0]
    n = 2 ** p
    r1 = r
    graph[r][e] = num
    num += 1
    # print graph

    if e < p:
        graph, num = handle_node(graph, r1, e + 1, num)
        rj = r + 2 ** (p - e - 1)
        if rj < nrow:
            graph, num = handle_node(graph, rj, e + 1, num)
    return graph, num


def main():
    power = 3
    n = 2 ** power
    k = 4
    frozenbits = np.zeros(n - k, dtype=int)
    frozenbitposition = np.array((0, 1, 2, 4), dtype=int)


    encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, reverse=False)
    decoder = PolarDecoder(n, k, frozenbitposition, frozenbits)

    bits = np.ones(k, dtype=int)
    # bits = np.array([1, 0, 1, 0], dtype=int)
    print "bits: ", bits
    evec = encoder.encode(bits)
    print "froz: ", encoder._insert_frozen_bits(bits)
    print "evec: ", evec
    # dvec = decoder.decode(evec)
    # print "dec:  ", dvec

    # llr = decoder._llr(4, evec, np.array([0, 0, 0]))
    # print "llr=", llr


    power = 2
    n = 2 ** power
    nodes = np.arange(n, dtype=int)
    bit_reversed = decoder._vector_bit_reversed(nodes, power)
    print bit_reversed

    graph = np.zeros((n, power + 1), dtype=int)
    print graph
    graph = traverse_graph(graph, bit_reversed)

    print "This was our graph. Now see result"
    print graph




if __name__ == '__main__':
    main()