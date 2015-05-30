__author__ = 'johannes'

import numpy as np
from common import PolarCommon

# for dev
from encoder import PolarEncoder
from matplotlib import pyplot as plt


class PolarDecoder(PolarCommon):
    def __init__(self, n, k, frozen_bit_position, frozenbits=None, apply_bit_reversal=True):
        self.info_bit_position_reversed = np.delete(np.arange(n), frozen_bit_position)
        PolarCommon.__init__(self, n, k, frozen_bit_position, frozenbits, apply_bit_reversal)

        self.error_probability = 0.1  # this is kind of a dummy value. usually chosen individually.
        self.bsc_lr = ((1 - self.error_probability) / self.error_probability, self.error_probability / (1 - self.error_probability))
        self.bsc_llrs = np.log(self.bsc_lr)


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
        return (la * exp_val) + lb

    def _lr_bit(self, bit):
        return self.bsc_lr[bit]

    def _lr_odd(self, la, lb):
        # la is upper branch and lb is lower branch
        return (la * lb + 1) / (la + lb)

    def _lr_even(self, la, lb, f):
        # la is upper branch and lb is lower branch, f is last decoded bit.
        return (la ** (1 - (2 * f))) * lb

    def _lr_bit_decision(self, lr):
        if lr < 1:
            return int(1)
        return int(0)

    def _get_even_indices_values(self, u_hat):
        # looks like overkill for some indexing, but zero and one based indexing mix-up gives you haedaches.
        return u_hat[1::2]

    def _get_odd_indices_values(self, u_hat):
        return u_hat[0::2]

    def _calculate_lrs(self, y, u):
        ue = self._get_even_indices_values(u)
        uo = self._get_odd_indices_values(u)
        ya = y[0:y.size//2]
        yb = y[(y.size//2):]
        # print y, '-->', u, 'size', y.size, ya, yb
        la = self._lr_decision_element(ya, (ue + uo) % 2)
        lb = self._lr_decision_element(yb, ue)
        return la, lb

    def _lr_decision_element(self, y, u):
        if y.size == 1:
            return self._lr_bit(y[0])
        if u.size % 2 == 0:  # use odd branch formula
            la, lb = self._calculate_lrs(y, u)
            return self._lr_odd(la, lb)
        else:
            ui = u[-1]
            la, lb = self._calculate_lrs(y, u[0:-1])
            return self._lr_even(la, lb, ui)

    def _lr_sc_decoder(self, y):
        # this is the standard SC decoder as derived from the formulas. It sticks to natural bit order.
        u = np.array([], dtype=int)
        for i in range(y.size):
            lr = self._lr_decision_element(y, u)
            ui = self._retrieve_bit_from_lr(lr, i)
            u = np.append(u, ui)
        return u

    def _retrieve_bit_from_lr(self, lr, pos):
        f_index = np.where(self.frozen_bit_position == pos)[0]
        if not f_index.size == 0:
            ui = self.frozenbits[f_index]
        else:
            ui = self._lr_bit_decision(lr)
        return ui

    def decode(self, data, is_packed=False):
        if len(data) is not self.N:
            raise ValueError("len(data)={0} is not equal to n={1}!".format(len(data), self.N))
        if is_packed:
            data = np.unpackbits(data)
        data = self._lr_sc_decoder(data)
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


def test_reverse_enc_dec():
    n = 16
    k = 8
    frozenbits = np.zeros(n - k)
    frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    bits = np.random.randint(2, size=k)
    encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, apply_bit_reversal=True)
    decoder = PolarDecoder(n, k, frozenbitposition, frozenbits, apply_bit_reversal=True)
    encoded = encoder.encode(bits)
    print 'encoded:', encoded
    rx = decoder.decode(encoded)
    print 'bits:', bits
    print 'rx  :', rx
    print (bits == rx).all()



def main():
    power = 3
    n = 2 ** power
    k = 4
    frozenbits = np.zeros(n - k, dtype=int)
    frozenbitposition = np.array((0, 1, 2, 4), dtype=int)
    frozenbitposition4 = np.array((0, 1), dtype=int)


    encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, apply_bit_reversal=False)
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
    evec[1] = 0
    deced = decoder._lr_sc_decoder(evec)
    print 'SC decoded:', deced


    power = 1
    n = 2 ** power
    nodes = np.arange(n, dtype=int)
    bit_reversed = decoder._vector_bit_reversed(nodes, power)
    print bit_reversed

    graph = np.zeros((n, power + 1), dtype=int)
    print graph
    graph = traverse_graph(graph, bit_reversed)

    print "This was our graph. Now see result"
    print graph

    test_reverse_enc_dec()



if __name__ == '__main__':
    main()