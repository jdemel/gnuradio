__author__ = 'johannes'

import numpy as np
from common import PolarCommon


class PolarEncoder(PolarCommon):
    def __init__(self, n, k, frozen_bit_position, frozenbits=None, reverse=True):
        PolarCommon.__init__(self, n, k, frozen_bit_position, frozenbits, reverse)
        self.G = self._gn(n, reverse)

    def _gn(self, n, reverse=True):
        # this matrix is called generator matrix
        if n == 1:
            return np.array([1, ])
        b = self._bn(n, reverse)
        f = self._fn(n)
        g = np.dot(b, f)
        return g

    def _bn(self, n, reverse=True):
        # this is a bit reversal matrix.
        if reverse is False:
            return np.identity(n, dtype=int)
        lw = int(np.log2(n))  # number of used bits
        indexes = [self._bit_reverse(i, lw) for i in range(n)]
        bn = np.zeros((n, n), type(n))
        for i, index in enumerate(indexes):
            bn[i][index] = 1
        return bn

    def _fn(self, n):
        # this matrix defines the actual channel combining.
        if n == 1:
            return np.array([1, ])
        f2 = np.array([[1, 0], [1, 1]], np.int)
        nump = int(np.log2(n)) - 1  # number of Kronecker products to calculate
        fn = f2
        for i in range(nump):
            fn = np.kron(fn, f2)
        return fn

    def get_gn(self):
        return self.G

    def _insert_frozen_bits(self, data):
        prototype = np.empty(self.N, dtype=int)
        prototype[self.frozen_bit_position] = self.frozenbits
        prototype[self.info_bit_position] = data
        return prototype

    def _encode_matrix(self, data):
        data = self._insert_frozen_bits(data)
        data = np.dot(data, self.G) % 2
        data = data.astype(dtype=int)
        return data

    def encode(self, data, is_packed=False):
        if len(data) is not self.K:
            raise ValueError("len(data)={0} is not equal to k={1}!".format(len(data), self.K))
        if is_packed:
            data = np.unpackbits(data)
        if np.max(data) > 1 or np.min(data) < 0:
            raise ValueError("can only encode bits!")
        data = self._encode_matrix(data)
        if is_packed:
            data = np.packbits(data)
        return data

    def _encode_efficient(self, vec):
        vec = self._insert_frozen_bits(vec)
        nstages = int(np.log2(self.N))

        if self.reverse:
            vec = vec[self.bit_reverse_positions]

        pos = np.arange(self.N, dtype=int)
        for i in range(nstages):
            splitted = np.reshape(pos, (2 ** (i + 1), -1))
            upper_branch = splitted[0::2].flatten()
            lower_branch = splitted[1::2].flatten()
            vec[upper_branch] = (vec[upper_branch] + vec[lower_branch]) % 2
        return vec


def compare_results(encoder, ntests, k):
    for n in range(ntests):
        bits = np.random.randint(2, size=k)
        menc = encoder.encode(bits)
        fenc = encoder._encode_efficient(bits)
        if (menc == fenc).all() == False:
            return False
    return True


def test_encoder_impls():
    print('comparing encoder implementations')
    ntests = 1000
    n = 16
    k = 8
    frozenbits = np.zeros(n - k)
    # frozenbitposition8 = np.array((0, 1, 2, 4), dtype=int)  # keep it!
    frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, reverse=True)
    print 'reverse=True , result:', compare_results(encoder, ntests, k)

    encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, reverse=False)
    print 'reverse=False, result:', compare_results(encoder, ntests, k)


def main():
    print "main in encoder"
    test_encoder_impls()


if __name__ == '__main__':
    main()