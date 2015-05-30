__author__ = 'johannes'


import numpy as np

'''
PolarCommon holds value checks and common initializer code for both Encoder and Decoder.
'''


class PolarCommon:
    def __init__(self, n, k, frozen_bit_position, frozenbits=None, reverse=True):
        if not self._is_power_of_two(n):
            raise ValueError("n={0} is not a power of 2!".format(n))
        if frozenbits is None:
            frozenbits = np.zeros(n - k, dtype=np.int)
        if not len(frozenbits) == n - k:
            raise ValueError("len(frozenbits)={0} is not equal to n-k={1}!".format(len(frozenbits), n - k))
        if not frozenbits.dtype == np.int:
            frozenbits = frozenbits.astype(dtype=int)
        if not len(frozen_bit_position) == (n - k):
            raise ValueError("len(frozen_bit_position)={0} is not equal to n-k={1}!".format(len(frozen_bit_position), n - k))
        if not frozen_bit_position.dtype == np.int:
            frozen_bit_position = frozen_bit_position.astype(dtype=int)

        self.frozen_bit_position_unreversed = frozen_bit_position
        if not reverse:
            numbits = np.int(np.round(np.log2(n)))
            frozen_bit_position = np.array([self._bit_reverse(i, numbits) for i in frozen_bit_position], dtype=int)
            frozen_bit_position = np.sort(frozen_bit_position)

        self.bit_reverse_positions = self._vector_bit_reversed(np.arange(n, dtype=int), int(np.log2(n)))

        self.reverse = reverse
        self.N = n
        self.K = k
        self.frozenbits = frozenbits
        self.frozen_bit_position = frozen_bit_position
        self.info_bit_position = np.delete(np.arange(self.N), self.frozen_bit_position)

        print 'frozenbits', self.frozenbits
        print 'frozenpos ', self.frozen_bit_position

    def _is_power_of_two(self, num):
        if type(num) != int:
            return False  # make sure we only compute integers.
        return num != 0 and ((num & (num - 1)) == 0)

    def _bit_reverse(self, value, n):
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

    def _vector_bit_reversed(self, vec, n):
        return np.array([self._bit_reverse(e, n) for e in vec], dtype=vec.dtype)

    def info_print(self):
        print "POLAR code ({0}, {1})".format(self.N, self.K)
