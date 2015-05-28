__author__ = 'johannes'

import numpy as np
from encoder import PolarEncoder
from decoder import PolarDecoder


def test_enc_dec_chain():
    ntests = 1
    n = 16
    k = 8
    frozenbits = np.zeros(n - k)
    frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    for i in range(ntests):
        bits = np.random.randint(2, size=k)
        encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, reverse=True)
        # encoder_f = PolarEncoder(n, k, frozenbitposition, frozenbits, reverse=False)
        decoder = PolarDecoder(n, k, frozenbitposition, frozenbits)
        encoded = encoder.encode(bits)
        rx = decoder.decode(encoded)
        print 'bits:', bits
        print 'rx  :', rx
        print (bits == rx).all()


def main():
    # n = 16
    # k = 8
    # frozenbits = np.zeros(n - k)
    # frozenbitposition8 = np.array((0, 1, 2, 4), dtype=int)
    # frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    # print frozenbitposition
    #
    # encoder = PolarEncoder(n, k, frozenbitposition, frozenbits)
    # encoder_f = PolarEncoder(n, k, frozenbitposition, frozenbits, reverse=False)
    # decoder = PolarDecoder(n, k, frozenbitposition, frozenbits)
    #
    # data = np.ones(k)
    # print encoder.encode(data)
    # print encoder_f.encode(data)
    #
    # # bits = np.ones(k)
    # # enc_vec = encoder.encode(bits)
    # #
    # # print enc_vec
    # #
    # # dec_vec = decoder.decode(enc_vec)
    # # print dec_vec
    #
    # n = 10
    # N = 2 ** n
    # print "for exp=", n, "block length=", N
    # seq = 128
    # # seq = 1023
    # print np.binary_repr(seq, n)
    # print encoder._bit_reverse(8, n)
    #
    # reversed = np.int(0)
    # rmask = np.int(1)
    # lmask = np.int(2 ** (n - 1))
    # for i in range(n // 2):
    #     shiftval = n - 1 - (i * 2)
    #     rshift = np.left_shift(np.bitwise_and(seq, rmask), shiftval)
    #     lshift = np.right_shift(np.bitwise_and(seq, lmask), shiftval)
    #     reversed = np.bitwise_or(reversed, rshift)
    #     reversed = np.bitwise_or(reversed, lshift)
    #     rmask = np.left_shift(rmask, 1)
    #     lmask = np.right_shift(lmask, 1)
    #
    # print "reved=", np.binary_repr(reversed, n)
    test_enc_dec_chain()




if __name__ == '__main__':
    main()