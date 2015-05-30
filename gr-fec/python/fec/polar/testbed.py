__author__ = 'johannes'

import numpy as np
from encoder import PolarEncoder
from decoder import PolarDecoder


def test_enc_dec_chain():
    ntests = 100
    n = 16
    k = 8
    frozenbits = np.zeros(n - k)
    frozenbitposition = np.array((0, 1, 2, 3, 4, 5, 8, 9), dtype=int)
    for i in range(ntests):
        bits = np.random.randint(2, size=k)
        encoder = PolarEncoder(n, k, frozenbitposition, frozenbits, apply_bit_reversal=True)
        decoder = PolarDecoder(n, k, frozenbitposition, frozenbits, apply_bit_reversal=True)
        encoded = encoder.encode(bits)
        rx = decoder.decode(encoded)
        # print 'bits:', bits
        # print 'rx  :', rx
        # print (bits == rx).all()
        if not (bits == rx).all():
            raise ValueError('Test failed, input and output differ', bits, '!=', rx)



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

    test_enc_dec_chain()




if __name__ == '__main__':
    main()