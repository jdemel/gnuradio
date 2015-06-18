#!/usr/bin/env python
#
# Copyright 2015 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
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

from gnuradio import gr, gr_unittest, blocks
import fec_swig as fec
from _qa_helper import _qa_helper
import numpy as np

from extended_encoder import extended_encoder
from extended_decoder import extended_decoder
from polar.encoder import PolarEncoder
from polar.decoder import PolarDecoder
from polar.helper_functions import get_frozen_bit_positions
# from polar.helper_functions import bit_reverse_vector


class test_polar_decoder_sc(gr_unittest.TestCase):

    def setUp(self):
        self.tb = gr.top_block()

    def tearDown(self):
        self.tb = None

    def test_001_setup(self):
        is_packed = False
        block_size = 16
        num_info_bits = 8
        frozen_bit_positions = np.arange(block_size - num_info_bits)
        frozen_bit_values = np.array([],)

        polar_decoder = fec.polar_decoder_sc.make(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values)

        self.assertEqual(block_size, polar_decoder.get_output_size())
        self.assertEqual(num_info_bits, polar_decoder.get_input_size())
        self.assertFloatTuplesAlmostEqual((float(num_info_bits) / block_size, ), (polar_decoder.rate(), ))
        self.assertFalse(polar_decoder.set_frame_size(10))

    def test_002_one_vector(self):
        print "test"
        is_packed = False
        block_size = 16
        num_info_bits = 8
        num_frozen_bits = block_size - num_info_bits
        frozen_bit_positions = get_frozen_bit_positions('polar', block_size, num_frozen_bits, 0.11)
        frozen_bit_values = np.array([0] * num_frozen_bits,)

        python_decoder = PolarDecoder(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values)

        data = np.ones(block_size, dtype=int)

        polar_decoder = fec.polar_decoder_sc.make(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values, is_packed)
        src = blocks.vector_source_f(data, False)
        dec_block = extended_decoder(polar_decoder, None)
        snk = blocks.vector_sink_b(1)

        self.tb.connect(src, dec_block)
        self.tb.connect(dec_block, snk)
        self.tb.run()

        res = np.array(snk.data()).astype(dtype=int)

        ref = python_decoder.decode(data)

        print(res)
        print(ref)


if __name__ == '__main__':
    gr_unittest.run(test_polar_decoder_sc)


