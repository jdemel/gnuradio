/* -*- c++ -*- */
/* 
 * Copyright 2015 Free Software Foundation, Inc.
 * 
 * This file is part of GNU Radio
 * 
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 * 
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gnuradio/io_signature.h>
#include <gnuradio/fec/polar_decoder_sc.h>

namespace gr {
  namespace fec {

    generic_decoder::sptr
    polar_decoder_sc::make(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions,
                           std::vector<char> frozen_bit_values, bool is_packed)
    {
      return generic_decoder::sptr(new polar_decoder_sc(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values, is_packed));
    }

    polar_decoder_sc::polar_decoder_sc(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed):
            d_block_size(block_size),
            d_block_power((int)log2(block_size)),
            d_is_packed(is_packed),
            d_num_info_bits(num_info_bits),
            d_frozen_bit_positions(frozen_bit_positions),
            d_frozen_bit_values(frozen_bit_values)
    {
    }



    polar_decoder_sc::~polar_decoder_sc()
    {
    }

    void
    polar_decoder_sc::generic_work(void* in_buffer, void* out_buffer)
    {
    }

  } /* namespace fec */
} /* namespace gr */

