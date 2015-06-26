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
#include <gnuradio/fec/polar_decoder_sc_list.h>

#include <cmath>

namespace gr
{
  namespace fec
  {

    generic_decoder::sptr
    polar_decoder_sc_list::make(int max_list_size, int block_size, int num_info_bits,
                                std::vector<int> frozen_bit_positions,
                                std::vector<char> frozen_bit_values, bool is_packed)
    {
      return generic_decoder::sptr(
          new polar_decoder_sc_list(max_list_size, block_size, num_info_bits, frozen_bit_positions,
                                    frozen_bit_values, is_packed));
    }

    polar_decoder_sc_list::polar_decoder_sc_list(int max_list_size, int block_size,
                                                 int num_info_bits,
                                                 std::vector<int> frozen_bit_positions,
                                                 std::vector<char> frozen_bit_values,
                                                 bool is_packed) :
            polar_decoder_common(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values,
                                 is_packed), d_max_list_size(max_list_size)
    {
      calculate_path_metric(5.0, 3.0, 0);
      calculate_path_metric(5.0, 3.0, 1);
      calculate_path_metric(5.0, -3.0, 0);
      calculate_path_metric(5.0, -3.0, 1);
      for(unsigned char c = 0; c < 2; c++){
        for(float f = -3.0; f < 4.0; f += 6.0){
          std::cout << (int) c << ", " << f << "\t --> " << calculate_path_metric(5.0, f, c) << std::endl;
        }
      }
    }

    polar_decoder_sc_list::~polar_decoder_sc_list()
    {
    }

    void
    polar_decoder_sc_list::generic_work(void* in_buffer, void* out_buffer)
    {
    }

    float
    polar_decoder_sc_list::calculate_path_metric(const float last_pm, const float llr,
                                                 const unsigned char u_hat) const
    {
//      if((u_hat == 0 && llr > 0.0f) || (u_hat == 1 && llr < 0.0f)){
      if(u_hat == (unsigned char) (0.5 * 1 - copysignf(1.0f, llr))){
        return last_pm;
      }
      return last_pm + abs(llr);
    }

  } /* namespace fec */
} /* namespace gr */

