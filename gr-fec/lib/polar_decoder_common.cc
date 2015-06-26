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
#include <gnuradio/fec/polar_decoder_common.h>
#include <volk/volk.h>

#include <cstdio>

namespace gr {
  namespace fec {

    polar_decoder_common::polar_decoder_common(int block_size, int num_info_bits,
                                               std::vector<int> frozen_bit_positions,
                                               std::vector<char> frozen_bit_values, bool is_packed) :
        polar_common(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values, is_packed),
        D_LLR_FACTOR(2.19722458f),
        d_frozen_bit_positions(frozen_bit_positions),
        d_frozen_bit_values(frozen_bit_values)
    {
    }

    polar_decoder_common::~polar_decoder_common()
    {
    }

    void
    polar_decoder_common::initialize_llr_vector(float* llrs, const float* input)
    {
      volk_32f_s32f_multiply_32f(llrs + block_size() * block_power(), input, D_LLR_FACTOR, block_size());
    }

    float
    polar_decoder_common::llr_odd(const float la, const float lb) const
    {
      return copysignf(1.0f, la) * copysignf(1.0f, lb) * std::min(fabs(la), fabs(lb));
    }

    float
    polar_decoder_common::llr_even(const float la, const float lb, const unsigned char f) const
    {
      switch(f){
        case 0:
          return lb + la;
        default:
          return lb - la;
      }
    }

    void
    polar_decoder_common::butterfly(float* llrs, int call_row, int stage, unsigned char* u, const int u_num)
    {
      if(!(block_power() > stage)){
        return;
      }

      const int stage_half_block_size = block_size() >> (stage + 1);
      if((call_row % (0x1 << (block_power() - stage))) >= stage_half_block_size){
        const int upper_right = call_row - stage_half_block_size;
        const unsigned char f = u[u_num - 1];
        llrs[call_row] = llr_even(llrs[block_size() + upper_right], llrs[block_size() + call_row], f);
        return;
      }

      const int lower_right = call_row + stage_half_block_size;

      unsigned char* u_half = u + block_size();

      odd_xor_even_values(u_half, u, u_num);
      butterfly(llrs + block_size(), call_row, stage + 1, u_half, u_num / 2);

      even_u_values(u_half, u, u_num);
      butterfly(llrs + block_size(), lower_right, stage + 1, u_half, u_num / 2);

      llrs[call_row] = llr_odd(llrs[block_size() + call_row], llrs[block_size() + lower_right]);
    }


    void
    polar_decoder_common::even_u_values(unsigned char* u_even, const unsigned char* u,
                                             const int u_num)
    {
      u++;
      for(int i = 1; i < u_num; i += 2){
        *u_even++ = *u;
        u += 2;
      }
    }

    void
    polar_decoder_common::odd_xor_even_values(unsigned char* u_xor, const unsigned char* u,
                                                   const int u_num)
    {
      for(int i = 1; i < u_num; i += 2){
        *u_xor++ = *u ^ *(u + 1);
        u += 2;
      }
    }

    void
    polar_decoder_common::extract_info_bits(unsigned char* output, const unsigned char* input) const
    {
      unsigned int frozenbit_num = 0;
      for(int i = 0; i < block_size(); i++){
//        std::cout << i << ": " << int(input[i]) << "\t";
        if(frozenbit_num < d_frozen_bit_positions.size() && d_frozen_bit_positions.at(frozenbit_num) == i){
//          std::cout << "frozen = TRUE!\n";
          frozenbit_num++;
        }
        else{
//          std::cout << "frozen = FALSE!\n";
          *output++ = *input;
        }
        input++;
      }
//      std::cout << std::endl;
    }

    void
    polar_decoder_common::print_pretty_llr_vector(const float* llr_vec) const
    {
      for(int row = 0; row < block_size(); row++) {
        std::cout << row << "->" << int(bit_reverse(row, block_power())) << ": ";
        for(int stage = 0; stage < block_power() + 1; stage++) {
          printf("%+4.2f, ", llr_vec[(stage * block_size()) + row]);
        }
        std::cout << std::endl;
      }
    }

  } /* namespace fec */
} /* namespace gr */

