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
#include <volk/volk.h>

#include <cmath>
#include <cstdio>

namespace gr
{
  namespace fec
  {

    generic_decoder::sptr
    polar_decoder_sc::make(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions,
                           std::vector<char> frozen_bit_values, bool is_packed)
    {
      return generic_decoder::sptr(
          new polar_decoder_sc(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values,
                               is_packed));
    }

    polar_decoder_sc::polar_decoder_sc(int block_size, int num_info_bits,
                                       std::vector<int> frozen_bit_positions,
                                       std::vector<char> frozen_bit_values, bool is_packed) :
        polar_decoder_common(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values, is_packed),
//        D_LLR_FACTOR(2.19722458f),
        d_frozen_bit_counter(0)
    {
      d_llr_vec = (float*) volk_malloc(sizeof(float) * block_size * (block_power() + 1), volk_get_alignment());
      d_u_hat_vec = (unsigned char*) volk_malloc(block_size * (block_power() + 1), volk_get_alignment());
    }

    polar_decoder_sc::~polar_decoder_sc()
    {
      volk_free(d_llr_vec);
      volk_free(d_u_hat_vec);
    }

    void
    polar_decoder_sc::generic_work(void* in_buffer, void* out_buffer)
    {
      const float *in = (const float*) in_buffer;
      unsigned char *out = (unsigned char*) out_buffer;

      initialize_llr_vector(d_llr_vec, in);
      sc_decode(d_llr_vec, d_u_hat_vec);
      extract_info_bits(out, d_u_hat_vec);
    }

    void
    polar_decoder_sc::sc_decode(float* llrs, unsigned char* u)
    {
      d_frozen_bit_counter = 0;
      for(int i = 0; i < block_size(); i++){
        int row = bit_reverse(i, block_power());
        butterfly(llrs, row, 0, u, i);
        u[i] = retrieve_bit_from_llr(llrs[row], i);
      }
    }

//    void
//    polar_decoder_sc::initialize_llr_vector(float* llrs, const float* input)
//    {
//      volk_32f_s32f_multiply_32f(llrs + block_size() * block_power(), input, D_LLR_FACTOR, block_size());
//    }
//
//    void
//    polar_decoder_sc::butterfly(float* llrs, int call_row, int stage, unsigned char* u, const int u_num)
//    {
//      if(!(block_power() > stage)){
//        return;
//      }
//
//      const int stage_half_block_size = block_size() >> (stage + 1);
//      if((call_row % (0x1 << (block_power() - stage))) >= stage_half_block_size){
//        const int upper_right = call_row - stage_half_block_size;
//        const unsigned char f = u[u_num - 1];
//        llrs[call_row] = llr_even(llrs[block_size() + upper_right], llrs[block_size() + call_row], f);
//        return;
//      }
//
//      const int lower_right = call_row + stage_half_block_size;
//
//      unsigned char* u_half = u + block_size();
//
//      odd_xor_even_values(u_half, u, u_num);
//      butterfly(llrs + block_size(), call_row, stage + 1, u_half, u_num / 2);
//
//      even_u_values(u_half, u, u_num);
//      butterfly(llrs + block_size(), lower_right, stage + 1, u_half, u_num / 2);
//
//      llrs[call_row] = llr_odd(llrs[block_size() + call_row], llrs[block_size() + lower_right]);
//    }
//
//
//    void
//    gr::fec::polar_decoder_sc::even_u_values(unsigned char* u_even, const unsigned char* u,
//                                             const int u_num)
//    {
//      u++;
//      for(int i = 1; i < u_num; i += 2){
//        *u_even++ = *u;
//        u += 2;
//      }
//    }
//
//    void
//    gr::fec::polar_decoder_sc::odd_xor_even_values(unsigned char* u_xor, const unsigned char* u,
//                                                   const int u_num)
//    {
//      for(int i = 1; i < u_num; i += 2){
//        *u_xor++ = *u ^ *(u + 1);
//        u += 2;
//      }
//    }

    unsigned char
    polar_decoder_sc::retrieve_bit_from_llr(float llr, const int pos)
    {
      if(d_frozen_bit_counter < d_frozen_bit_positions.size() && pos == d_frozen_bit_positions.at(d_frozen_bit_counter)){
        return d_frozen_bit_values.at(d_frozen_bit_counter++);
      }
      return llr_bit_decision(llr);
    }
  } /* namespace fec */
} /* namespace gr */
