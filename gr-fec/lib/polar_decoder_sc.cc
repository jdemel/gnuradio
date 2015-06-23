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
        polar_common(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values, is_packed),
        D_LLR_FACTOR(2.19722458f),
            d_frozen_bit_positions(frozen_bit_positions),
            d_frozen_bit_values(frozen_bit_values)
    {
      d_llr_vec = (float*) volk_malloc(sizeof(float) * block_size * (block_power() + 1), volk_get_alignment());
    }

    polar_decoder_sc::~polar_decoder_sc()
    {
      volk_free(d_llr_vec);
    }

    void
    polar_decoder_sc::generic_work(void* in_buffer, void* out_buffer)
    {
      const float *in = (const float*) in_buffer;
      unsigned char *out = (unsigned char*) out_buffer;

      initialize_llr_vector(d_llr_vec, in);

      unsigned char* u = (unsigned char*) volk_malloc(block_size(), volk_get_alignment());
      sc_decode(d_llr_vec, u);

//      print_pretty_llr_vector();
      extract_info_bits(out, u);

      volk_free(u);
    }

    void
    polar_decoder_sc::sc_decode(float* llrs, unsigned char* u)
    {
//      d_node_counter = 0;
      d_frozen_bit_counter = 0;
      for(int i = 0; i < block_size(); i++){
        int row = bit_reverse(i, block_power());
        butterfly(llrs, row, 0, u, i);
        u[i] = retrieve_bit_from_llr(llrs[row], i);
//        print_pretty_llr_vector();
      }
    }

    float
    polar_decoder_sc::llr_odd(float la, float lb) const
    {
      return copysignf(1.0f, la) * copysignf(1.0f, lb) * std::min(fabs(la), fabs(lb));
    }

    float
    polar_decoder_sc::llr_even(float la, float lb, unsigned char f) const
    {
      switch(f){
        case 0:
          return lb + la;
        default:
          return lb - la;
      }
    }

    void
    polar_decoder_sc::initialize_llr_vector(float* llrs, const float* input)
    {
      volk_32f_s32f_multiply_32f(llrs + block_size() * block_power(), input, D_LLR_FACTOR, block_size());
    }

    void
    polar_decoder_sc::butterfly(float* llrs, int call_row, int stage, const unsigned char* u, const int u_num)
    {
//      for(int i = 0; i < stage; i++){std::cout << "\t";}
//      std::cout << d_node_counter + 1 <<  "\tcall_row: " << call_row << "\tstage: " << stage << "\t"<< u_num << "\t";
//      for(int i = 0; i < u_num; i++){
//        std::cout << int(u[i]) << ", ";
//      }
//      std::cout << std::endl;
//      d_node_counter++;
      if(!(block_power() > stage)){
//        std::cout << "call_row: " << call_row << "\tstage: " << stage << "\tRETURN" << std::endl;
        return;
      }
      const int stage_offset = block_size() * stage;
      const int next_offset = block_size() * (stage + 1);

      const int stage_half_block_size = block_size() >> (stage + 1);
      if((call_row % (0x1 << (block_power() - stage))) >= stage_half_block_size){
        const int upper_right = call_row - stage_half_block_size;
        const unsigned char f = u[u_num - 1];
        llrs[stage_offset + call_row] = llr_even(llrs[next_offset + upper_right], llrs[next_offset + call_row], f);
        return;
      }

      const int lower_right = call_row + stage_half_block_size;

      unsigned char* u_half = (unsigned char*) volk_malloc(u_num / 2, volk_get_alignment());
      odd_xor_even_values(u_half, u, u_num);

      butterfly(llrs, call_row, stage + 1, u_half, u_num / 2);
      even_u_values(u_half, u, u_num);
      butterfly(llrs, lower_right, stage + 1, u_half, u_num / 2);
      volk_free(u_half);

      llrs[stage_offset + call_row] = llr_odd(llrs[next_offset + call_row], llrs[next_offset + lower_right]);
    }


    void
    gr::fec::polar_decoder_sc::even_u_values(unsigned char* u_even, const unsigned char* u,
                                             const int u_num)
    {
      int r = 0;
      for(int i = 1; i < u_num; i += 2, r++){
        u_even[r] = u[i];
      }
    }

    void
    gr::fec::polar_decoder_sc::odd_xor_even_values(unsigned char* u_xor, const unsigned char* u,
                                                   const int u_num)
    {
      int r = 0;
      for(int i = 1; i < u_num; i += 2){
        u_xor[r] = u[i] ^ u[i - 1];
        r++;
      }
    }

    void
    polar_decoder_sc::print_pretty_llr_vector() const
    {
      for(int row = 0; row < block_size(); row++){
        std::cout << row << "->" << int(bit_reverse(row, block_power())) << ": ";
        for(int stage = 0; stage < block_power() + 1; stage++){
          printf("%+4.2f, ", d_llr_vec[(stage * block_size()) + row]);
        }
        std::cout << std::endl;
      }
    }


    unsigned char
    polar_decoder_sc::retrieve_bit_from_llr(float llr, const int pos)
    {
      if(d_frozen_bit_counter < d_frozen_bit_positions.size() && pos == d_frozen_bit_positions.at(d_frozen_bit_counter)){
        return d_frozen_bit_values.at(d_frozen_bit_counter++);
      }
      return llr_bit_decision(llr);
    }

    void
    polar_decoder_sc::extract_info_bits(unsigned char* output, const unsigned char* input)
    {
      unsigned int frozenbit_num = 0;
      int infobit_num = 0;
      for(int i = 0; i < block_size(); i++){
//        std::cout << i << ": " << int(input[i]) << "\t";
        if(frozenbit_num < d_frozen_bit_positions.size() && d_frozen_bit_positions.at(frozenbit_num) == i){
//          std::cout << "frozen = TRUE!\n";
          frozenbit_num++;
          continue;
        }
        else{
//          std::cout << "frozen = FALSE!\n";
          *(output + infobit_num) = *(input + i);
          infobit_num++;
        }
      }
//      std::cout << std::endl;
    }

  } /* namespace fec */
} /* namespace gr */
