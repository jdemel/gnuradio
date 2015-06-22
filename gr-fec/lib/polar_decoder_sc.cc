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
      d_llr_vec = (float*) volk_malloc(sizeof(float) * block_size, volk_get_alignment());
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
      std::cout << "init'ed!\n";
      for(int i = 0; i < block_size(); i++){
        std::cout << d_llr_vec[i] << ", ";
      }
      std::cout << std::endl;

      unsigned char* u = (unsigned char*) volk_malloc(block_size(), volk_get_alignment());
      sc_decode(d_llr_vec, u);
      std::cout << "first llr: " << d_llr_vec[0] << std::endl;
      for(int i = 0; i < block_size(); i++){
        std::cout << d_llr_vec[i] << ", ";
      }
      std::cout << std::endl;
      extract_info_bits(out, u);

//      for(int i = 0; i < num_info_bits(); i++){
//        *out++ = llr_bit_decision(*in++);
//      }
      std::cout << "work DONE!\n";
//      memcpy(out, in, num_info_bits());
      volk_free(u);
    }

    void
    polar_decoder_sc::sc_decode(float* llrs, unsigned char* u)
    {
      d_node_counter = 0;
      for(int i = 0; i < block_size(); i++){
        int row = bit_reverse(i, block_power());
        butterfly(llrs, row, 0, u, i);
        u[i] = llr_bit_decision(llrs[row]);
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
      volk_32f_s32f_multiply_32f(llrs, input, D_LLR_FACTOR, block_size());
    }

    void
    polar_decoder_sc::butterfly(float* llrs, int call_row, int stage, const unsigned char* u, const int u_num)
    {
      for(int i = 0; i < stage; i++){std::cout << "\t";}
      std::cout << d_node_counter <<  "\tcall_row: " << call_row << "\tstage: " << stage << "\t";
      for(int i = 0; i < u_num; i++){
        std::cout << int(u[i]) << ", ";
      }
      std::cout << std::endl;
      d_node_counter++;
      if(!(block_power() > stage)){
//        std::cout << "call_row: " << call_row << "\tstage: " << stage << "\tRETURN" << std::endl;
        return;
      }



      const int stage_half_block_size = block_size() >> (stage + 1);
      if((call_row % (0x1 << (block_power() - stage))) >= stage_half_block_size){
        const int upper_right = call_row - stage_half_block_size;
        const unsigned char f = u[u_num - 1];
        llrs[call_row] = llr_even(llrs[upper_right], llrs[call_row], f);
        return;
      }

      const int lower_right = call_row + stage_half_block_size;

      unsigned char* u_half = (unsigned char*) volk_malloc(u_num / 2, volk_get_alignment());
      odd_xor_even_values(u_half, u, u_num);

      butterfly(llrs, call_row, stage + 1, u, u_num / 2);
      even_u_values(u_half, u, u_num);
      butterfly(llrs, lower_right, stage + 1, u, u_num / 2);
      volk_free(u_half);

//      float la = llrs[call_row];
//      float lb = llrs[lower_right];
      llrs[call_row] = llr_odd(llrs[call_row], llrs[lower_right]);
//      std::cout << "call_row: " << call_row << "\tstage: " << stage << "\tla: " << la << "\tlb: " << lb << std::endl;
    }


    void
    gr::fec::polar_decoder_sc::even_u_values(unsigned char* u_even, const unsigned char* u,
                                             const int u_num)
    {
//      std::cout << u_num << "\teven vals: ";
      int r = 0;
      for(int i = 1; i < u_num; i += 2, r++){
        u_even[r] = u[i];
//        std::cout << int(u_even[r]) << ", ";
      }
//      std::cout << std::endl;
    }

    void
    gr::fec::polar_decoder_sc::odd_xor_even_values(unsigned char* u_xor, const unsigned char* u,
                                                   const int u_num)
    {
//      std::cout << u_num << "\txor  vals: ";
      int r = 0;
      for(int i = 1; i < u_num; i += 2, r++){
        u_xor[r] = u[i] ^ u[i - 1];
//        std::cout << int(u_xor[r]) << ", ";
      }
//      std::cout << std::endl;
    }

    void
    polar_decoder_sc::extract_info_bits(unsigned char* output, const unsigned char* input)
    {
      unsigned int frozenbit_num = 0;
      int infobit_num = 0;
      for(int i = 0; i < block_size(); i++){
        std::cout << int(input[i]) << ", ";
        if(frozenbit_num < d_frozen_bit_positions.size() && d_frozen_bit_positions.at(frozenbit_num) == i){
          frozenbit_num++;
          continue;
        }
        else{
          *(output + infobit_num) = *(input + i);
        }
      }
      std::cout << std::endl;
    }

  } /* namespace fec */
} /* namespace gr */
