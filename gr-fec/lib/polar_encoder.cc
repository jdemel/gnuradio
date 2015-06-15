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
#include <gnuradio/fec/polar_encoder.h>
#include <cmath>
#include <stdexcept>
#include <volk/volk.h>

namespace gr {
  namespace fec {

    generic_encoder::sptr
    polar_encoder::make(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions,
                        std::vector<char> frozen_bit_values)
    {
      return generic_encoder::sptr(new polar_encoder(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values));
    }

    polar_encoder::polar_encoder(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values):
        d_block_size(block_size),
        d_num_info_bits(num_info_bits),
        d_frozen_bit_positions(frozen_bit_positions),
        d_frozen_bit_values(frozen_bit_values)
    {
      int n = (int) log2(float(block_size));
      d_block_power = n;
      if(pow(2, n) != block_size){
        throw std::runtime_error("block_size MUST be a power of 2!");
      }
      unsigned int num_frozen_bits = d_block_size - d_num_info_bits;
      if(num_frozen_bits != d_frozen_bit_positions.size()){
        std::cout << num_frozen_bits << "num pos:" << d_frozen_bit_positions.size();
        throw std::runtime_error("number of frozen bit positions must equal block_size - num_info_bits");
      }

      while(d_frozen_bit_values.size() < num_frozen_bits){
        d_frozen_bit_values.push_back(0);
      }

//      for(unsigned int i = 0; i < num_frozen_bits; i++){
//        std::cout << "position " << d_frozen_bit_positions[i] << ", value " << int(d_frozen_bit_values[i]) << std::endl;
//      }

    }

    polar_encoder::~polar_encoder()
    {
    }

    void
    polar_encoder::generic_work(void* in_buffer, void* out_buffer)
    {
      const char *in = (const char *) in_buffer;
      char *out = (char *) out_buffer;
      char* uncoded_arr = (char*) volk_malloc(sizeof(char) * d_block_size, volk_get_alignment());
      insert_frozen_bits(uncoded_arr, in);
      bit_reverse_vector(out, uncoded_arr);
      encode_vector(out);

//      for(int i = 0; i < d_block_size; i++){
//        out[i] = bit_reverse(long(in[i]), d_block_power);
//      }
    }

    void
    polar_encoder::bit_reverse_vector(char* target, const char* input)
    {
      for(int i = 0; i < d_block_size; i++){
        target[bit_reverse(long(i), d_block_power)] = input[i];
      }
    }

    void
    polar_encoder::encode_vector(char* target)
    {
      for(int stage = 0; stage < d_block_power; stage++){
        int n_branches = pow(2, stage);
        int branch_elements = d_block_size / (2 * n_branches);
        for(int branch = 0; branch < n_branches; branch++){
          for(int e = 0; e < branch_elements; e++){
            int pos = branch * branch_elements * 2 + e;
//            std::cout << "branch_elements" << branch_elements << ", stage=" << stage << ", branch=" << branch << ", pos=" << pos << std::endl;
            target[pos] ^= target[pos + branch_elements];
          }

        }

      }
    }

    long
    polar_encoder::bit_reverse(long value, int active_bits)
    {
      long r = 0;
      for(int i = 0; i < active_bits; i++){
        r <<= 1;
        r |= value & 1;
        value >>= 1;
      }
      return r;
    }

    void
    polar_encoder::insert_frozen_bits(char* target, const char* input)
    {
      int frozen_num = 0;
      int num_frozen_bits = d_block_size - d_num_info_bits;
      int info_num = 0;
      for(int i = 0; i < d_block_size; i++){
        if(frozen_num < num_frozen_bits && d_frozen_bit_positions.at(frozen_num) == i){
          target[i] = d_frozen_bit_values.at(frozen_num);
          frozen_num++;
        }
        else{
          target[i] = input[info_num];
          info_num++;
        }
      }
    }

  } /* namespace fec */
} /* namespace gr */


