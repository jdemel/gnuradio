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

#include <gnuradio/blocks/pack_k_bits.h>
#include <gnuradio/blocks/unpack_k_bits.h>

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

      int k = 8;
      d_unpacker = new gr::blocks::kernel::unpack_k_bits(k);
      d_packer = new gr::blocks::kernel::pack_k_bits(k);
    }

    polar_encoder::~polar_encoder()
    {
      delete d_unpacker;
      delete d_packer;
    }

    void
    polar_encoder::generic_work(void* in_buffer, void* out_buffer)
    {
      const char *in = (const char *) in_buffer;
      char *out = (char *) out_buffer;
      char* uncoded_arr = (char*) volk_malloc(sizeof(char) * d_block_size, volk_get_alignment());
      insert_frozen_bits(uncoded_arr, in);
      bit_reverse_vector(out, uncoded_arr);
//      encode_vector(out);

      unsigned char* temp = (unsigned char*) volk_malloc(d_block_size, volk_get_alignment());
      int num_bytes = d_block_size >> 3;
      d_packer->pack(temp, (unsigned char*) out, num_bytes);
      encode_vector_packed(temp);
      d_unpacker->unpack((unsigned char*) out, temp, num_bytes);
      volk_free(temp);


      volk_free(uncoded_arr);
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
            target[pos] ^= target[pos + branch_elements];
          }
        }
      }
    }

    void
    polar_encoder::encode_vector_packed(unsigned char* target)
    {
      encode_vector_packed_subbyte(target);
      encode_vector_packed_interbyte(target);
    }

    void
    polar_encoder::encode_vector_packed_subbyte(unsigned char* target)
    {
      int shift_vals[3] = {1, 2, 4};
      int bit_masks[3] = {0xaa, 0xcc, 0xF0};
      int num_bytes_per_block = d_block_size >> 3;
      for(int stage = 0; stage < 3; stage++){
        unsigned char* it = target;
        for(int byte = 0; byte < num_bytes_per_block; byte++){
          *it ^= bit_masks[stage] & (*it << shift_vals[stage]);
          it++;
        }
      }
    }

    void
    polar_encoder::encode_vector_packed_interbyte(unsigned char* target)
    {
      int branch_byte_size = 1;
      unsigned char* pos;
      int n_branches = d_block_size >> 4;
      for(int stage = 3; stage < d_block_power; stage++){
        pos = target;
        for(int branch = 0; branch < n_branches; branch++){
          for(int block = 0; block < branch_byte_size; block++){
            *(pos + block) ^= *(pos + block + branch_byte_size);
          }
          pos += branch_byte_size << 1;
        }
        n_branches >>= 1;
        branch_byte_size <<= 1;
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


