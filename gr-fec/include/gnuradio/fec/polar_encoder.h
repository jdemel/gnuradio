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


#ifndef INCLUDED_FEC_POLAR_ENCODER_H
#define INCLUDED_FEC_POLAR_ENCODER_H

#include <gnuradio/fec/api.h>
#include <gnuradio/fec/generic_encoder.h>

namespace gr {
  namespace fec {

    /*!
     * \brief POLAR encoder
     *
     */
    class FEC_API polar_encoder : public generic_encoder
    {
    public:
      static generic_encoder::sptr make(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values);

      ~polar_encoder();



      // FECAPI
      void generic_work(void *in_buffer, void *out_buffer);
      double rate(){return (1.0 * get_input_size() / get_output_size());};
      int get_input_size(){return d_num_info_bits;};
      int get_output_size(){return d_block_size;};
      bool set_frame_size(unsigned int frame_size){return false;};

    private:
      polar_encoder(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values);
      int d_block_size; // depending on paper called 'N' or 'm'
      int d_block_power;
      int d_num_info_bits; // mostly abbreviated by 'K'
      std::vector<int> d_frozen_bit_positions;
      std::vector<char> d_frozen_bit_values;

      void insert_frozen_bits(char* target, const char* input);
      void bit_reverse_vector(char* target, const char* input);
      void encode_vector(char* target);

      // helper functions
      long bit_reverse(long value, int active_bits);

    };

  } // namespace fec
} // namespace gr

#endif /* INCLUDED_FEC_POLAR_ENCODER_H */

