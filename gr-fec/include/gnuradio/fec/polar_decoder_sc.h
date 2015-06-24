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


#ifndef INCLUDED_FEC_POLAR_DECODER_SC_H
#define INCLUDED_FEC_POLAR_DECODER_SC_H

#include <gnuradio/fec/api.h>
#include <gnuradio/fec/generic_decoder.h>
#include <gnuradio/fec/polar_common.h>


namespace gr {
  namespace fec {

    /*!
     * \brief Standard successive cancellation decoder for POLAR codes
     * It expects float input with bits mapped 1 --> -1, 0 --> 1
     * Or: f = 1.0 - 2.0 * bit
     *
     */
    class FEC_API polar_decoder_sc : public generic_decoder, public polar_common
    {
    public:
      static generic_decoder::sptr make(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed = false);
      ~polar_decoder_sc();

      // FECAPI
      void generic_work(void *in_buffer, void *out_buffer);
      double rate(){return (1.0 * get_input_size() / get_output_size());};
      int get_input_size(){return block_size() / (is_packed() ? 8 : 1);};
      int get_output_size(){return num_info_bits() / (is_packed() ? 8 : 1);};
      bool set_frame_size(unsigned int frame_size){return false;};

    private:
      polar_decoder_sc(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed);
      const float D_LLR_FACTOR;
      std::vector<int> d_frozen_bit_positions;
      std::vector<int> d_info_bit_positions;
      std::vector<char> d_frozen_bit_values;

      unsigned int d_frozen_bit_counter;

      float* d_llr_vec;
      void print_pretty_llr_vector() const;
      void initialize_llr_vector(float* llrs, const float* input);
      unsigned char* d_u_hat_vec;

      float llr_odd(const float la, const float lb) const;
      float llr_even(float la, float lb, unsigned char f) const;
      unsigned char llr_bit_decision(float llr) const {return (llr < 0.0f) ? 1 : 0;};
      unsigned char retrieve_bit_from_llr(float llr, const int pos);
      void sc_decode(float* llrs, unsigned char* u);
      void butterfly(float* llrs, int call_row, int stage, unsigned char* u, const int u_num);
      void even_u_values(unsigned char* u_even, const unsigned char* u, const int u_num);
      void odd_xor_even_values(unsigned char* u_xor, const unsigned char* u, const int u_num);

      void extract_info_bits(unsigned char* output, const unsigned char* input);



    };

  } // namespace fec
} // namespace gr

#endif /* INCLUDED_FEC_POLAR_DECODER_SC_H */

