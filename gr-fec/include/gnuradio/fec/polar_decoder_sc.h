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

// Forward declaration for those objects. SWIG doesn't like them to be #include'd.
namespace gr {
  namespace blocks {
    namespace kernel {
      class pack_k_bits;
      class unpack_k_bits;
    }
  }
}

namespace gr {
  namespace fec {

    /*!
     * \brief Standard successive cancellation decoder for POLAR codes
     *
     */
    class FEC_API polar_decoder_sc : public generic_decoder
    {
    public:
      static generic_decoder::sptr make(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed = false);
      ~polar_decoder_sc();

      // FECAPI
      void generic_work(void *in_buffer, void *out_buffer);
      double rate(){return (1.0 * get_input_size() / get_output_size());};
      int get_input_size(){return d_num_info_bits / (d_is_packed ? 8 : 1);};
      int get_output_size(){return d_block_size / (d_is_packed ? 8 : 1);};
      bool set_frame_size(unsigned int frame_size){return false;};

    private:
      polar_decoder_sc(int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed);
      int d_block_size; // depending on paper called 'N' or 'm'
      int d_block_power;
      bool d_is_packed;
      int d_num_info_bits; // mostly abbreviated by 'K'
      std::vector<int> d_frozen_bit_positions;
      std::vector<int> d_info_bit_positions;
      std::vector<char> d_frozen_bit_values;
    };

  } // namespace fec
} // namespace gr

#endif /* INCLUDED_FEC_POLAR_DECODER_SC_H */

