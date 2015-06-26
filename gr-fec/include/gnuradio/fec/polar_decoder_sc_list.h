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


#ifndef INCLUDED_POLAR_FEC_DECODER_SC_LIST_H
#define INCLUDED_POLAR_FEC_DECODER_SC_LIST_H

#include <gnuradio/fec/api.h>
#include <gnuradio/fec/polar_decoder_common.h>

namespace gr {
  namespace fec {

    /*!
     * \brief implements a successive cancellation list decoder for polar codes
     *
     */
    class FEC_API polar_decoder_sc_list : public polar_decoder_common
    {
    public:
      static generic_decoder::sptr make(int max_list_size, int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed = false);
      ~polar_decoder_sc_list();

      // FECAPI
      void generic_work(void *in_buffer, void *out_buffer);

    private:
      polar_decoder_sc_list(int max_list_size, int block_size, int num_info_bits, std::vector<int> frozen_bit_positions, std::vector<char> frozen_bit_values, bool is_packed);
      int d_max_list_size;

      // just a class to hold all the necessary info for the list of paths
      class path
      {
      public:
        path(int block_size, int block_power);
        ~path();
        void update_metrics(const int u_num, const int pos);
        void set_ui(const unsigned char ui, const int u_num){u_hat_vec[u_num] = ui; path_metric = ui ? path_metric1 : path_metric0;};
        float calculate_path_metric(const float last_pm, const float llr, const unsigned char u_hat) const;
        void duplicate_path(const path* original_path, const int block_size, const int block_power);
        float* llr_vec;
        unsigned char* u_hat_vec;
        float path_metric;
        float path_metric0;
        float path_metric1;
        bool is_active;
      };
      typedef boost::shared_ptr<path> path_sptr;
      std::vector<path_sptr> d_path_list;
      unsigned int d_frozen_bit_counter;
      unsigned int d_active_path_counter;
      void activate_path(int old_path_num, int new_path_num);
      void kill_path(int num);



      void decode_list();
      void update_active_paths(int u_num);
      void calculate_next_llr_in_paths(int u_num);
      void calculate_next_llr(path_sptr current_path, int u_num);




    };

  } // namespace fec
} // namespace gr

#endif /* INCLUDED_POLAR_FEC_DECODER_SC_LIST_H */

