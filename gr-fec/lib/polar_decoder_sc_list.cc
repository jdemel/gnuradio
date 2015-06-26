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
#include <gnuradio/fec/polar_decoder_sc_list.h>
#include <volk/volk.h>

#include <cmath>

namespace gr
{
  namespace fec
  {

    generic_decoder::sptr
    polar_decoder_sc_list::make(int max_list_size, int block_size, int num_info_bits,
                                std::vector<int> frozen_bit_positions,
                                std::vector<char> frozen_bit_values, bool is_packed)
    {
      return generic_decoder::sptr(
          new polar_decoder_sc_list(max_list_size, block_size, num_info_bits, frozen_bit_positions,
                                    frozen_bit_values, is_packed));
    }

    polar_decoder_sc_list::polar_decoder_sc_list(int max_list_size, int block_size,
                                                 int num_info_bits,
                                                 std::vector<int> frozen_bit_positions,
                                                 std::vector<char> frozen_bit_values,
                                                 bool is_packed) :
            polar_decoder_common(block_size, num_info_bits, frozen_bit_positions, frozen_bit_values,
                                 is_packed), d_max_list_size(max_list_size), d_frozen_bit_counter(0)
    {
      for(int i = 0; i < max_list_size; i++){
        d_path_list.push_back(path_sptr(new path(block_size, block_power())));
      }
    }

    polar_decoder_sc_list::~polar_decoder_sc_list()
    {
      d_path_list.clear();
    }

    void
    polar_decoder_sc_list::generic_work(void* in_buffer, void* out_buffer)
    {
      const float *in = (const float*) in_buffer;
      unsigned char *out = (unsigned char*) out_buffer;
      std::cout << "generic_work WORK\n";
      initialize_llr_vector(d_path_list[0]->llr_vec, in);
      d_path_list[0]->is_active = true;
      decode_list();
      print_pretty_llr_vector(d_path_list[0]->llr_vec);
      extract_info_bits(out, d_path_list[0]->u_hat_vec);
      std::cout << "generic_work finished\n";
    }

    void
    polar_decoder_sc_list::decode_list()
    {

      for(int i = 0; i < block_size(); i++){
        calculate_next_llr_in_paths(i);
        update_active_paths(i);
      }
    }

    void
    polar_decoder_sc_list::calculate_next_llr_in_paths(int u_num)
    {
      for(unsigned int i = 0; i < d_path_list.size(); i++){
        if(d_path_list[i]->is_active){
          calculate_next_llr(d_path_list[i], u_num);
        }
      }
    }

    void
    polar_decoder_sc_list::update_active_paths(int u_num)
    {
      const int row = bit_reverse(u_num, block_power());


      // 1. if frozen bit, update with known value
      if(d_frozen_bit_counter < d_frozen_bit_positions.size() && u_num == d_frozen_bit_positions.at(d_frozen_bit_counter)){
        unsigned char frozen_bit = d_frozen_bit_values[d_frozen_bit_counter];
        for(unsigned int i = 0; i < d_path_list.size(); i++){
          if(d_path_list[i]->is_active){
            d_path_list[i]->set_ui(frozen_bit, u_num);
          }
        }
      }

      // 2. info bit, not all paths in use
//      else if(d_active_path_counter < d_max_list_size){
//        for(int i = 0; i < d_path_list.size(); i++){
//          if(d_path_list[i]->is_active){
//            for(int e = 0; e < d_path_list.size(); e++){
//              if(!d_path_list[e]->is_active){
//                activate_path(i, e);
//                d_path_list[i]->set_ui(0, u_num);
//                d_path_list[e]->set_ui(1, u_num);
//              }
//            }
//          }
//        }
//      }
      d_path_list[0]->set_ui(llr_bit_decision(d_path_list[0]->llr_vec[row]), u_num);
//      float metric_sum = 0;
//      int active_path_counter = 0;
//      for(int i = 0; i < d_path_list.size(); i++){
//        metric_sum += d_path_list[i]->path_metric0;
//        metric_sum += d_path_list[i]->path_metric1;
//        active_path_counter++;
//      }
//      metric_sum /= 2 * active_path_counter;

    }

    void
    polar_decoder_sc_list::activate_path(int old_path_num, int new_path_num)
    {
      d_path_list[new_path_num]->duplicate_path(d_path_list[old_path_num].get(), block_size(), block_power());
      d_active_path_counter++;
    }

    void
    polar_decoder_sc_list::kill_path(int num)
    {
      d_path_list[num]->is_active = false;
      d_active_path_counter--;
    }

    void
    polar_decoder_sc_list::calculate_next_llr(path_sptr current_path, int u_num)
    {
      int row = bit_reverse(u_num, block_power());
      butterfly(current_path->llr_vec, row, 0, current_path->u_hat_vec, u_num);
      current_path->update_metrics(u_num, row);
    }



    polar_decoder_sc_list::path::path(int block_size, int block_power):
        path_metric(0.0f),
        path_metric0(0.0f),
        path_metric1(0.0f),
        is_active(false)
    {
      llr_vec = (float*) volk_malloc(sizeof(float) * block_size * (block_power + 1), volk_get_alignment());
      u_hat_vec = (unsigned char*) volk_malloc(sizeof(unsigned char) * block_size * (block_power + 1), volk_get_alignment());
    }

    polar_decoder_sc_list::path::~path()
    {
      volk_free(llr_vec);
      volk_free(u_hat_vec);
    }

    void
    polar_decoder_sc_list::path::update_metrics(const int u_num, const int pos)
    {
      path_metric0 = calculate_path_metric(path_metric, *(llr_vec + pos), 0);
      path_metric1 = calculate_path_metric(path_metric, *(llr_vec + pos), 1);
    }

    float
    polar_decoder_sc_list::path::calculate_path_metric(const float last_pm, const float llr,
                                                 const unsigned char u_hat) const
    {
//      if((u_hat == 0 && llr > 0.0f) || (u_hat == 1 && llr < 0.0f)){
      if(u_hat == (unsigned char) (0.5 * 1 - copysignf(1.0f, llr))){
        return last_pm;
      }
      return last_pm + abs(llr);
    }

    void
    polar_decoder_sc_list::path::duplicate_path(const path* original_path, const int block_size, const int block_power)
    {
      memcpy(llr_vec, original_path->llr_vec, sizeof(float) * block_size * (block_power + 1));
      memcpy(u_hat_vec, original_path->u_hat_vec, sizeof(unsigned char) * block_size * (block_power + 1));
      path_metric = original_path->path_metric;
      is_active = true;
    }
  } /* namespace fec */
} /* namespace gr */


