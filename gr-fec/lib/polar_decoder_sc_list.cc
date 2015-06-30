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
#include <algorithm>

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
                                 is_packed), d_max_list_size(max_list_size), d_frozen_bit_counter(0),
                                 d_active_path_counter(0)
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
//      std::cout << "generic_work WORK\n";
      initialize_llr_vector(d_path_list[0]->llr_vec, in);
      activate_path(0, 0);
      decode_list();
//      for(unsigned int i = 0; i < d_path_list.size(); i++){
//        std::cout << "\npath num " << i << std::endl;
//        print_pretty_llr_vector(d_path_list[i]->llr_vec);
//      }
//      print_pretty_llr_vector(d_path_list[0]->llr_vec);
      extract_info_bits(out, d_path_list[find_survivor()]->u_hat_vec);
//      std::cout << "generic_work finished\n";
    }

    int
    polar_decoder_sc_list::find_survivor() const
    {
      float path_val = 200000000.0f;
      int survivor = 0;
      for(unsigned int i =0 ; i < d_path_list.size(); i++){
        if(d_path_list[i]->is_active && d_path_list[i]->path_metric < path_val){
          survivor = i;
          path_val = d_path_list[i]->path_metric;
        }
      }
      std::cout << "path survivor: " << survivor << std::endl;
      return survivor;
    }

    void
    polar_decoder_sc_list::decode_list()
    {

      for(int i = 0; i < block_size(); i++){
        std::cout << "\n\ndecode n = " << i << std::endl;
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
//          std::cout << "calculate_next_llr_paths, npath = " << i << std::endl;
//          print_pretty_llr_vector(d_path_list[i]->llr_vec);
        }
      }
    }

    void
    polar_decoder_sc_list::update_active_paths(const int u_num)
    {
//      const int row = bit_reverse(u_num, block_power());
      std::cout << "update_active_paths, u_num = " << u_num << ", frozenbit_num = " << d_frozen_bit_counter << std::endl;

      // 1. if frozen bit, update with known value
      if(d_frozen_bit_counter < d_frozen_bit_positions.size() && u_num == d_frozen_bit_positions.at(d_frozen_bit_counter)){
        update_with_frozenbit(u_num);
      }

      // 2. info bit, not all paths in use
      else if(d_active_path_counter < d_max_list_size){
        add_active_paths(u_num);
      }

      // 3. info bit, maximum active paths reached.
      else{
        select_best_paths(u_num);
      }
    }


    void
    polar_decoder_sc_list::select_best_paths(const int u_num)
    {
      std::vector<float> metrics;
      metrics.reserve(2 * block_size());
      for(unsigned int i = 0; i < d_path_list.size(); i++){
        metrics.push_back(d_path_list[i]->path_metric0);
        metrics.push_back(d_path_list[i]->path_metric1);
      }
      std::sort(metrics.begin(), metrics.end());
      const float median = (metrics[d_max_list_size - 1] + metrics[d_max_list_size]) / 2;
      for(unsigned int i = 0; i < metrics.size(); i++){
        std::cout << metrics[i] << ", ";
      }

      std::cout << "\nselect_best_paths, u_num = " << u_num << ", median = " << median << std::endl;
//      unsigned int kill_paths = 0;
//      unsigned int equal_paths = 0;
//      unsigned int duplicate_paths = 0;
//
//      for(unsigned int i = 0; i < d_path_list.size(); i++){
//        if(d_path_list[i]->path_metric0 > median && d_path_list[i]->path_metric1 > median){
//          kill_paths++;
//        }
//        else if(d_path_list[i]->path_metric0 == median && d_path_list[i]->path_metric1 == median){
//          equal_paths++;
//        }
//        else if(d_path_list[i]->path_metric0 < median && d_path_list[i]->path_metric1 < median){
//          duplicate_paths++;
//        }
//      }

      for(unsigned int i = 0; i < d_path_list.size(); i++){
        if(d_path_list[i]->path_metric0 >= median && d_path_list[i]->path_metric1 >= median){
          kill_path(i);
        }
      }

//      unsigned int n_equal_kill = std::min(0, d_active_path_counter - (d_max_list_size / 2));
//      for(unsigned int i = 0; i < d_path_list.size(); i++){
//        if(d_path_list[i]->path_metric0 == median && d_path_list[i]->path_metric1 == median){
//          kill_path(i);
//        }
//      }

      for(unsigned int i = 0; i < d_path_list.size(); i++){
        if(d_path_list[i]->path_metric0 < median && d_path_list[i]->path_metric1 < median){
          duplicate_and_set_path(i, u_num);
        }
        else{
          const unsigned char u_hat = (d_path_list[i]->path_metric0 < d_path_list[i]->path_metric1) ? 0 : 1;
          d_path_list[i]->set_ui(u_hat, u_num);
        }
      }
    }

    void
    polar_decoder_sc_list::add_active_paths(const int u_num)
    {
      // make sure you know all the formerly active paths beforehand.
      std::vector<int> active_paths;
      for(unsigned int i = 0; i < d_path_list.size(); i++){
        if(d_path_list[i]->is_active){
          active_paths.push_back(i);
        }
      }

      // duplicate and set new
      for(unsigned int i = 0; i < active_paths.size(); i++){
        const int o_active = active_paths[i];
        duplicate_and_set_path(o_active, u_num);
      }
    }

    void
    polar_decoder_sc_list::duplicate_and_set_path(const int nactive, const int u_num)
    {
      if(! (d_active_path_counter < d_max_list_size - 1)){
        const unsigned char u_hat = (d_path_list[nactive]->path_metric0 < d_path_list[nactive]->path_metric1) ? 0 : 1;
        d_path_list[nactive]->set_ui(u_hat, u_num);
      }
      int n_active_pos = 0;
      while(d_path_list[n_active_pos]->is_active){
        n_active_pos++;
      }
      activate_path(nactive, n_active_pos);
      d_path_list[nactive]->set_ui(0, u_num);
      d_path_list[n_active_pos]->set_ui(1, u_num);
    }

    void
    polar_decoder_sc_list::activate_path(int old_path_num, int new_path_num)
    {
      d_path_list[new_path_num]->duplicate_path(d_path_list[old_path_num].get(), block_size(), block_power());
      d_active_path_counter++;
//      print_pretty_llr_vector(d_path_list[new_path_num]->llr_vec);
    }

    void
    polar_decoder_sc_list::kill_path(int num)
    {
      std::cout << "killpath = " << num << std::endl;
      d_path_list[num]->is_active = false;
      d_active_path_counter--;
    }

    void
    polar_decoder_sc_list::update_with_frozenbit(const int u_num)
    {
      unsigned char frozen_bit = d_frozen_bit_values[d_frozen_bit_counter];

//      std::cout << "update_frozen_bit u_num = " << u_num << ", with " << int(frozen_bit) << std::endl;
      for(unsigned int i = 0; i < d_path_list.size(); i++){
        if(d_path_list[i]->is_active){
          d_path_list[i]->set_ui(frozen_bit, u_num);
        }
      }
      d_frozen_bit_counter++;
    }

    void
    polar_decoder_sc_list::calculate_next_llr(path_sptr current_path, int u_num)
    {
      int row = bit_reverse(u_num, block_power());
      std::cout << "\t\t\t\tcalc_next_llr = " << u_num << "->" << row << ", metric = " << current_path->path_metric << std::endl;
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
      memset(llr_vec, 0, sizeof(float) * block_size * (block_power + 1));
      u_hat_vec = (unsigned char*) volk_malloc(sizeof(unsigned char) * block_size * (block_power + 1), volk_get_alignment());
      memset(u_hat_vec, 0, sizeof(unsigned char) * block_size * (block_power + 1));
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
      std::cout << "metrics updated! " << *(llr_vec + pos) << "\t0: " << path_metric0 << "\t1: " << path_metric1 << std::endl;
    }

    float
    polar_decoder_sc_list::path::calculate_path_metric(const float last_pm, const float llr,
                                                 const unsigned char u_hat) const
    {
//      if((u_hat == 0 && llr > 0.0f) || (u_hat == 1 && llr < 0.0f)){
      if(u_hat == (unsigned char) (0.5 * 1 - copysignf(1.0f, llr))){
        return last_pm;
      }
      return last_pm + fabs(llr);
    }

    void
    polar_decoder_sc_list::path::set_ui(const unsigned char ui, const int u_num)
    {
      u_hat_vec[u_num] = ui;
      path_metric = ui ? path_metric1 : path_metric0;
      std::cout << "set_ui( " << int(ui) << ", " << u_num << " ) metrics: 0: " << path_metric0 << ", 1: " << path_metric1 << " --> " << path_metric << std::endl;
    }

    void
    polar_decoder_sc_list::path::duplicate_path(const path* original_path, const int block_size, const int block_power)
    {
      memcpy(llr_vec, original_path->llr_vec, sizeof(float) * block_size * (block_power + 1));
      memcpy(u_hat_vec, original_path->u_hat_vec, sizeof(unsigned char) * block_size * (block_power + 1));
      path_metric = original_path->path_metric;
      path_metric0 = original_path->path_metric0;
      path_metric1 = original_path->path_metric1;
      is_active = true;
    }
  } /* namespace fec */
} /* namespace gr */


