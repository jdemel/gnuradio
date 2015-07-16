#!/usr/bin/env python
#
# Copyright 2015 Free Software Foundation, Inc.
#
# This file is part of GNU Radio
#
# GNU Radio is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# GNU Radio is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Radio; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.

# turn this folder into a Python module

import channel_construction as cc
from channel_construction_bec import bhattacharyya_bounds


def load_frozen_bits_info(is_prototype, channel, block_size, num_info_bits, design_snr, mu):
    num_frozen_bits = block_size - num_info_bits
    z_params = bhattacharyya_bounds(design_snr, block_size)
    data_set = {
        'positions': cc.get_frozen_bit_indices_from_z_parameters(z_params, num_frozen_bits),
        'values': [0, ] * num_frozen_bits,
        'block_size': block_size,
        'num_info_bits': num_info_bits,
        'num_frozenbits': num_frozen_bits,
        'design_snr': design_snr,
        'channel': channel,
        'mu': mu,
        }

    if is_prototype:
        return data_set

    if channel == 'BSC':
        z_params = cc.load_z_parameters(block_size, design_snr, mu)
        data_set['positions'] = cc.get_frozen_bit_indices_from_z_parameters(z_params, num_frozen_bits)

    return data_set