# This file is part of "probitlcm" which is released under GPL v3.
#
# Copyright (c) 2022-2024 Eric Alan Wayman <ewayman2@illinois.edu>.
#
# This program is FLO (free/libre/open) software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import pandas as pd
import numpy as np
from probitlcm import _core

def load_single_value(dir_path, fname):
    # tmpmat = np.empty((1, 1))
    # tmpmat[0, 0] = ear
    fpath = dir_path.joinpath(fname)
    ear_mat = _core.load_arma_mat_np(str(fpath))
    ear = ear_mat.item()
    return ear

def save_single_value(ear, dir_path, fname):
    tmpmat = np.empty((1, 1))
    tmpmat[0, 0] = ear
    fpath = dir_path.joinpath(fname)
    _core.save_arma_mat_np(tmpmat, str(fpath))

def convert_table_to_html(table, modify_decimal):
    if modify_decimal is False:
        return pd.DataFrame(table).to_html()
    else:
        return pd.DataFrame(table).to_html(float_format='{:10.2f}'.format)
