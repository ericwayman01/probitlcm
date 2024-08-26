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

from probitlcm import _core

import argparse, json, pathlib

import sys
if sys.version_info[1] < 11:
    import toml
else:
    import tomllib as toml

import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib.pyplot as plt


def run_stuff(dataset_dir, data_analysis_path, N_files, chainnum):
    # set up required quantities
    chainnum_dirname = f"chain_{chainnum:03}"    
    chain_results_path = data_analysis_path.joinpath(chainnum_dirname)
    # load Y data
    fpath = dataset_dir.joinpath("responses.txt")
    Y_obs = _core.load_arma_umat_np(str(fpath))
    Y_obs_int = Y_obs.astype(int)
    # load M_j_s
    jsonfilepath = dataset_dir.joinpath("M_j_s.json")
    M_j_s = json.loads(jsonfilepath.read_bytes())
    M_j_s = np.asarray(M_j_s, dtype=np.uint)
    # do work
    N = Y_obs.shape[0]
    J = Y_obs.shape[1]
    ## dimensions
    # N = 3960
    # J = 17
    # N_files = 6
    ## program params
    size_of_chunk = 100
    size_of_chunk_zb = size_of_chunk - 1
    S_thinned = N_files * size_of_chunk
    ## data structure
    Ymat_pred = np.empty((N, J, S_thinned), dtype=np.uint)
    # load Y pred draws
    for i in range(1, N_files + 1):
        min_idx_val = size_of_chunk * (i - 1)
        max_idx_val = size_of_chunk * i
        fpath = chain_results_path.joinpath(f"Ymat_pred_chunk_{i:03}.txt")
        Ymat_pred[:, :, min_idx_val:max_idx_val] = _core.load_arma_ucube_np(
            str(fpath))
    Ymat_pred_int = Ymat_pred.astype(int)
    
    counts_Ymat = _core.aggregateAllCounts_np(Y_obs, M_j_s)
    counts_Ymat = counts_Ymat.ravel().astype(int)
    counts_Ymat_pred = np.empty((S_thinned, len(counts_Ymat)), dtype=int)
    tmparray = np.empty((N, J), dtype=np.uint)
    for s in range(1, S_thinned + 1):
        tmparray = np.copy(Ymat_pred[:, :, s - 1])
        counts_Ymat_pred[s - 1] = _core.aggregateAllCounts_np(
            tmparray, M_j_s).ravel().astype(int)
    ## distance: one vector minus the other, absolute value, sum
    mann_whitney_xvar = np.empty((S_thinned))
    for s in range(1, S_thinned + 1):
        mann_whitney_xvar[s - 1] = np.sum(np.abs(counts_Ymat_pred[s - 1] -
                                                 counts_Ymat))
    mann_whitney_yvar = np.empty(int(S_thinned * (S_thinned - 1) / 2))
    idx = 1
    for s in range(1, S_thinned + 1):
        for t in range(s + 1, S_thinned + 1):
            mann_whitney_yvar[idx - 1] = np.sum(
                np.abs(counts_Ymat_pred[s - 1] - counts_Ymat_pred[t - 1]))
            idx += 1
    # calc p_value
    statistic, p_value = mannwhitneyu(mann_whitney_xvar, mann_whitney_yvar,
                                      alternative="greater")
    stat_mat = np.empty((1))
    stat_mat[0] = statistic
    p_value_mat = np.empty((1))
    p_value_mat[0] = p_value
    _core.save_arma_mat_np(stat_mat, str(
        chain_results_path.joinpath("bartholomew_stat.txt")))
    _core.save_arma_mat_np(p_value_mat, str(
        chain_results_path.joinpath("bartholomew_p_value.txt")))
    bins = np.linspace(min(np.min(mann_whitney_xvar),
                           np.min(mann_whitney_yvar)),
                       max(np.max(mann_whitney_xvar),
                           np.max(mann_whitney_yvar)), 100)
    plt.hist(mann_whitney_xvar, bins, alpha=0.5, label='x', density=True,
             stacked=True)
    plt.hist(mann_whitney_yvar, bins, alpha=0.5, label='y', density=True,
             stacked=True)
    plt.legend(loc='upper right')
    plt.savefig(chain_results_path.joinpath("bartholomew_histogram.png"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", required=True)
    parser.add_argument("--setup_num", required=True)
    parser.add_argument("--chainnum", required=True)
    parser.add_argument("--dataset", required=True) # e.g. stepbd
    parser.add_argument("--N_files", required=True)
    args = parser.parse_args()
    N_files = int(args.N_files)
    chainnum = int(args.chainnum)
    # same format for config file    
    run_dir = pathlib.Path.cwd()
    # load config file
    with open('config_data_analysis.toml') as fileObj:
        config = toml.load(fileObj)
    if args.environ == "laptop":
        process_dir = config['laptop_process_dir']
    elif args.environ == "cluster":
        process_dir = config['cluster_process_dir']
    process_dir = pathlib.Path(process_dir)
    # set up path to dataset_dir
    dataset = args.dataset
    setup_num = int(args.setup_num)
    jsonfilename_stem = f"setup_{setup_num:04}"
    # set up path to dataset_dir
    dataset_dir = run_dir.joinpath(dataset)
    # set up data_analysis_path (like the process_dir of scenario_launch)
    data_analysis_path = process_dir.joinpath(dataset, jsonfilename_stem)
    run_stuff(dataset_dir, data_analysis_path, N_files, chainnum)
