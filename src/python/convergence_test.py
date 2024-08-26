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

import pathlib

import numpy as np

# This function "geweke" was copied from PyMC v3.7.
# That software was distributed under the Apache License, Version 2.0.
# The copyrights for that release are as follows:
# Copyright (c) 2006 Christopher J. Fonnesbeck (Academic Free License)
# Copyright (c) 2007-2008 Christopher J. Fonnesbeck, Anand Prabhakar Patil, David Huard (Academic Free License)
# Copyright (c) 2009-2017 The PyMC developers (see contributors to pymc-devs on GitHub)
def geweke(x, first=.1, last=.5, intervals=20):
    R"""Return z-scores for convergence diagnostics.

    Compare the mean of the first % of series with the mean of the last % of
    series. x is divided into a number of segments for which this difference is
    computed. If the series is converged, this score should oscillate between
    -1 and 1.

    Parameters
    ----------
    x : array-like
      The trace of some stochastic parameter.
    first : float
      The fraction of series at the beginning of the trace.
    last : float
      The fraction of series at the end to be compared with the section
      at the beginning.
    intervals : int
      The number of segments.

    Returns
    -------
    scores : list [[]]
      Return a list of [i, score], where i is the starting index for each
      interval and score the Geweke score on the interval.

    Notes
    -----

    The Geweke score on some series x is computed by:

      .. math:: \frac{E[x_s] - E[x_e]}{\sqrt{V[x_s] + V[x_e]}}

    where :math:`E` stands for the mean, :math:`V` the variance,
    :math:`x_s` a section at the start of the series and
    :math:`x_e` a section at the end of the series.

    References
    ----------
    Geweke (1992)
    """

    if np.ndim(x) > 1:
        return [geweke(y, first, last, intervals) for y in np.transpose(x)]

    # Filter out invalid intervals
    for interval in (first, last):
        if interval <= 0 or interval >= 1:
            raise ValueError(
                "Invalid intervals for Geweke convergence analysis",
                (first,
                 last))
    if first + last >= 1:
        raise ValueError(
            "Invalid intervals for Geweke convergence analysis",
            (first,
             last))

    # Initialize list of z-scores
    zscores = []

    # Last index value
    end = len(x) - 1

    # Start intervals going up to the <last>% of the chain
    last_start_idx = (1 - last) * end

    # Calculate starting indices
    start_indices = np.arange(0, int(last_start_idx), step=int(
        (last_start_idx) / (intervals - 1)))

    # Loop over start indices
    for start in start_indices:
        # Calculate slices
        first_slice = x[start: start + int(first * (end - start))]
        last_slice = x[int(end - last * (end - start)):]

        z = first_slice.mean() - last_slice.mean()
        z /= np.sqrt(first_slice.var() + last_slice.var())

        zscores.append([start, z])

    if intervals is None:
        return np.array(zscores[0])
    else:
        return np.array(zscores)

# my functions

def check_results_list(results_list, results_string):
    for idx, result in enumerate(results_list):
        zscores = None
        zscores = result[:, 1]
        zscore = zscores[0]
        if (zscore < -2.0) or (zscore > 2.0):
            str_to_add = f"failed: {idx}. zscore: "
            results_string += str_to_add
        results_string += str(zscore)
        results_string += "\n"
    return results_string

def check_results_dict(results_dict, results_string):
    for key, result in results_dict.items():
        zscores = None
        zscores = result[:, 1]
        zscore = zscores[0]
        if (zscore < -2.0) or (zscore > 2.0):
            str_to_add = f"failed: {key}. zscore: "
            results_string += str_to_add
        results_string += str(zscore)
        results_string += "\n"
    return results_string

def check_beta(replic_path, results_string, burnin):
    draws_beta = _core.load_arma_cube_np(
        str(replic_path.joinpath("draws_beta.txt")))
    draws_delta = _core.load_arma_ucube_np(
        str(replic_path.joinpath("draws_delta.txt")))
    n_rows = draws_beta.shape[0]
    n_cols = draws_beta.shape[1]
    results_dict = dict()
    for x in range(0, n_rows):
        for y in range(0, n_cols):
            my_key = (x, y)
            results_dict[my_key] = geweke(
                draws_beta[x, y, burnin:], 0.1, 0.5, 2)
    results_string += "checking beta results:\n"
    results_string = check_results_dict(results_dict, results_string)
    results_string += "\n"
    return results_string

def check_gammas(replic_path, results_string, data, burnin):
    K = data["K"]
    results_list = list()
    for k in range(1, K + 1):
        if data["L_k_s"][k - 1] > 2:
            fname = f"draws_gamma_{k}.txt"
            gamma_fullmat = _core.load_arma_mat_np(
                str(replic_path.joinpath(fname)))
            gamma_submat = np.delete(gamma_fullmat, [0, 1, -1], 1)
            n_cols = gamma_submat.shape[1]
            for i in range(n_cols):
                results_list.append(geweke(gamma_submat[burnin:, i],
                                           0.1, 0.5, 2))
    results_string += "checking gamma results:\n"
    results_string = check_results_list(results_list, results_string)
    results_string += "\n"
    return results_string

def check_omega(replic_path, results_string, burnin):
    fname = "draws_omega.txt"
    draws_omega = _core.load_arma_mat_np(
        str(replic_path.joinpath(fname))).ravel()
    draws_omega = draws_omega[burnin:]
    results_list = list()
    results_list.append(geweke(draws_omega, 0.1, 0.5, 2))
    results_string += "checking omega results:\n"
    results_string = check_results_list(results_list, results_string)
    results_string += "\n"
    return results_string

def check_simple_matrix_param(replic_path, results_string, param_name,
                              burnin):
    fname = f"draws_{param_name}.txt"
    draws_param = _core.load_arma_cube_np(str(replic_path.joinpath(fname)))
    n_rows = draws_param[:, :, 0].shape[0]
    n_cols = draws_param[:, :, 0].shape[1]
    results_dict = dict()
    for x in range(0, n_rows):
        for y in range(0, n_cols):
            my_key = (x, y)
            results_dict[my_key] = geweke(draws_param[x, y, burnin:],
                                          0.1, 0.5, 2)
    results_string += f"checking {param_name} results:\n"
    results_string = check_results_dict(results_dict, results_string)
    results_string += "\n"
    return results_string
    
def run_geweke_test(replic_path, data):
    burnin = data["burnin"]
    results_string = ""
    results_string = check_beta(replic_path, results_string,
                                burnin)
    results_string = check_gammas(replic_path, results_string,
                                  data, burnin)
    results_string = check_omega(replic_path, results_string, burnin)
    results_string = check_simple_matrix_param(replic_path,
                                               results_string,
                                               "Rmat",
                                               burnin)
    results_string = check_simple_matrix_param(replic_path,
                                               results_string,
                                               "lambda",
                                               burnin)
    fpath = replic_path.joinpath("convergence_test_results.txt")
    with open(fpath, 'w') as f:
        f.writelines(results_string)
