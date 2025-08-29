# This file is part of "probitlcm" which is released under GPL v3.
#
# Copyright (c) 2022-2025 Eric Alan Wayman <ericwaymanpublications@mathworks.org>.
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

import json
import math
import pathlib

import jinja2
import numpy as np
import pandas as pd

from probitlcm import report_helpers
from probitlcm import _core

# Note: In present form, the HTML reports for all values of statistic_type
#     have the same filename.

def calc_mcse_cube(my_cube, number_of_replics):
    return np.std(my_cube, axis=2) / math.sqrt(number_of_replics)

def calculate_stats_for_beta_and_delta(schema_file_path, scenario_path,
                                       number_of_replics, data,
                                       statistic_type):
    # dict to store results
    results_dict = dict()
    # set up data structures from data
    L_k_s = data["L_k_s"]
    H = data["H"]
    per_effect_M_j_s = data["per_effect_M_j_s"]
    q = data["q"]
    J = q * len(per_effect_M_j_s)
    mycube_beta = np.empty((H, J, number_of_replics))
    mycube_delta = np.empty((H, J, number_of_replics))
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        fpath = scenario_path.joinpath(replicnum_dirname, "stat_beta.txt")
        mycube_beta[:, :, replicnum - 1] = _core.load_arma_mat_np(str(fpath))
        fpath = scenario_path.joinpath(replicnum_dirname, "stat_delta.txt")
        mycube_delta[:, :, replicnum - 1] = _core.load_arma_mat_np(str(fpath))
    # calculate beta stats and delta stats
    # mae can mean "median" or "mean" absolute error
    mae_beta = None
    mae_delta = None
    if statistic_type == "mean":
        mae_beta = np.mean(mycube_beta, axis=2) # avg_beta_stat
        mae_delta = np.mean(mycube_delta, axis=2)
    elif statistic_type == "median":
        mae_beta = np.median(mycube_beta, axis=2)
        mae_delta = np.median(mycube_delta, axis=2)
    avg_of_mae_beta = np.mean(mae_beta)
    avg_of_mae_delta = np.mean(mae_delta)
    fname = f"avg_of_mae_beta_{statistic_type}.txt"
    report_helpers.save_single_value(avg_of_mae_beta, scenario_path, fname)
    fname = f"avg_of_mae_delta_{statistic_type}.txt"
    report_helpers.save_single_value(avg_of_mae_delta, scenario_path,
                                     fname)
    mae_beta_html = report_helpers.convert_table_to_html(
        mae_beta, True)
    mae_delta_html = report_helpers.convert_table_to_html(
        mae_delta, True)
    # save monte carlo errors
    avg_mcse_beta = np.mean(calc_mcse_cube(mycube_beta, number_of_replics))
    avg_mcse_delta = np.mean(calc_mcse_cube(mycube_delta, number_of_replics))
    fname = "avg_mcse_beta.txt"
    report_helpers.save_single_value(avg_mcse_beta, scenario_path,
                                     fname)
    fname = "avg_mcse_delta.txt"
    report_helpers.save_single_value(avg_mcse_delta, scenario_path,
                                     fname)
    return mae_beta_html, avg_of_mae_beta, mae_delta_html, avg_of_mae_delta

def calculate_theta_stats(schema_file_path, scenario_path,
                          number_of_replics, data, statistic_type):
    # set up data structures from data
    per_effect_M_j_s = data["per_effect_M_j_s"]
    q = data["q"]
    J = q * len(per_effect_M_j_s)
    H_K = data["H_K"]
    fpath = scenario_path.joinpath("M_j_s.txt")
    M_j_s = _core.load_arma_umat_np(str(fpath))
    M_j_s = np.transpose(M_j_s).ravel()
    theta_j_avg_stats = np.empty((J))
    # do theta cubes one at a time
    for j in range(1, J + 1):
        theta_cube = np.empty((H_K, M_j_s[j - 1], number_of_replics))
        # load data into cube
        for replicnum in range(1, number_of_replics + 1):
            replicnum_dirname = f"replic_{replicnum:03}"
            fpath = scenario_path.joinpath(
                replicnum_dirname, f'stat_theta_j_{j:03}.txt')
            theta_cube[:, :, replicnum - 1] = _core.load_arma_mat_np(
                str(fpath))
        mae_thetas = None
        if statistic_type == "mean":
            mae_thetas = np.mean(theta_cube, axis=2)
        elif statistic_type == "median":
            mae_thetas = np.median(theta_cube, axis=2)
        theta_j_avg_stats[j - 1] = np.mean(mae_thetas)
    avg_of_avg_theta_stat = np.mean(theta_j_avg_stats)
    fname = f"avg_of_avg_theta_stat_{statistic_type}.txt"
    report_helpers.save_single_value(avg_of_avg_theta_stat,
                                     scenario_path, fname)
    avg_theta_stat_html = report_helpers.convert_table_to_html(
        theta_j_avg_stats, True)
    # save monte carlo errors
    avg_mcse_theta = np.mean(calc_mcse_cube(theta_cube, number_of_replics))
    fname = "avg_mcse_theta.txt"
    report_helpers.save_single_value(avg_mcse_theta, scenario_path,
                                     fname)
    return avg_theta_stat_html, avg_of_avg_theta_stat

def calculate_matrix_param_stats(schema_file_path, scenario_path,
                                 number_of_replics, data, statistic_type,
                                 param_name, n_rows):
    K = data["K"]
    mycube_param = np.empty((n_rows, K, number_of_replics))
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        fname = "stat_" + param_name + ".txt"
        fpath = scenario_path.joinpath(replicnum_dirname, fname)
        mycube_param[:, :, replicnum - 1] = _core.load_arma_mat_np(str(fpath))
    mae_param = None
    if statistic_type == "mean":
        mae_param = np.mean(mycube_param, axis=2)
    elif statistic_type == "median":
        mae_param = np.median(mycube_param, axis=2)
    avg_of_mae_param = np.mean(mae_param)
    fname = f"avg_of_mae_{param_name}_{statistic_type}.txt"
    report_helpers.save_single_value(avg_of_mae_param,
                                     scenario_path, fname)
    mae_param_html = report_helpers.convert_table_to_html(
        mae_param, True)
    # save monte carlo errors
    avg_mcse_theta = np.mean(calc_mcse_cube(mycube_param, number_of_replics))
    fname = f"avg_mcse_{param_name}.txt"
    report_helpers.save_single_value(avg_mcse_theta, scenario_path,
                                     fname)
    return mae_param_html, avg_of_mae_param

# note: this is almost the same as the theta logic except for the L_k_s check
def calculate_gamma_stats(schema_file_path, scenario_path,
                          number_of_replics, data, statistic_type):
    # dict to store results
    results_dict = dict()
    # set up data structures from data
    L_k_s = data["L_k_s"]
    K = data["K"]
    gamma_avg_stats = np.empty(K)
    avg_gamma_mcse_values = np.empty(K)
    for k in range(1, K + 1):
        if L_k_s[k - 1] == 2:
            continue
        else:
            L_k = L_k_s[k - 1]
            # note that (L_k + 1) - 3 = L_k - 2
            gamma_cube = np.empty((1, L_k - 2, number_of_replics))
            for replicnum in range(1, number_of_replics + 1):
                replicnum_dirname = f"replic_{replicnum:03}"
                fname = f"stat_gamma_{k}.txt"
                fpath = scenario_path.joinpath(replicnum_dirname, fname)
                gamma_cube[:, :, replicnum - 1] = _core.load_arma_mat_np(
                    str(fpath))
            mae_gammas = None
            if statistic_type == "mean":
                mae_gammas = np.mean(gamma_cube, axis=2)
            elif statistic_type == "median":
                mae_gammas = np.median(gamma_cube, axis=2)
            avg_mae_gamma_k = np.mean(mae_gammas)
            gamma_avg_stats[k - 1] = avg_mae_gamma_k
            avg_mcse_gamma_k = np.mean(calc_mcse_cube(
                gamma_cube, number_of_replics))
            avg_gamma_mcse_values[k - 1] = avg_mcse_gamma_k
        fname = f"avg_of_mae_gamma_{k}_{statistic_type}.txt"
        report_helpers.save_single_value(avg_mae_gamma_k,
                                         scenario_path, fname)
    gamma_avg_of_avg_of_mae = np.mean(gamma_avg_stats)
    fname = f"avg_of_avg_of_mae_gamma_{statistic_type}.txt"
    report_helpers.save_single_value(gamma_avg_of_avg_of_mae,
                                     scenario_path,
                                     fname)
    # save monte carlo errors
    fname = f"avg_mcse_gamma.txt"
    avg_mcse_gamma = np.mean(avg_gamma_mcse_values)
    report_helpers.save_single_value(avg_mcse_gamma, scenario_path,
                                     fname)
    return gamma_avg_stats, gamma_avg_of_avg_of_mae

def calc_avg_of_a_metric_value(scenario_path, number_of_replics,
                               fname):
    sum_of_metric_values = 0
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        fpath = scenario_path.joinpath(replicnum_dirname, fname)
        replic_metric_value_mat = _core.load_arma_mat_np(str(fpath))
        replic_metric_value = replic_metric_value_mat.item()
        sum_of_metric_values += replic_metric_value
    avg_of_metric_values = sum_of_metric_values / number_of_replics
    return avg_of_metric_values

def do_beta_delta_subset_calcs(schema_file_path, scenario_path,
                               number_of_replics, data,
                               statistic_type):
    results_dict = dict()
    # set up data structures from data
    L_k_s = data["L_k_s"]
    H = data["H"]
    per_effect_M_j_s = data["per_effect_M_j_s"]
    q = data["q"]
    J = q * len(per_effect_M_j_s)
    mycube_beta = np.empty((H, J, number_of_replics))
    mycube_delta = np.empty((H, J, number_of_replics))
    # load cube
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        fpath = scenario_path.joinpath(replicnum_dirname, "stat_beta.txt")
        mycube_beta[:, :, replicnum - 1] = _core.load_arma_mat_np(str(fpath))
        fpath = scenario_path.joinpath(replicnum_dirname, "stat_delta.txt")
        mycube_delta[:, :, replicnum - 1] = _core.load_arma_mat_np(str(fpath))
    mae_beta = None
    mae_delta = None
    if statistic_type == "mean":
        mae_beta = np.mean(mycube_beta, axis=2) # avg_beta_stat
        mae_delta = np.mean(mycube_delta, axis=2)
    elif statistic_type == "median":
        mae_beta = np.median(mycube_beta, axis=2)
        mae_delta = np.median(mycube_delta, axis=2)
    # do masking
    fpath = scenario_path.joinpath("datagen_params", "datagen_delta.txt")
    datagen_delta = _core.load_arma_umat_np(str(fpath))
    num_to_filter = 0
    idx_datagen_delta = (datagen_delta == num_to_filter)
    datagen_delta_with_nans = np.where(idx_datagen_delta, datagen_delta, np.nan)
    mae_delta_0 = np.ma.array(mae_delta,
                              mask=np.isnan(datagen_delta_with_nans))
    mae_beta_0 = np.ma.array(mae_beta,
                             mask=np.isnan(datagen_delta_with_nans))
    num_to_filter = 1
    idx_datagen_delta = (datagen_delta == num_to_filter)
    datagen_delta_with_nans = np.where(idx_datagen_delta, datagen_delta, np.nan)
    mae_delta_1 = np.ma.array(mae_delta,
                              mask=np.isnan(datagen_delta_with_nans))
    mae_beta_1 = np.ma.array(mae_beta,
                             mask=np.isnan(datagen_delta_with_nans))
    avg_mae_delta_0 = np.mean(mae_delta_0)
    avg_mae_beta_0 = np.mean(mae_beta_0)
    avg_mae_delta_1 = np.mean(mae_delta_1)
    avg_mae_beta_1 = np.mean(mae_beta_1)
    fname = f"avg_of_mae_delta_0_{statistic_type}.txt"
    report_helpers.save_single_value(avg_mae_delta_0, scenario_path, fname)
    fname = f"avg_of_mae_beta_0_{statistic_type}.txt"
    report_helpers.save_single_value(avg_mae_beta_0, scenario_path, fname)
    fname = f"avg_of_mae_delta_1_{statistic_type}.txt"
    report_helpers.save_single_value(avg_mae_delta_1, scenario_path, fname)
    fname = f"avg_of_mae_beta_1_{statistic_type}.txt"
    report_helpers.save_single_value(avg_mae_beta_1, scenario_path, fname)
    return avg_mae_delta_0, avg_mae_beta_0, avg_mae_delta_1, avg_mae_beta_1

### both reports
def generate_report(schema_file_path, other_json_files_path,
                    scenario_path, number_of_replics,
                    statistic_type):
    # load json data for use in functions
    data = json.loads(schema_file_path.read_bytes())
    # load other json files
    jsonfile_src_path = other_json_files_path.joinpath("01_fixed_vals.json")
    data_more = json.loads(jsonfile_src_path.read_bytes())
    data.update(data_more)
    # continue
    T = data["T"]
    effects_list = list()
    fpath = scenario_path.joinpath("effects_table.txt")
    effects_table = _core.load_arma_umat_np(str(fpath))
    effects_list = [str(xx) for xx in effects_table]
    # add H and H_K to data (data consists mostly of contents of the json file)
    H = len(effects_list)
    H_K = np.prod(data['L_k_s'])
    data["H"] = H 
    data["H_K"] = H_K
    if data["covariates"] == "age_assignedsex":
        data["D"] = 3
    results_dict = dict()
    mae_beta_html, avg_of_mae_beta, \
        mae_delta_html, avg_of_mae_delta = calculate_stats_for_beta_and_delta(
            schema_file_path, scenario_path, number_of_replics, data,
            statistic_type)
    results_dict["mae_beta_html"] = mae_beta_html
    results_dict["avg_of_mae_beta"] = avg_of_mae_beta
    results_dict["mae_delta_html"] = mae_delta_html
    results_dict["avg_of_mae_delta"] = avg_of_mae_delta
    avg_theta_stat_html, avg_of_avg_theta_stat = calculate_theta_stats(
        schema_file_path,
        scenario_path,
        number_of_replics,
        data, statistic_type)
    results_dict["avg_theta_stat_html"] = avg_theta_stat_html
    results_dict["avg_of_avg_theta_stat"] = avg_of_avg_theta_stat
    mae_param_html, avg_of_mae_param = calculate_matrix_param_stats(
        schema_file_path,
        scenario_path,
        number_of_replics,
        data,
        statistic_type,
        "lambda",
        data["D"])
    results_dict["mae_lambda_html"] = mae_param_html
    results_dict["avg_of_mae_lambda"] = avg_of_mae_param
    mae_param_html, avg_of_mae_param = calculate_matrix_param_stats(
        schema_file_path,
        scenario_path,
        number_of_replics,
        data,
        statistic_type,
        "Rmat",
        data["K"])
    results_dict["mae_Rmat_html"] = mae_param_html
    results_dict["avg_of_mae_Rmat"] = avg_of_mae_param
    # note that gamma_avg_of_mae_dict may be empty
    gamma_avg_stats, gamma_avg_of_avg_of_mae = calculate_gamma_stats(
        schema_file_path,
        scenario_path,
        number_of_replics,
        data,
        statistic_type)
    results_dict["gamma_avg_stats"] = gamma_avg_stats.tolist()
    results_dict["gamma_avg_of_avg_of_mae"] = gamma_avg_of_avg_of_mae
    avg_of_class_recovery_metric = calc_avg_of_a_metric_value(
        scenario_path, number_of_replics, "stat_class_recovery.txt")
    report_helpers.save_single_value(avg_of_class_recovery_metric,
                            scenario_path, "avg_of_class_recovery_metric.txt")
    # rename these later
    ## note that these are for averages only
    avg_mae_delta_0, avg_mae_beta_0, \
        avg_mae_delta_1, avg_mae_beta_1  = do_beta_delta_subset_calcs(
            schema_file_path, scenario_path, number_of_replics, data,
            statistic_type)
    results_dict["avg_mae_delta_0"] = avg_mae_delta_0
    results_dict["avg_mae_beta_0"] = avg_mae_beta_0
    results_dict["avg_mae_delta_1"] = avg_mae_delta_1
    results_dict["avg_mae_beta_1"] = avg_mae_beta_1
    # render template and write to file
    jinja_env = jinja2.Environment(loader=jinja2.PackageLoader(
        "probitlcm", "templates"))
    template_fname = "scenario_crosssec_template.html"
    template = jinja_env.get_template(template_fname)
    html_out = template.render(
        results_dict=results_dict,
        avg_of_class_recovery_metric=avg_of_class_recovery_metric,
        T=T,
        title="ourtitle")
    fname = f"report_{statistic_type}.html"
    report_path = scenario_path.joinpath(fname)
    with open(report_path, 'wb') as file_:
        file_.write(html_out.encode("utf-8"))
    print("finished writing report")
