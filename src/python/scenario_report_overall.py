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

import pathlib
import jinja2
import json
import numpy as np
import pandas as pd

from probitlcm import report_helpers
from probitlcm import _core

def calculate_stats_for_beta_and_delta(schema_file_path, scenario_path,
                                             number_of_replics, data):
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
    mae_beta = np.mean(mycube_beta, axis=2) # avg_beta_stat
    mae_delta = np.mean(mycube_delta, axis=2)
    avg_of_mae_beta = np.mean(mae_beta)
    avg_of_mae_delta = np.mean(mae_delta)
    fname = "avg_of_mae_beta.txt" # average mae of beta
    report_helpers.save_single_value(avg_of_mae_beta, scenario_path, fname)
    fname = f"avg_of_mae_delta.txt" # average mae of delta
    report_helpers.save_single_value(avg_of_mae_delta, scenario_path,
                                     fname)
    results_dict["mae_beta"] = report_helpers.convert_table_to_html(
        mae_beta, True)
    results_dict["avg_of_mae_beta"] = avg_of_mae_beta
    results_dict["mae_delta"] = report_helpers.convert_table_to_html(
        mae_delta, True)
    results_dict["avg_of_mae_delta"] = avg_of_mae_delta    
    return results_dict

def calculate_theta_stats(schema_file_path, scenario_path,
                          number_of_replics, data):
    # dict to store results
    results_dict = dict()
    # set up data structures from data
    per_effect_M_j_s = data['per_effect_M_j_s']
    q = data['q']
    J = q * len(per_effect_M_j_s)
    H_K = data['H_K']
    fpath = scenario_path.joinpath('M_j_s.txt')
    M_j_s = _core.load_arma_umat_np(str(fpath))
    M_j_s = np.transpose(M_j_s)
    # form of theta matrix will be: rows are replications, each j column is
    #     the sum of the theta_j_stat matrix for that replication
    mymat_theta = np.empty((number_of_replics, J))    
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        for j in range(1, J + 1):
            fpath = scenario_path.joinpath(
                replicnum_dirname, f'stat_theta_j_{j:03}_sum.txt')
            mymat_theta[replicnum - 1, j - 1] = _core.load_arma_mat_np(
                str(fpath))
    # each column must be summed,
    # divided by the number of replications, divided by H_K
    # and then divided by the appropriate M_j
    mymat_theta_colsums = np.sum(mymat_theta, axis = 0)
    avg_theta_stat = mymat_theta_colsums / number_of_replics
    avg_theta_stat = avg_theta_stat / M_j_s
    avg_theta_stat = avg_theta_stat / H_K
    avg_of_avg_theta_stat = np.mean(avg_theta_stat)
    fname = 'avg_of_avg_theta_stat.txt'
    report_helpers.save_single_value(avg_of_avg_theta_stat, scenario_path, fname)
    results_dict["avg_theta_stat"] = report_helpers.convert_table_to_html(
        avg_theta_stat, True)
    results_dict["avg_of_avg_theta_stat"] = avg_of_avg_theta_stat
    return results_dict

def calculate_matrix_param_stats(schema_file_path, scenario_path,
                                 number_of_replics, data,
                                 param_name, n_rows):
    results_dict = dict()
    K = data["K"]
    mycube_param = np.empty((n_rows, K, number_of_replics))
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        fname = "stat_" + param_name + ".txt"
        fpath = scenario_path.joinpath(replicnum_dirname, fname)
        mycube_param[:, :, replicnum - 1] = _core.load_arma_mat_np(str(fpath))
    mae_param = np.mean(mycube_param, axis=2)
    avg_of_mae_param = np.mean(mae_param)
    fname = "avg_of_mae_" + param_name + ".txt"
    report_helpers.save_single_value(avg_of_mae_param,
                                     scenario_path, fname)
    dict_key = "mae_" + param_name
    results_dict[dict_key] = report_helpers.convert_table_to_html(
        mae_param, True)
    dict_key = "avg_of_mae_" + param_name
    results_dict[dict_key] = avg_of_mae_param
    return results_dict

def calculate_gamma_stats(schema_file_path, scenario_path,
                          number_of_replics, data):
    # dict to store results
    results_dict = dict()
    # set up data structures from data
    L_k_s = data["L_k_s"]
    K = data["K"]
    gamma_cubes_dict = dict()
    gamma_avg_of_mae_dict = dict()
    for k in range(1, K + 1):
        if L_k_s[k - 1] == 2:
            continue
        else:
            L_k = L_k_s[k - 1]
            # each value of the dict is a row vector
            # note that (L_k + 1) - 3 = L_k - 2
            gamma_cubes_dict[k] = np.empty((1, L_k - 2, number_of_replics))
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        for k in gamma_cubes_dict.keys():
            fname = f"stat_gamma_{k}.txt"
            fpath = scenario_path.joinpath(replicnum_dirname, fname)
            gamma_cubes_dict[k][:, :, replicnum - 1] = _core.load_arma_mat_np(
                str(fpath))
    for k in gamma_cubes_dict.keys():
        mae_gamma_k = np.mean(gamma_cubes_dict[k], axis=2)
        gamma_avg_of_mae_dict[k] = np.mean(mae_gamma_k)
        fname = f"avg_of_mae_gamma_{k}.txt"
        report_helpers.save_single_value(gamma_avg_of_mae_dict[k],
                                scenario_path, fname)
    return gamma_avg_of_mae_dict

def calc_gamma_avg_of_avg_of_mae(scenario_path,
                                 gamma_avg_of_mae_dict):
    if len(gamma_avg_of_mae_dict.keys()) is 0:
        return np.empty((0)) # an empty vector of length zero
    else:
        gamma_avg_of_avg_stats_vector = np.empty(
            len(gamma_avg_of_mae_dict.keys()))
        i = 0
        for k in gamma_avg_of_mae_dict.keys():
            gamma_avg_of_avg_stats_vector[i] = gamma_avg_of_mae_dict[k]
            i = i + 1
        # gamma_avg_of_avg_of_avg_stats will be a 1 by 1 array
        gamma_avg_of_avg_of_mae = np.mean(gamma_avg_of_avg_stats_vector)
        fname = f"avg_of_avg_of_mae_gamma.txt"
        report_helpers.save_single_value(gamma_avg_of_avg_of_mae,
                                         scenario_path,
                                         fname)
        return gamma_avg_of_avg_of_mae

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
    
def generate_report(schema_file_path, other_json_files_path,
                    scenario_path, number_of_replics):
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
    vars_dict = dict()
    results_dict = calculate_stats_for_beta_and_delta(
        schema_file_path, scenario_path, number_of_replics, data)
    results_dict_theta = calculate_theta_stats(schema_file_path,
                                           scenario_path,
                                           number_of_replics,
                                           data)
    results_dict_lambda = calculate_matrix_param_stats(schema_file_path,
                                                       scenario_path,
                                                       number_of_replics,
                                                       data,
                                                       "lambda",
                                                       data["D"])
    results_dict_Rmat = calculate_matrix_param_stats(schema_file_path,
                                                  scenario_path,
                                                  number_of_replics,
                                                  data,
                                                  "Rmat",
                                                  data["K"])
    # note that gamma_avg_of_mae_dict may be empty
    gamma_avg_of_mae_dict = calculate_gamma_stats(schema_file_path,
                                                        scenario_path,
                                                        number_of_replics,
                                                        data)
    gamma_avg_of_avg_of_mae = calc_gamma_avg_of_avg_of_mae(
        scenario_path, gamma_avg_of_mae_dict)
    # avg_of_metric_values = calc_avg_of_a_metric_value(scenario_path,
    #                                                   number_of_replics,
    #                                                   'metric_value.txt')
    avg_of_class_recovery_metric = calc_avg_of_a_metric_value(
        scenario_path, number_of_replics, "stat_class_recovery.txt")
    report_helpers.save_single_value(avg_of_class_recovery_metric,
                            scenario_path, "avg_of_class_recovery_metric.txt")
    # rename these later
    beta_delta_subsets = dict()
    beta_delta_subsets["avg_of_avg_of_ae_delta_0"] = calc_avg_of_a_metric_value(
        scenario_path, number_of_replics, "stat_delta_corresp_0.txt")
    report_helpers.save_single_value(
        beta_delta_subsets["avg_of_avg_of_ae_delta_0"],
        scenario_path, "avg_of_avg_of_ae_delta_0.txt")
    beta_delta_subsets["avg_of_avg_of_ae_delta_1"] = calc_avg_of_a_metric_value(
        scenario_path, number_of_replics, "stat_delta_corresp_1.txt")
    report_helpers.save_single_value(
        beta_delta_subsets["avg_of_avg_of_ae_delta_1"],
        scenario_path, "avg_of_avg_of_ae_delta_1.txt")
    beta_delta_subsets["avg_of_avg_of_ae_beta_0"] = calc_avg_of_a_metric_value(
        scenario_path, number_of_replics, "stat_beta_corresp_0.txt")
    report_helpers.save_single_value(
        beta_delta_subsets["avg_of_avg_of_ae_beta_0"],
        scenario_path, "avg_of_avg_of_ae_beta_0.txt")
    beta_delta_subsets["avg_of_avg_of_ae_beta_1"] = calc_avg_of_a_metric_value(
        scenario_path, number_of_replics, "stat_beta_corresp_1.txt")
    report_helpers.save_single_value(
        beta_delta_subsets["avg_of_avg_of_ae_beta_1"],
        scenario_path, "avg_of_avg_of_ae_beta_1.txt")
    # render template and write to file
    jinja_env = jinja2.Environment(loader=jinja2.PackageLoader(
        "probitlcm", "templates"))
    template_fname = "scenario_crosssec_template.html"
    template = jinja_env.get_template(template_fname)
    html_out = template.render(
        vars_dict=vars_dict,
        results_dict=results_dict,
        results_dict_theta=results_dict_theta,
        results_dict_lambda=results_dict_lambda,
        results_dict_Rmat=results_dict_Rmat,
        gamma_avg_of_mae_dict=gamma_avg_of_mae_dict,
        gamma_avg_of_avg_of_mae=gamma_avg_of_avg_of_mae,
        # avg_of_metric_values=avg_of_metric_values,
        beta_delta_subsets=beta_delta_subsets,
        avg_of_class_recovery_metric=avg_of_class_recovery_metric,
        T=T,
        title="ourtitle")
    report_path = scenario_path.joinpath("report.html")
    with open(report_path, 'wb') as file_:
        file_.write(html_out.encode("utf-8"))
    print("finished writing report.html")

