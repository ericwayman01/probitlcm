# This file is part of "probitlcm" which is released under GPL v3.
#
# Copyright (c) 2022-2024 Eric Alan Wayman <ewayman2@illinois.edu>.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import pandas as pd
import jinja2

from probitlcm import report_helpers

from probitlcm import _core

# Define array calculation functions

def check_interval_membership(x, interval_lb, interval_ub):
    if x >= interval_lb and x <= interval_ub:
        return 1
    else:
        return 0
v_check_interval_membership = np.vectorize(check_interval_membership)

def count_matching_delta_entries(datagen_delta, draws_delta, burnin):
   num_matches_table = np.empty((datagen_delta.shape[0],
                                 datagen_delta.shape[1]))
   for i in range(0, datagen_delta.shape[0]):
      for j in range(0, datagen_delta.shape[1]):
         if datagen_delta[i, j] == 0:
            num_matches_table[i, j] = np.count_nonzero(
                draws_delta[i, j, burnin:] == 0)
         else:
            num_matches_table[i, j] = np.count_nonzero(
                draws_delta[i, j, burnin:])
   return num_matches_table

# Define report generating functions

def calculate_beta_statistics(draws_beta, datagen_beta, avg_beta,
                              scenario_path, replic_path, burnin,
                              effects_list):
    ## beta_recov_table
    # compare
    quantile_lower_bound = 0.05
    quantile_upper_bound = 0.95
    # calculate quantiles of empirical cdf for each coordinate of beta
    quantiles = np.quantile(draws_beta[:, :, burnin:],
                            (quantile_lower_bound, quantile_upper_bound), 2)
    quantiles_lb = quantiles[0, :, :]
    quantiles_ub = quantiles[1, :, :]
    # check if data generating beta falls into the above quantiles or not
    beta_recov_table = v_check_interval_membership(datagen_beta,
                                                   quantiles_lb, quantiles_ub)
    beta_recov_count = np.count_nonzero(beta_recov_table)
    beta_recov_size = beta_recov_table.size
    # datagenbeta
    datagen_beta_df = pd.DataFrame(datagen_beta)
    datagen_beta_df.index = effects_list
    beta_dict = dict()
    beta_dict["datagen_beta"] = datagen_beta_df.to_html(
        float_format='{:10.2f}'.format)
    beta_dict["beta_recov_table"] = report_helpers.convert_table_to_html(
        beta_recov_table, False)
    beta_dict["beta_recov_count"] = beta_recov_count
    beta_dict["beta_recov_size"] = beta_recov_size
    betahat_avg_error = np.abs(avg_beta - datagen_beta)
    fname = "stat_beta.txt"
    fpath = replic_path.joinpath(fname)
    _core.save_arma_mat_np(betahat_avg_error, str(fpath))
    beta_dict["betahat_avg_error"] = report_helpers.convert_table_to_html(
        betahat_avg_error, True)
    return beta_dict    

def calculate_delta_statistics(draws_delta, datagen_delta, avg_delta,
                               scenario_path, replic_path, burnin,
                               effects_list):
    # perepare count matching delta entries table
    num_matches_table = count_matching_delta_entries(datagen_delta,
                                                     draws_delta, burnin)
    # calculate abs(delta_hat - delta)
    deltahat = np.int64(avg_delta > 0.5)
    deltahat_avg_error = np.where(deltahat == datagen_delta, 1, 0)
    deltahat_avg_error = deltahat_avg_error.astype(np.float64)
    # save stat
    fname = "stat_delta.txt"
    fpath = replic_path.joinpath(fname)
    _core.save_arma_mat_np(deltahat_avg_error, str(fpath))
    delta_dict = dict()
    delta_dict["num_matches_table"] = pd.DataFrame(num_matches_table).to_html()
    delta_dict["datagen_delta"] = pd.DataFrame(datagen_delta).to_html()
    delta_dict["deltahat_avg_error"] = report_helpers.convert_table_to_html(
        deltahat_avg_error, False)
    return delta_dict

def calc_delta_beta_avg_of_ae_subset(datagen_delta, avg_delta,
                                     datagen_beta, avg_beta,
                                     num_to_filter):
    idx_datagen_delta = (datagen_delta == num_to_filter)
    # calculate delta_avg_mae
    datagen_delta = np.where(idx_datagen_delta, datagen_delta, np.nan)
    # avg_delta_corresp = np.where(idx_datagen_delta, avg_delta, np.nan)
    datagen_delta_masked = np.ma.array(datagen_delta,
                                       mask=np.isnan(datagen_delta))
    avg_delta_masked = np.ma.array(avg_delta, mask=np.isnan(datagen_delta))
    deltahat_corresp = np.int64(avg_delta_masked > 0.5)
    delta_corresp_result = (deltahat_corresp == datagen_delta_masked)
    delta_avg_mae = np.sum(
        delta_corresp_result) / delta_corresp_result.count()
    # calculate beta_avg_mae
    datagen_beta_masked = np.ma.array(datagen_beta,
                                      mask=np.isnan(datagen_delta))
    avg_beta_masked = np.ma.array(avg_beta, mask=np.isnan(datagen_delta))
    avg_deviat_beta = np.abs(datagen_beta_masked - avg_beta_masked)
    beta_avg_mae = np.sum(avg_deviat_beta) / avg_deviat_beta.count()
    return (delta_avg_mae, beta_avg_mae)

def calculate_beta_and_delta_stats(scenario_path,
                                   scenario_datagen_params_path,
                                   replic_path,
                                   burnin, effects_list):
    # load everything
    ## load delta
    fpath = replic_path.joinpath("draws_delta.txt")
    draws_delta = _core.load_arma_ucube_np(str(fpath))
    fpath = scenario_datagen_params_path.joinpath("datagen_delta.txt")
    datagen_delta = _core.load_arma_umat_np(str(fpath))
    ## load beta
    fpath = replic_path.joinpath("draws_beta.txt")
    draws_beta = _core.load_arma_cube_np(str(fpath))
    fpath = scenario_datagen_params_path.joinpath("datagen_beta.txt")
    datagen_beta = _core.load_arma_mat_np(str(fpath))
    fpath = replic_path.joinpath("average_beta.txt")
    avg_beta = _core.load_arma_mat_np(str(fpath))
    fpath = replic_path.joinpath("average_delta.txt")
    avg_delta = _core.load_arma_mat_np(str(fpath))
    ## beta dict
    beta_dict = calculate_beta_statistics(draws_beta, datagen_beta, avg_beta,
                                          scenario_path, replic_path, burnin,
                                          effects_list)
    ## delta dict

    delta_dict = calculate_delta_statistics(draws_delta, datagen_delta,
                                            avg_delta,
                                            scenario_path, replic_path, burnin,
                                            effects_list)
    # calculate avg_of_ae's
    delta_avg_of_ae_0, beta_avg_of_ae_0 = calc_delta_beta_avg_of_ae_subset(
        datagen_delta, avg_delta, datagen_beta, avg_beta, 0)
    delta_avg_of_ae_1, beta_avg_of_ae_1 = calc_delta_beta_avg_of_ae_subset(
        datagen_delta, avg_delta, datagen_beta, avg_beta, 1)
    # save ears
    report_helpers.save_single_value(delta_avg_of_ae_0, replic_path,
                            "stat_delta_corresp_0.txt")
    report_helpers.save_single_value(beta_avg_of_ae_0, replic_path,
                            "stat_beta_corresp_0.txt")
    report_helpers.save_single_value(delta_avg_of_ae_1, replic_path,
                            "stat_delta_corresp_1.txt")
    report_helpers.save_single_value(beta_avg_of_ae_1, replic_path,
                            "stat_beta_corresp_1.txt")    
    # append avg ae's
    delta_dict["delta_avg_of_ae_0"] = delta_avg_of_ae_0
    delta_dict["delta_avg_of_ae_1"] = delta_avg_of_ae_1
    beta_dict["beta_avg_of_ae_0"] = beta_avg_of_ae_0
    beta_dict["beta_avg_of_ae_1"] = beta_avg_of_ae_1
    return (delta_dict, beta_dict)

def calculate_theta_statistics(scenario_path, scenario_datagen_params_path,
                               replic_path, burnin, J):
    theta_mat_stats_list = list()
    for j in range(1, J + 1):
        # load avg mat and true mat
        fpath = replic_path.joinpath(f'theta_j_mat_avg_{j:03}.txt')
        theta_j_hat = _core.load_arma_mat_np(str(fpath))
        fpath = scenario_datagen_params_path.joinpath(
            f'datagen_theta_j_mat_{j:03}.txt')
        datagen_theta_j = _core.load_arma_mat_np(str(fpath))
        # subtract them and take the abs value
        theta_j_stat = np.abs(theta_j_hat - datagen_theta_j)
        fname = f'stat_theta_j_{j:03}.txt'
        fpath = replic_path.joinpath(fname)
        _core.save_arma_mat_np(theta_j_stat, str(fpath))
        theta_mat_stats_list.append(report_helpers.convert_table_to_html(
            theta_j_stat, True))
        # sum up values for each J and store them (so when we do the
        #     overall report, we don't have to know the dimensions of each
        #     theta_j_stat matrix to build each cube
        theta_j_stat_sum = np.sum(theta_j_stat)
        fname = f'stat_theta_j_{j:03}_sum.txt'
        report_helpers.save_single_value(theta_j_stat_sum, replic_path, fname)
    # return the list to be printed in the template
    return theta_mat_stats_list

def calculate_simple_average_param_stat(scenario_datagen_params_path,
                                        replic_path,
                                        burnin, param_name):
    # load draws of lambda and datagen lambda
    fname = f'average_{param_name}.txt'
    fpath = replic_path.joinpath(fname)
    avg_param = _core.load_arma_mat_np(str(fpath))
    fname = f'datagen_{param_name}.txt'
    fpath = scenario_datagen_params_path.joinpath(fname)
    datagen_param = _core.load_arma_mat_np(str(fpath))
    param_stat = np.abs(avg_param - datagen_param)
    fname = f'stat_{param_name}.txt'
    fpath = replic_path.joinpath(fname)
    _core.save_arma_mat_np(param_stat, str(fpath))        
    return report_helpers.convert_table_to_html(param_stat, True)

def calculate_gamma_statistics(scenario_datagen_params_path, replic_path,
                               burnin, K, L_k_s):
    gamma_stats_dict = dict()
    fname = ''
    p = ''
    for k in range(1, K + 1):
        if L_k_s[k - 1] == 2:
            continue
        else:
            L_k = L_k_s[k - 1]
            # load row vector for gamma empirical average and datagen gamma
            fname = f"average_gamma_{k}.txt"
            fpath = replic_path.joinpath(fname)
            avg_param = _core.load_arma_mat_np(str(fpath))
            fname = f"datagen_gamma_{k}.txt"
            fpath = scenario_datagen_params_path.joinpath(fname)
            datagen_param = _core.load_arma_mat_np(str(fpath))
            # remove first two and last dimensions
            avg_param = np.delete(avg_param, (0, 1, L_k), 1)
            datagen_param = np.delete(datagen_param, (0, 1, L_k), 1)
            # param_stat is a row vector
            param_stat = np.abs(avg_param - datagen_param)
            fname = f"stat_gamma_{k}.txt"
            fpath = replic_path.joinpath(fname)
            _core.save_arma_mat_np(param_stat, str(fpath))        
            gamma_stats_dict[k] = report_helpers.convert_table_to_html(
                param_stat, True)
    return gamma_stats_dict

def calculate_class_count_stats(scenario_path, replic_path):
    # load tables
    fname = "class_counts.txt"
    fpath = replic_path.joinpath(fname)
    class_counts = _core.load_arma_umat_np(str(fpath))
    fname = "datagen_class_numbers.txt"
    fpath = replic_path.joinpath(fname)
    class_datagen = _core.load_arma_umat_np(str(fpath)).ravel()
    # now form tables
    C = class_counts.shape[1]
    class_report_table = np.empty((C, C))
    for true_class in range(C):
        # find all rows (subjects) for which the class_datagen is a particular
        #     class ("true class"), then take the argmax over axis = 1 (each
        #     row). this finds, for each subject, the mode class. that's joe.
        set_of_rows = np.argmax(class_counts[class_datagen == true_class, :],
                                axis=1)
        for c in range(C):
            class_report_table[true_class, c] = np.count_nonzero(
                set_of_rows == c)
    class_recovery_metric = np.trace(class_report_table) / \
        np.sum(class_report_table)
    fname = "stat_class_recovery.txt"
    report_helpers.save_single_value(class_recovery_metric, replic_path, fname)
    class_report_table = report_helpers.convert_table_to_html(
                class_report_table, True)
    return (class_report_table, class_recovery_metric)

def pretty_print_2d_ndarray(myarray):
   print(np.array_str(myarray, precision=3, suppress_small=True))

def generate_report(scenario_path,
                    scenario_datagen_params_path,
                    replic_path, data, replicnum):
    burnin = data["burnin"]
    chain_length_after_burnin = data["chain_length_after_burnin"]
    T = data["T"]
    order = data["order"]
    K = data["K"]
    effects_list = list()
    fpath = str(scenario_path.joinpath("effects_table.txt"))
    effects_table = _core.load_arma_umat_np(fpath)
    effects_list = [str(xx) for xx in effects_table]
    J = len(data["per_effect_M_j_s"]) * data["q"]
    delta_dict, beta_dict = calculate_beta_and_delta_stats(
        scenario_path, scenario_datagen_params_path, replic_path,
        burnin, effects_list)
    theta_mat_stats_list = calculate_theta_statistics(
        scenario_path, scenario_datagen_params_path, replic_path, burnin, J)
    class_report_table, class_recovery_metric = calculate_class_count_stats(
        scenario_path, replic_path)
    lambda_stat = calculate_simple_average_param_stat(
        scenario_datagen_params_path, replic_path, burnin, "lambda")
    Rmat_stat = calculate_simple_average_param_stat(
        scenario_datagen_params_path, replic_path, burnin, "Rmat")
    # note that the gamma_stats_dict may be empty (with
    #     no gamma_stat_{k}.npy having been saved either)
    gamma_stats_dict = calculate_gamma_statistics(scenario_datagen_params_path,
                                                  replic_path,
                                                  burnin,
                                                  data["K"], data["L_k_s"])
    # fpath = replic_path.joinpath('metric_value.txt')
    # metric_value_mat = _core.load_arma_mat_np(str(fpath))
    # metric_value = metric_value_mat.item()
    vars_dict = dict()
    vars_dict["replicnum"] = replicnum
    jinja_env = jinja2.Environment(loader=jinja2.PackageLoader(
        "probitlcm", "templates"))
    template_fname = "replication_crosssec_template.html"
    template = jinja_env.get_template(template_fname)
    html_out = template.render(vars_dict=vars_dict,
                               beta_dict=beta_dict,
                               delta_dict=delta_dict,
                               theta_mat_stats_list=theta_mat_stats_list,
                               lambda_stat=lambda_stat,
                               Rmat_stat=Rmat_stat,
                               gamma_stats_dict=gamma_stats_dict,
                               class_report_table=class_report_table,
                               class_recovery_metric=class_recovery_metric,
                               # metric_value=metric_value,
                               T=T,
                               title="thetitle")
    report_path = replic_path.joinpath("report.html")
    with open(report_path, 'wb') as file_:
        file_.write(html_out.encode("utf-8"))
