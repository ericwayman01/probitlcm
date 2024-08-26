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

# for template file import
import importlib.resources as pkg_resources
from probitlcm import templates
from probitlcm import report_helpers
from probitlcm import _core
from probitlcm.scenario_report_single_replic import check_interval_membership

import pandas as pd
import numpy as np
import jinja2

v_check_interval_membership = np.vectorize(check_interval_membership)

def build_table_of_credible_intervals(quantiles_lb, quantiles_ub):
    x_len = quantiles_lb.shape[0]
    y_len = quantiles_lb.shape[1]
    list_of_lists = list()
    for x in range(x_len):
        myrow = list()
        for y in range(y_len):
            myrow.append((quantiles_lb[x, y], quantiles_ub[x, y]))
        list_of_lists.append(myrow)
    return pd.DataFrame(list_of_lists)

def calculate_theta_statistics(chain_results_path, J):
    theta_j_hats_list = list()
    for j in range(1, J + 1):
        # load avg mat and true mat
        fpath = replic_path.joinpath(f'theta_j_mat_avg_{j:03}.txt')
        theta_j_hat = _core.load_arma_mat_np(str(fpath))
        theta_mat_stats_list.append(report_helpers.convert_table_to_html(
            theta_j_stat, True))
    return theta_j_hats_list

def calculate_gamma_statistics(chain_results_path, K, L_k_s):
    gamma_estimates_dict = dict()
    fname = ''
    p = ''
    for k in range(1, K + 1):
        if L_k_s[k - 1] == 2:
            continue
        else:
            L_k = L_k_s[k - 1]
            # load row vector for gamma empirical average
            fname = f'average_gamma_{k}.txt'
            fpath = chain_results_path.joinpath(fname)
            avg_param = _core.load_arma_mat_np(str(fpath))
            # remove first two and last dimensions
            avg_param = np.delete(avg_param, (0, 1, L_k), 1)
            # param_stat is a row vector
            gamma_estimates_dict[k] = report_helpers.convert_table_to_html(
                avg_param, True)
    return gamma_estimates_dict

def do_estimates_work_beta_delta(effects_dict, dataset_dir,
                                 data_analysis_path,
                                 chain_results_path, burnin):
    fname = "draws_beta.txt"
    fpath = chain_results_path.joinpath(fname)
    draws_beta = _core.load_arma_cube_np(str(fpath))
    fname = "average_beta.txt"
    fpath = chain_results_path.joinpath(fname)
    average_beta = _core.load_arma_mat_np(str(fpath))
    # compare
    quantile_lower_bound = 0.025
    quantile_upper_bound = 0.975
    # calculate quantiles of empirical cdf for each coordinate of beta
    quantiles = np.quantile(draws_beta[:, :, burnin:],
                            (quantile_lower_bound, quantile_upper_bound), 2)
    quantiles_lb = quantiles[0, :, :]
    quantiles_ub = quantiles[1, :, :]
    beta_recov_table = v_check_interval_membership(0,
                                                   quantiles_lb, quantiles_ub)
    my_mask = np.ma.masked_where(beta_recov_table == 1, average_beta)
    df = pd.DataFrame(my_mask)
    df = df.replace(np.nan, "")
    # load colnames
    responses_df = pd.read_csv(dataset_dir.joinpath("responses.csv"))
    responses_colnames = list(responses_df.columns)
    # add them to table
    df.columns = responses_colnames
    # split into two tables
    # df_one = df.iloc[:, 0:9]
    # df_two = df.iloc[:, 9:]
    # df_one.insert(0, "effect", effects_dict["effects_list"])
    # df_two.insert(0, "effect", effects_dict["effects_list"])
    df.insert(0, "effect", effects_dict["effects_list"])
    # save result
    # path_results_tex = chain_results_path.joinpath(
    #     "results_beta_delta_with_sparsity_1.tex")
    # with open(path_results_tex, 'w') as f:
    #     f.write(df_one.to_latex(float_format="{:.2f}".format))
    # path_results_tex = chain_results_path.joinpath(
    #     "results_beta_delta_with_sparsity_2.tex")
    # with open(path_results_tex, 'w') as f:
    #     f.write(df_two.to_latex(float_format="{:.2f}".format))
    # save sparse table
    df.to_pickle(
        chain_results_path.joinpath(
            "results_beta_delta_with_sparsity_pandas_df.xz"),
        compression="xz")
    # build and save credible intervals table
    df = build_table_of_credible_intervals(quantiles_lb, quantiles_ub)
    df.to_pickle(
        chain_results_path.joinpath(
            "results_beta_credible_intervals_pandas_df.xz"),
        compression="xz")

def do_estimates_work_lambda(effects_dict,
                             param_estimates_dict, chain_results_path,
                             burnin):
    param_name = "lambda";
    fname = "average_" + param_name + ".txt"
    fpath = chain_results_path.joinpath(fname)
    my_df = pd.DataFrame(_core.load_arma_mat_np(str(fpath)))
    my_df = my_df.round(decimals = 2)
    param_estimates_dict[param_name] = my_df.to_html(
        float_format='{:10.2f}'.format)
    fname = "results_" + param_name + ".tex"
    path_results_tex = chain_results_path.joinpath(fname)
    with open(path_results_tex, 'w') as f:
        f.write(my_df.to_latex(float_format="{:.2f}".format))
    # build table of credible intervals
    fname = "draws_lambda.txt"
    fpath = chain_results_path.joinpath(fname)
    draws_lambda = _core.load_arma_cube_np(str(fpath))
    quantile_lower_bound = 0.025
    quantile_upper_bound = 0.975
    quantiles = np.quantile(draws_lambda[:, :, burnin:],
                            (quantile_lower_bound, quantile_upper_bound), 2)
    quantiles_lb = quantiles[0, :, :]
    quantiles_ub = quantiles[1, :, :]
    df = build_table_of_credible_intervals(quantiles_lb, quantiles_ub)
    df.to_pickle(
        chain_results_path.joinpath(
            "results_lambda_credible_intervals_pandas_df.xz"),
        compression="xz")

def do_estimates_work(param_name, effects_dict,
                      param_estimates_dict, chain_results_path,
                      use_colnames, responses_colnames):
    fname = "average_" + param_name + ".txt"
    fpath = chain_results_path.joinpath(fname)
    my_df = pd.DataFrame(_core.load_arma_mat_np(str(fpath)))
    if use_colnames is True:
        my_df.columns = responses_colnames
    # custom per-parameter settings
    if param_name in ["beta", "delta"]:
        my_df.insert(0, "effect", effects_dict["effects_list"])
    my_df = my_df.round(decimals = 2)
    param_estimates_dict[param_name] = my_df.to_html(
        float_format='{:10.2f}'.format)
    fname = "results_" + param_name + ".tex"
    path_results_tex = chain_results_path.joinpath(fname)
    with open(path_results_tex, 'w') as f:
        f.write(my_df.to_latex(float_format="{:.2f}".format))

def generate_report(dataset_dir, data_analysis_path, chain_results_path,
                    data, chainnum):
    burnin = data["burnin"]
    chain_length_after_burnin = data["chain_length_after_burnin"]
    T = data["T"]
    order = data["order"]
    K = data["K"]
    effects_dict = dict()
    fpath = data_analysis_path.joinpath("effects_table.txt")
    effects_table = _core.load_arma_umat_np(str(fpath))
    effects_dict["effects_list"] = [str(xx) for xx in effects_table]
    H = len(effects_dict["effects_list"])
    ## for transition model
    H_otr = H
    data["H"] = H
    data["H_otr"] = H_otr
    # create dicts for stuff
    vars_dict = dict()
    param_estimates_dict = dict()
    burnin = data["burnin"]
    # load column names from responses.csv
    responses_df = pd.read_csv(dataset_dir.joinpath('responses.csv'))
    responses_colnames = list(responses_df.columns)
    # load effects table for labels
    # do some estimates work
    do_estimates_work("beta", effects_dict,
                      param_estimates_dict, chain_results_path,
                      True, responses_colnames)
    do_estimates_work("delta", effects_dict,
                      param_estimates_dict, chain_results_path,
                      True, responses_colnames)
    do_estimates_work_beta_delta(effects_dict, dataset_dir, data_analysis_path,
                                 chain_results_path, burnin)
    do_estimates_work_lambda(effects_dict,
                             param_estimates_dict, chain_results_path,
                             burnin)
    do_estimates_work("Rmat", effects_dict,
                      param_estimates_dict, chain_results_path,
                      False, responses_colnames)
    J = responses_df.shape[1]
    L_k_s = data["L_k_s"]
    gamma_estimates_dict = calculate_gamma_statistics(chain_results_path,
                                                      K, L_k_s)
    # render template and write to file
    jinja_env = jinja2.Environment(loader=jinja2.PackageLoader(
        "probitlcm", "templates"))
    template = jinja_env.get_template("chain_template.html")
    html_out = template.render(vars_dict=vars_dict,
                               param_estimates_dict=param_estimates_dict,
                               gamma_estimates_dict=gamma_estimates_dict)
    report_path = chain_results_path.joinpath("report.html")
    with open(report_path, 'wb') as file_:
        file_.write(html_out.encode("utf-8"))
