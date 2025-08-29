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

## standard library
import argparse
import json
import pathlib
import copy
import tomllib

## other
import numpy as np
import pandas as pd

from probitlcm import _core
from probitlcm import report_helpers

def move_column_inplace(df, col, pos):
    col = df.pop(col)
    df.insert(pos, col.name, col)

def create_and_save_tex_file(scenarios_info_collated, simulation_path,
                             fname, table_num, T):
    df = pd.DataFrame(scenarios_info_collated)
    df.style.hide(axis="index")
    # move "N" column to first position
    move_column_inplace(df, "N", 0)
    move_column_inplace(df, "J", 1)
    move_column_inplace(df, "K", 2)
    move_column_inplace(df, "L_k_s", 3)
    move_column_inplace(df, "rho", 4)
    # rename columns
    df.rename(columns={"rho": "\\(\\rho\\)",
                       "L_k_s": "\\(L\\)",
                       "K": "\\(K\\)",
                       "J": "\\(J\\)",
                       "N": "\\(N\\)"},
              inplace=True)
    # rename certain columns
    if table_num == 1:
        df.rename(columns={"gamma": "\\(\\gamma\\)",
                           "theta": "\\(\\eta\\)",
                           "Rmat": "\\(R\\)",
                           "lambda": "\\(\\lambda\\)"},
                  inplace = True)
    elif table_num == 2:
        df.rename(columns={"beta": "\\(\\beta\\)",
                           "delta": "\\(\\delta\\)",
                           "delta_0": "\\(\\delta_0\\)",
                           "delta_1": "\\(\\delta_1\\)",
                           "beta_0": "\\(\\beta_0\\)",
                           "beta_1": "\\(\\beta_1\\)"},
                  inplace = True)
    simulation_results_file_path_1 = simulation_path.joinpath(fname)
    with open(simulation_results_file_path_1, 'w') as f:
        f.write(df.style.hide(axis="index").format(precision=3).to_latex())

# simulation_name should have the word "simulation" at the beginning,
#     e.g. "simulation_cross_sec"
def build_report_all_scenarios(environ, sim_info_dir_name,
                               total_num_of_scenarios,
                               statistic_type):
    # load config file
    with open("config_simulation.toml", "rb") as fileObj:
        config = tomllib.load(fileObj)
    if args.environ == "laptop":
        process_dir = config['laptop_process_dir']
    elif args.environ == "cluster":
        process_dir = config['cluster_process_dir']
    process_dir = pathlib.Path(process_dir)
    #sim_info_dir_name = pathlib.Path(args.sim_info_dir_name)
    simulation_path = process_dir.joinpath(sim_info_dir_name)
    run_dir = pathlib.Path.cwd()
    other_json_files_path = run_dir.joinpath(sim_info_dir_name, "json_files")
    jsonfile_src_path = other_json_files_path.joinpath("01_fixed_vals.json")
    fixed_vals_json = json.loads(jsonfile_src_path.read_bytes())
    T = fixed_vals_json["T"]
    quants_from_scenario_json = {"N", "rho", "L_k_s", "K", "J"}
    # set up various metrics
    # build table 1
    dict_of_per_scenario_metrics_1 = {
        "theta": f"avg_of_avg_theta_stat_{statistic_type}.txt",
        "Rmat": f"avg_of_mae_Rmat_{statistic_type}.txt",
        "lambda": f"avg_of_mae_lambda_{statistic_type}.txt",
        "alpha n": f"avg_of_class_recovery_metric.txt"
    }
    list_for_table_1 = list()
    for scenario_num in range(1, total_num_of_scenarios + 1):
        # build the dict we want
        ## first get the fixed quantities
        jsonfilename_stem = f"scenario_{scenario_num:04}"
        scenario_path = simulation_path.joinpath(jsonfilename_stem)
        jsonfilename = jsonfilename_stem + ".json"
        jsonfilepath = scenario_path.joinpath(jsonfilename)
        data = json.loads(jsonfilepath.read_bytes())
        dict_table_row = {i: data[i] for i in quants_from_scenario_json}
        ## now get the metric values
        ### only load gamma if it is relevant to particular scenario
        if any(l > 2 for l in dict_table_row["L_k_s"]):
            metric_value = report_helpers.load_single_value(
                scenario_path, f"avg_of_avg_of_mae_gamma_{statistic_type}.txt")
            dict_table_row["gamma"] = metric_value
        else:
            dict_table_row["gamma"] = ""
        for key, value in dict_of_per_scenario_metrics_1.items():
            metric_value = report_helpers.load_single_value(scenario_path,
                                                            value)
            dict_table_row[key] = metric_value
        list_for_table_1.append(dict_table_row)
    create_and_save_tex_file(list_for_table_1,
                             simulation_path,
                             f"simulation_results_1_{statistic_type}.tex",
                             1, T)
    # build table 2
    dict_of_per_scenario_metrics_2 = {
        "beta": f"avg_of_mae_beta_{statistic_type}.txt",
        "delta": f"avg_of_mae_delta_{statistic_type}.txt",
        "delta_0": f"avg_of_mae_delta_0_{statistic_type}.txt",
        "delta_1": f"avg_of_mae_delta_1_{statistic_type}.txt",
        "beta_0": f"avg_of_mae_beta_0_{statistic_type}.txt",
        "beta_1": f"avg_of_mae_beta_1_{statistic_type}.txt"
    }
    list_for_table_2 = list()
    for scenario_num in range(1, total_num_of_scenarios + 1):
        # build the dict we want
        ## first get the fixed quantities
        jsonfilename_stem = f"scenario_{scenario_num:04}"
        scenario_path = simulation_path.joinpath(jsonfilename_stem)
        jsonfilename = jsonfilename_stem + ".json"
        jsonfilepath = scenario_path.joinpath(jsonfilename)
        data = json.loads(jsonfilepath.read_bytes())
        dict_table_row = {i: data[i] for i in quants_from_scenario_json}
        ## now get the metric values        
        for key, value in dict_of_per_scenario_metrics_2.items():
            metric_value = report_helpers.load_single_value(
                scenario_path, value)
            dict_table_row[key] = metric_value
        list_for_table_2.append(dict_table_row)
    create_and_save_tex_file(list_for_table_2,
                             simulation_path,
                             f"simulation_results_2_{statistic_type}.tex",
                             2, T)

def build_report_msce(environ, sim_info_dir_name,
                      total_num_of_scenarios):
    # load config file
    with open("config_simulation.toml", "rb") as fileObj:
        config = tomllib.load(fileObj)
    if args.environ == "laptop":
        process_dir = config['laptop_process_dir']
    elif args.environ == "cluster":
        process_dir = config['cluster_process_dir']
    process_dir = pathlib.Path(process_dir)
    #sim_info_dir_name = pathlib.Path(args.sim_info_dir_name)
    simulation_path = process_dir.joinpath(sim_info_dir_name)
    run_dir = pathlib.Path.cwd()
    other_json_files_path = run_dir.joinpath(sim_info_dir_name, "json_files")
    jsonfile_src_path = other_json_files_path.joinpath("01_fixed_vals.json")
    fixed_vals_json = json.loads(jsonfile_src_path.read_bytes())
    T = fixed_vals_json["T"]
    quants_from_scenario_json = {"N", "rho", "L_k_s", "K", "J"}
    # set up various metrics
    param_names_list = ["theta", "beta", "delta", "gamma",
                        "lambda", "Rmat"]
    dict_of_per_scenario_metrics = dict()
    for param_name in param_names_list:
        dict_of_per_scenario_metrics[param_name] = f"avg_mcse_{param_name}.txt"
    list_for_table = list()
    for scenario_num in range(1, total_num_of_scenarios + 1):
        # build the dict we want
        ## first get the fixed quantities
        jsonfilename_stem = f"scenario_{scenario_num:04}"
        scenario_path = simulation_path.joinpath(jsonfilename_stem)
        jsonfilename = jsonfilename_stem + ".json"
        jsonfilepath = scenario_path.joinpath(jsonfilename)
        data = json.loads(jsonfilepath.read_bytes())
        dict_table_row = {i: data[i] for i in quants_from_scenario_json}
        list_for_table.append(dict_table_row)
        ## now get the metric values        
        for key, value in dict_of_per_scenario_metrics.items():
            metric_value = report_helpers.load_single_value(
                scenario_path, value)
            dict_table_row[key] = metric_value
        list_for_table.append(dict_table_row)
    # make data frame and save as latex table
    df = pd.DataFrame(list_for_table)
    df.style.hide(axis="index")
    # move "N" column to first position
    move_column_inplace(df, "N", 0)
    move_column_inplace(df, "J", 1)
    move_column_inplace(df, "K", 2)
    move_column_inplace(df, "L_k_s", 3)
    move_column_inplace(df, "rho", 4)
    # rename columns
    df.rename(columns={"rho": "\\(\\rho\\)",
                       "L_k_s": "\\(L\\)",
                       "K": "\\(K\\)",
                       "J": "\\(J\\)",
                       "N": "\\(N\\)"},
              inplace=True)
    df.rename(columns={"theta": "\\(\\eta\\)",
                       "beta": "\\(\\beta\\)",
                       "delta": "\\(\\delta\\)",
                       "gamma": "\\(\\gamma\\)",
                       "lambda": "\\(\\lambda\\)",
                       "Rmat": "\\(R\\)"},
              inplace = True)
    fname = "simulation_results_mcse.tex"
    fpath = simulation_path.joinpath(fname)
    with open(fpath, 'w') as f:
        f.write(df.style.hide(axis="index").format(precision=3).to_latex())

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", choices = ["laptop", "cluster"],
                        required=True)
    parser.add_argument("--sim_info_dir_name")
    parser.add_argument('--total_num_of_scenarios', required=True)
    parser.add_argument("--statistic_type",
                        choices = ["mean", "median"],required=True)
    args = parser.parse_args()
    total_num_of_scenarios = int(args.total_num_of_scenarios)
    build_report_all_scenarios(args.environ, args.sim_info_dir_name,
                               total_num_of_scenarios,
                               args.statistic_type)
    build_report_msce(args.environ, args.sim_info_dir_name,
                      total_num_of_scenarios)
    print("Generated overall reports.")
