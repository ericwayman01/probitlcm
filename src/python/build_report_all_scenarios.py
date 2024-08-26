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

## standard library
import argparse
import json
import pathlib
import copy

## other
import numpy as np
import pandas as pd

import sys
if sys.version_info[1] < 11:
    import toml
else:
    import tomllib as toml

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
                           "avg_of_avg_of_ae_delta_0": "\\(\\delta_0\\)",
                           "avg_of_avg_of_ae_delta_1": "\\(\\delta_1\\)",
                           "avg_of_avg_of_ae_beta_0": "\\(\\beta_0\\)",
                           "avg_of_avg_of_ae_beta_1": "\\(\\beta_1\\)"},
                  inplace = True)
    simulation_results_file_path_1 = simulation_path.joinpath(fname)
    with open(simulation_results_file_path_1, 'w') as f:
        f.write(df.style.hide(axis="index").format(precision=3).to_latex())

# simulation_name should have the word "simulation" at the beginning,
#     e.g. "simulation_cross_sec"
def build_report_all_scenarios(environ, sim_info_dir_name,
                               total_num_of_scenarios):
    # load config file
    with open('config_simulation.toml') as fileObj:
        config = toml.load(fileObj)
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
        "theta": "avg_of_avg_theta_stat.txt",
        "Rmat": "avg_of_mae_Rmat.txt"
        # "classrecov": "avg_of_class_recovery_metric.txt",
    }
    dict_of_per_scenario_metrics_1["lambda"] = "avg_of_mae_lambda.txt"
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
                scenario_path, 'avg_of_avg_of_mae_gamma.txt')
            dict_table_row["gamma"] = metric_value
        else:
            dict_table_row["gamma"] = ""
        for key, value in dict_of_per_scenario_metrics_1.items():
            metric_value = report_helpers.load_single_value(scenario_path,
                                                            value)
            dict_table_row[key] = metric_value
        list_for_table_1.append(dict_table_row)
    create_and_save_tex_file(list_for_table_1,
                             simulation_path, 'simulation_results_1.tex',
                             1, T)
    # build table 2
    dict_of_per_scenario_metrics_2 = {
        "beta": "avg_of_mae_beta.txt",
        "delta": "avg_of_mae_delta.txt",
        "avg_of_avg_of_ae_delta_0": "avg_of_avg_of_ae_delta_0.txt",
        "avg_of_avg_of_ae_delta_1": "avg_of_avg_of_ae_delta_1.txt",
        "avg_of_avg_of_ae_beta_0": "avg_of_avg_of_ae_beta_0.txt",
        "avg_of_avg_of_ae_beta_1": "avg_of_avg_of_ae_beta_1.txt"
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
                             simulation_path, "simulation_results_2.tex",
                             2, T)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", choices = ["laptop", "cluster"],
                        required=True)
    parser.add_argument("--sim_info_dir_name")
    parser.add_argument('--total_num_of_scenarios', required=True)
    args = parser.parse_args()
    total_num_of_scenarios = int(args.total_num_of_scenarios)
    build_report_all_scenarios(args.environ, args.sim_info_dir_name,
                               total_num_of_scenarios)
    print("Generated overall reports.")
