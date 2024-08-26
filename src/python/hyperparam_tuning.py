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
import shutil
import argparse
import subprocess
import json
import itertools

import time

import sys
if sys.version_info[1] < 11:
    import toml
else:
    import tomllib as toml

from probitlcm.scenario_launch import scenario_launch_setup
from probitlcm.scenario_launch import scenario_launch_laptop

from probitlcm.data_analysis_run import data_analysis_run_laptop
from probitlcm.data_analysis_run import data_analysis_setup

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", choices = ["laptop"],
                        required=True)
    parser.add_argument("--dir_name", required=True)
    parser.add_argument("--type_of_run",
                        choices = ["simulation", "data_analysis"],
                        required=True)
    parser.add_argument("--scenario_or_setup_number", required=True)
    # parse arguments
    args = parser.parse_args()
    environ = args.environ
    dir_name = args.dir_name
    type_of_run = args.type_of_run
    scenario_or_setup_number = args.scenario_or_setup_number
    if type_of_run == "simulation":
        # run scenario_launch_setup
        scenario_path, jsonfilename_stem, process_dir, \
            other_json_files_path, number_of_replics, \
            scenarionumber_zb = scenario_launch_setup(environ,
                                                      dir_name,
                                                      scenario_or_setup_number)
        # load hyperparams
        print(other_json_files_path)
        fpath = other_json_files_path.joinpath(
            "01_list_hyperparams_possible_vals.json")
        dict_hyperparams_possible_vals = json.loads(fpath.read_text())
        # note: right now this only supports a single hyperparameter
        # run replications
        hyperparam_tuning = True
        tuning_path = process_dir.joinpath("hyperparam_tuning")
        tuning_path.mkdir(parents=True, exist_ok=True)
        fname = "02_current_tuning_hyperparam_vals.json"
        fpath = other_json_files_path.joinpath(fname)
        thekey = list(dict_hyperparams_possible_vals.keys())[0]
        thevalues = dict_hyperparams_possible_vals[thekey]
        # for loop
        for val in thevalues:
            mydict = dict()
            mydict[thekey] = val
            with fpath.open(mode='w') as f:
                json.dump(mydict, f, indent=4)
            number_of_replics = 1
            # now launch scenario
            scenario_launch_laptop(jsonfilename_stem, process_dir,
                                   other_json_files_path,
                                   number_of_replics, scenarionumber_zb,
                                   hyperparam_tuning,
                                   tuning_path)
            # delete the replic_path before returning to top of loop
            replicnum = 1
            replicnum_dirname = f"replic_{replicnum:03}"
            replic_path = scenario_path.joinpath(replicnum_dirname)
            shutil.rmtree(replic_path)    
    else: # type_of_run == "data_analysis"
        setup_num, dataset_dir, \
            process_dir, data_analysis_path, \
            number_of_chains = data_analysis_setup(
                environ, dir_name, scenario_or_setup_number)
        fpath = dataset_dir.joinpath(
            "01_list_hyperparams_possible_vals.json")
        dict_hyperparams_possible_vals = json.loads(fpath.read_text())
        hyperparam_tuning = True
        tuning_path = process_dir.joinpath("hyperparam_tuning")
        print("data_analysis_path:", data_analysis_path)
        tuning_path.mkdir(parents=True, exist_ok=True)
        fname = "02_current_tuning_hyperparam_vals.json"
        fpath = dataset_dir.joinpath(fname)
        thekey = list(dict_hyperparams_possible_vals.keys())[0]
        thevalues = dict_hyperparams_possible_vals[thekey]
        # for loop
        for val in thevalues:
            mydict = dict()
            mydict[thekey] = val
            with fpath.open(mode='w') as f:
                json.dump(mydict, f, indent=4)
            number_of_chains = 1
            # now launch data_analysis
            data_analysis_run_laptop(setup_num, dataset_dir,
                                     data_analysis_path, number_of_chains,
                                     hyperparam_tuning, tuning_path)            
            # delete the chain_path before returning to top of loop
            chainnum = 1
            chainnum_dirname = f"chain_{chainnum:03}"
            chain_path = data_analysis_path.joinpath(chainnum_dirname)
            shutil.rmtree(chain_path)
