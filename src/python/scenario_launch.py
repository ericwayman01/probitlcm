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
import pathlib
import shutil
import argparse
import subprocess
import json
import itertools

## other modules
import sys
if sys.version_info[1] < 11:
    import toml
else:
    import tomllib as toml

import numpy as np

## my modules
from probitlcm import scenario_run_single_replic
from probitlcm import report_helpers
from probitlcm import run_helpers

from probitlcm import _core

def scenario_launch_setup(environ, sim_info_dir_name,
                          scenarionumber):
    # set up paths
    run_dir = pathlib.Path.cwd()
    # load config file
    with open("config_simulation.toml") as fileObj:
        config = toml.load(fileObj)
    if environ == "laptop":
        process_dir = config["laptop_process_dir"]
    elif environ == "cluster":
        process_dir = config["cluster_process_dir"]
    process_dir = pathlib.Path(process_dir)
    sim_info_dir_name = pathlib.Path(sim_info_dir_name)
    process_dir = process_dir.joinpath(sim_info_dir_name)
    process_dir.mkdir(parents=True, exist_ok=True)
    ## run all replications for this particular scenario
    scenarionumber = int(scenarionumber)
    # zero-based scenario number for generating seeds
    scenarionumber_zb = scenarionumber - 1
    jsonfilename_stem = f"scenario_{scenarionumber:04}"
    number_of_replics = config['number_of_replics']
    seed_value_scenario = number_of_replics * scenarionumber_zb
    # copy all files over to proper process_dir
    param_files_path = run_dir.joinpath(sim_info_dir_name, "scenario_files",
                                        jsonfilename_stem)
    scenario_path = process_dir.joinpath(jsonfilename_stem)
    shutil.copytree(param_files_path, scenario_path)
    # get other json files path
    other_json_files_path = run_dir.joinpath(sim_info_dir_name, "json_files")
    return(scenario_path, jsonfilename_stem, process_dir,
           other_json_files_path, number_of_replics, scenarionumber_zb)

def scenario_launch_laptop(jsonfilename_stem,
                           process_dir, other_json_files_path,
                           number_of_replics,
                           scenario_number_zb,
                           hyperparam_tuning=False,
                           tuning_path=""):
    """A docstring for the scenario_run_laptop function."""
    # jsonfilename_stem = jsonfilename.split(".")[0]
    scenario_path = process_dir.joinpath(jsonfilename_stem)
    scenario_datagen_params_path = scenario_path.joinpath("datagen_params")
    for replicnum in range(1, number_of_replics + 1):
        scenario_run_single_replic.main(jsonfilename_stem,
                                        other_json_files_path,
                                        scenario_path,
                                        scenario_datagen_params_path,
                                        scenario_number_zb,
                                        replicnum,
                                        number_of_replics,
                                        hyperparam_tuning,
                                        tuning_path)

def scenario_launch_cluster(jsonfilename_stem,
                            process_dir, other_json_files_path,
                            number_of_replics,
                            scenario_number_zb):
    # note that for environ cluster, must check that all runs are complete
    #     before generating final report
    scenario_path = process_dir.joinpath(jsonfilename_stem) 
    scenario_datagen_params_path = scenario_path.joinpath("datagen_params")
    for replicnum in range(1, number_of_replics + 1):
        fstring_part_one = f"JSONFILENAMESTEM={jsonfilename_stem}"
        fstring_part_one_pt_five = f"OTHERJSONFILESPATH={other_json_files_path}"
        fstring_part_two = f"SCENARIOPATH={scenario_path}"
        fstring_part_two_pt_five = f"SCENARIODATAGENPARAMSPATH={scenario_datagen_params_path}"
        fstring_part_three = f"SCENARIONUMBERZB={scenario_number_zb}"        
        fstring_part_four = f"REPLICNUM={replicnum}"
        fstring_part_five = f"NUMBEROFREPLICS={number_of_replics}"
        input_string = "--export" + "=" + fstring_part_one + "," \
            + fstring_part_one_pt_five + "," \
            + fstring_part_two + "," + fstring_part_two_pt_five + "," \
            + fstring_part_three \
            + "," + fstring_part_four + "," + fstring_part_five
        subprocess.run(["sbatch", input_string, "simulation.slurm"])

## main

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", choices = ["laptop", "cluster"],
                        required=True)
    parser.add_argument("--sim_info_dir_name", required=True)
    parser.add_argument("--scenarionumber", required=True)
    # parse args
    args = parser.parse_args()
    environ = args.environ
    sim_info_dir_name = args.sim_info_dir_name
    scenarionumber = args.scenarionumber
    # run scenario_launch_setup
    scenario_path, jsonfilename_stem, process_dir, \
        other_json_files_path, number_of_replics, \
        scenarionumber_zb = scenario_launch_setup(environ,
                                                  sim_info_dir_name,
                                                  scenarionumber)
    # run replications
    hyperparam_tuning = False
    tuning_path = process_dir.joinpath("hyperparam_tuning")
    if environ == "laptop":
        scenario_launch_laptop(jsonfilename_stem, process_dir,
                               other_json_files_path,
                               number_of_replics, scenarionumber_zb,
                               hyperparam_tuning,
                               tuning_path)
    elif environ == "cluster":
        scenario_launch_cluster(jsonfilename_stem, process_dir,
                                other_json_files_path,
                                number_of_replics, scenarionumber_zb)
