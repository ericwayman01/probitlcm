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
import pathlib

## other modules
import sys
if sys.version_info[1] < 11:
    import toml
else:
    import tomllib as toml

## my modules
from probitlcm import scenario_report_overall


def check_if_all_replics_finished(scenario_path, number_of_replics):
    # build list_of_replics to check
    list_of_replics = list()
    for replicnum in range(1, number_of_replics + 1):
        replicnum_dirname = f"replic_{replicnum:03}"
        replic_path = scenario_path.joinpath(replicnum_dirname,
                                                     "done.txt")
        list_of_replics.append(replic_path)
    completed_replics = 0
    for done_file in list_of_replics:
        if done_file.exists():
            completed_replics += 1
    if completed_replics == number_of_replics:
        return True
    else:
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", choices = ["laptop", "cluster"],
                        required=True)
    parser.add_argument("--sim_info_dir_name", required=True)
    parser.add_argument("--scenarionumber", required=True)
    args = parser.parse_args()

    # run_dir = pathlib.Path.cwd()
    # load config file
    with open('config_simulation.toml') as fileObj:
        config = toml.load(fileObj)
    if args.environ == "laptop":
        process_dir = config['laptop_process_dir']
    elif args.environ == "cluster":
        process_dir = config['cluster_process_dir']
    process_dir = pathlib.Path(process_dir)
    sim_info_dir_name = pathlib.Path(args.sim_info_dir_name)
    process_dir = process_dir.joinpath(sim_info_dir_name)

    scenario_number = int(args.scenarionumber)
    jsonfilename_stem = f"scenario_{scenario_number:04}"
    jsonfilename = jsonfilename_stem + ".json"
    scenario_path = process_dir.joinpath(jsonfilename_stem)
    schema_file_path = scenario_path.joinpath(jsonfilename)
    # find other_json_files_path
    run_dir = pathlib.Path.cwd()
    other_json_files_path = run_dir.joinpath(sim_info_dir_name, "json_files")
    # zero-based scenario number for generating seeds

    number_of_replics = config['number_of_replics']

    # generate report for scenario
    all_replics_finished = check_if_all_replics_finished(
        scenario_path, number_of_replics)
    if all_replics_finished is False:
        print("Not all replics are finished, so nothing done.")
    else:
        scenario_report_overall.generate_report(schema_file_path,
                                                other_json_files_path,
                                                scenario_path,
                                                number_of_replics)
        print("Generated final statistics and report.")
