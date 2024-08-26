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

import argparse
import pathlib

import sys
if sys.version_info[1] < 11:
    import toml
else:
    import tomllib as toml

def find_convergence_failures(environ, sim_info_dir_name,
                              total_num_of_scenarios):
    with open('config_simulation.toml') as fileObj:
        config = toml.load(fileObj)
    if args.environ == "laptop":
        process_dir = config['laptop_process_dir']
    elif args.environ == "cluster":
        process_dir = config['cluster_process_dir']
    process_dir = pathlib.Path(process_dir)
    number_of_replics = config['number_of_replics']
    # process the results
    results_string = ""
    for scenario_number in range(1, total_num_of_scenarios + 1):
        for replicnum in range(1, number_of_replics + 1):
            scenario_dirname = f"scenario_{scenario_number:04}"
            replicnum_dirname = f"replic_{replicnum:03}"
            fname = "convergence_test_results.txt"
            fpath = process_dir.joinpath(sim_info_dir_name,
                                         scenario_dirname,
                                         replicnum_dirname,
                                         fname)
            with fpath.open() as f:
                if "failed" in f.read():
                    fpath_str = str(fpath)
                    str_to_add = f"failure in: {fpath_str}"
                    results_string += str_to_add
                    results_string += "\n"
        print("finished scenario number ", str(scenario_number))
    fpath = process_dir.joinpath(sim_info_dir_name,
                                 "convergence_test_failures.txt")
    with open(fpath, 'w') as f:
        f.writelines(results_string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--environ", choices = ["laptop", "cluster"],
                        required=True)
    parser.add_argument("--sim_info_dir_name", required=True)
    parser.add_argument("--total_num_of_scenarios", required=True)
    args = parser.parse_args()
    find_convergence_failures(args.environ, args.sim_info_dir_name,
          int(args.total_num_of_scenarios))
