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

# imports

## python extension module
from probitlcm import _core

## other modules from me
from probitlcm import scenario_report_single_replic
from probitlcm import convergence_test

## standard library
from pathlib import Path
import argparse
import json

## package
import numpy as np

# functions
def main(jsonfilename_stem, other_json_files_path,
         scenario_path, scenario_datagen_params_path,
         scenario_number_zb, replicnum,
         number_of_replics,
         hyperparam_tuning=False,
         tuning_path=""):
    # create directory for this particular replication
    replicnum_dirname = f"replic_{replicnum:03}"
    replic_path = scenario_path.joinpath(replicnum_dirname)
    replic_path.mkdir()
    # run simulation
    _core.run_replication(jsonfilename_stem,
                          scenario_number_zb,
                          replicnum,
                          number_of_replics,
                          str(other_json_files_path),
                          str(scenario_path),
                          str(scenario_datagen_params_path),
                          str(replic_path),
                          hyperparam_tuning,
                          str(tuning_path))
    # build report
    jsonfilename = jsonfilename_stem + ".json"
    jsonfile_src_path  = scenario_path.joinpath(jsonfilename)
    data = json.loads(jsonfile_src_path.read_bytes())
    # load other json files
    jsonfile_src_path = other_json_files_path.joinpath("01_fixed_vals.json")
    data_more = json.loads(jsonfile_src_path.read_bytes())
    data.update(data_more)
    # generate report
    if hyperparam_tuning == False:
        scenario_report_single_replic.generate_report(
            scenario_path, scenario_datagen_params_path,
            replic_path, data, replicnum)
        convergence_test.run_geweke_test(replic_path,
                                         data)
        # delete draws_*.txt
        for p in replic_path.glob("draws_*.txt"):
            p.unlink()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--jsonfilename_stem")
    parser.add_argument("--other_json_files_path")
    parser.add_argument("--scenario_path")
    parser.add_argument("--scenario_datagen_params_path")
    parser.add_argument("--scenario_number_zb")
    parser.add_argument("--replicnum")
    parser.add_argument("--number_of_replics")
    args = parser.parse_args()
    jsonfilename_stem = args.jsonfilename_stem
    other_json_files_path = Path(args.other_json_files_path)
    scenario_path = Path(args.scenario_path)
    scenario_datagen_params_path = Path(args.scenario_datagen_params_path)
    scenario_number_zb = int(args.scenario_number_zb)
    replicnum = int(args.replicnum)
    number_of_replics = int(args.number_of_replics)
    main(jsonfilename_stem, other_json_files_path,
         scenario_path, scenario_datagen_params_path,
         scenario_number_zb, replicnum, number_of_replics)
