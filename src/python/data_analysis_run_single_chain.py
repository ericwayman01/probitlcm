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

# imports

## python extension module
from probitlcm import _core

## other modules from me
from probitlcm import data_analysis_report_single_chain
from probitlcm import convergence_test

## standard library
import argparse
from pathlib import Path
import json


def main(setup_num, dataset_dir, data_analysis_path, chainnum,
         hyperparam_tuning=False, tuning_path=""):
    # create directory for this particular chain
    chainnum_dirname = f"chain_{chainnum:03}"    
    chain_results_path = data_analysis_path.joinpath(chainnum_dirname)
    chain_results_path.mkdir()
    # run chain
    _core.run_data_analysis_chain(setup_num,
                                  str(dataset_dir),
                                  str(data_analysis_path),
                                  str(chain_results_path),
                                  chainnum,
                                  hyperparam_tuning,
                                  str(tuning_path))
    # build report
    jsonfilename = f"setup_{setup_num:04}.json"
    jsonfile_src_path  = data_analysis_path.joinpath(jsonfilename)
    data = json.loads(jsonfile_src_path.read_bytes())
    # load other json files
    jsonfile_src_path = dataset_dir.joinpath("01_fixed_vals.json")
    data_more = json.loads(jsonfile_src_path.read_bytes())
    data.update(data_more)
    # generate report
    if hyperparam_tuning == False:        
        data_analysis_report_single_chain.generate_report(
            dataset_dir, data_analysis_path, chain_results_path,
            data, chainnum)
        convergence_test.run_geweke_test(chain_results_path,
                                         data)
        # delete draws_*.txt
        # for p in chain_results_path.glob("draws_*.txt"):
        #     p.unlink()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--setup_num")
    parser.add_argument("--dataset_dir")
    parser.add_argument("--data_analysis_path")
    parser.add_argument("--chainnum")
    args = parser.parse_args()
    setup_num = int(args.setup_num)
    dataset_dir = Path(args.dataset_dir)
    data_analysis_path = Path(args.data_analysis_path)
    chainnum = int(args.chainnum)
    main(setup_num, dataset_dir, data_analysis_path, chainnum)
