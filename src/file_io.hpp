/*
 * This file is part of "probitlcm" which is released under GPL v3.
 *
 * Copyright (c) 2022-2024 Eric Alan Wayman <ewayman2@illinois.edu>.
 *
 * This program is FLO (free/libre/open) software: you can redistribute
 * it and/or modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HEADERS_FILE_IO_HPP_
#define HEADERS_FILE_IO_HPP_

#include "struct_datagenvals.hpp"
#include "struct_mcmcdraws.hpp"
#include "struct_othervals.hpp"

#include <armadillo>

#include <string>

// template functions

template <typename TName>
void write_output_cube(OtherVals & othervals,
                       arma::Cube<TName> & draws_cube,
                       std::string parameter_name) {
    std::string filename;
    filename = othervals.replic_path +
        "/" + "draws_" + parameter_name + ".txt";
    draws_cube.save(filename, arma::arma_ascii);
}


template <typename TName>
void write_output_cube_average(OtherVals & othervals,
                               arma::Cube<TName> & draws_cube,
                               std::string parameter_name,
                               arma::uword chain_length_after_burnin) {
    std::string output_directory("./output_results");
    std::string filename;

    arma::Cube<TName> results_to_save;
    results_to_save = draws_cube.tail_slices(chain_length_after_burnin);
    arma::mat mat_of_average(draws_cube.n_rows, draws_cube.n_cols);
    mat_of_average.zeros();
    // uword n_slices = draws_cube.n_slices;
    arma::mat current_slice(draws_cube.n_rows, draws_cube.n_cols);
    for (arma::uword n = 0; n < chain_length_after_burnin; ++n) {
        current_slice = arma::conv_to<arma::mat>::from(
            results_to_save.slice(n));
        // mat_of_average += results_to_save.slice(n);
        mat_of_average += current_slice;
    }
    mat_of_average /= chain_length_after_burnin;
    filename = othervals.replic_path + "/" + "average_"
        + parameter_name + ".txt";
    mat_of_average.save(filename, arma::arma_ascii);
}

// helper functions

std::string pad_string_with_zeros(int, std::string,
                                  std::string);

// major functions

void write_datagen_values(const OtherVals &,
                          const DatagenVals &);

void write_mcmc_output(MCMCDraws &, DatagenVals &,
                       OtherVals &, bool, bool, int);

#endif
