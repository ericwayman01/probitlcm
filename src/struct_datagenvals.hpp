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

#ifndef HEADERS_STRUCT_DATAGENVALS_HPP_
#define HEADERS_STRUCT_DATAGENVALS_HPP_

#include <armadillo>

struct DatagenVals {
    // data
    arma::mat Xmat;
    // params
    arma::umat Ymat;
    arma::mat Ymat_star;
    arma::mat beta;
    arma::umat delta;
    arma::field<arma::mat> kappa; // really a row vector
    arma::umat alpha;
    arma::mat alpha_star;
    arma::mat alpha_star_expa;
    arma::field<arma::mat> gamma; // really a field of row vectors
    arma::field<arma::mat> gamma_expa; // really a field of row vectors
    arma::mat lambda;
    arma::mat lambda_expa;
    arma::mat Rmat;
    arma::mat Sigma;
    double omega;
    // misc values
    double age_rescale_factor;
    // derived parameters
    arma::field<arma::mat> theta_j_mats;
};

#endif
