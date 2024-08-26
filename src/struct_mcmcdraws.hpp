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

#ifndef HEADERS_STRUCT_MCMCDRAWS_HPP_
#define HEADERS_STRUCT_MCMCDRAWS_HPP_

#include <armadillo>

struct MCMCDraws {
    // params
    arma::ucube Ymat_pred_chunk; // only in use for data analysis runs
    arma::mat Ymat_star_current;
    arma::mat Ymat_star_previous;
    arma::cube beta;
    arma::Cube<arma::uword> delta;
    arma::field<arma::mat> kappa;
    arma::Mat<arma::uword> alpha_current;
    arma::Mat<arma::uword> alpha_previous;
    arma::mat alpha_star_current;
    arma::mat alpha_star_previous;
    arma::cube Rmat;
    arma::cube lambda;
    arma::cube lambda_expa;
    arma::vec omega;
    // expanded params
    arma::field<arma::mat> gamma;
    arma::field<arma::mat> gamma_expa;
    arma::mat alpha_star_expa_current;
    arma::mat alpha_star_expa_previous;
    arma::cube Sigma;
    // functions of derived parameters
    arma::field<arma::mat> theta_j_mats_sums;
    arma::umat class_counts;
    arma::umat total_class_counts_per_draw;
};

#endif
