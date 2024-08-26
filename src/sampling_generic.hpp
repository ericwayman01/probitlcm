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

#ifndef HEADERS_SAMPLING_GENERIC_HPP_
#define HEADERS_SAMPLING_GENERIC_HPP_

#include <armadillo>

double calc_log_of_normal_cdf_difference(double, double);

double log_sum_exp(double arr[], size_t count);

double inverse_of_normal_cdf(const double, const double, const double);

double sample_normal_truncated(const double, const double,
                               const double, const double);

double sample_exponential_truncated(const double,
                                    const double,
                                    const double);

arma::mat sample_matrix_variate_normal(arma::mat,
                                       arma::mat, arma::mat,
                                       arma::uword, arma::uword);

arma::mat sample_matrix_variate_normal_indep_rows(arma::mat,
                                                  arma::mat,
                                                  arma::uword,
                                                  arma::uword);

double sample_beta_distn(const double, const double);

arma::uword sample_bernoulli_distn(const double);

void sample_ordinal_values(arma::Mat<arma::uword> &,
                           const arma::mat &,
                           const arma::field<arma::mat> &,
                           const arma::Col<arma::uword> &,
                           const arma::uword,
                           const arma::uword,
                           const arma::uword);

void sample_ordinal_values_newer(arma::umat &,
                                 const arma::mat &,
                                 const arma::field<arma::mat> &,
                                 const arma::uvec &,
                                 const arma::uword,
                                 const arma::uword);

#endif
