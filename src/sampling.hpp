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

#ifndef HEADERS_SAMPLING_HPP_
#define HEADERS_SAMPLING_HPP_

#include "struct_datagenvals.hpp"
#include "struct_mcmcdraws.hpp"
#include "struct_othervals.hpp"

#include <armadillo>

#include <vector>

// Sigma-related

struct SigmaRelatVals {
    arma::field<arma::mat> Sigma_neg_k_neg_k_s;
    arma::field<arma::mat> Sigma_k_neg_k_s;
    arma::field<arma::mat> Sigma_neg_k_k_s;
    arma::vec Sigma_k_k_s;
};

// functions

void sample_ystar_and_kappa(const OtherVals &,
                            MCMCDraws &,
                            const arma::uword &);

std::vector<double> calc_part_of_conds(
            const arma::Mat<arma::uword> &,
            arma::vec, arma::uword, arma::uword);

double calc_L(std::vector<double> &);

bool check_cond1(std::vector<double> &);

void sample_beta_and_delta(DatagenVals & datagenvals,
                           const OtherVals &,
                           MCMCDraws &,
                           arma::uword);

std::vector<double> calc_cond_mean_and_sqrt_crosssec(
            SigmaRelatVals,
            const arma::rowvec &,
            const arma::mat &,
            const arma::mat &,
            const arma::uvec &,
            arma::uword, arma::uword, arma::uword);

std::vector<double> calc_cond_mean_and_sqrt_longit(
            SigmaRelatVals,
            const arma::rowvec &,
            const arma::rowvec &,
            const arma::rowvec &,
            const arma::mat &,
            const arma::mat &,
            const arma::uvec &,
            arma::uword);

arma::uword draw_alpha_nk(const OtherVals &,
                          MCMCDraws &,
                          std::vector<double>,
                          const arma::Row<arma::uword> &,
                          double,
                          double,
                          const arma::uword);

arma::uword draw_alpha_nk_t(const OtherVals &,
                            MCMCDraws &,
                            std::vector<double>,
                            const arma::mat &,
                            const arma::umat &,
                            const arma::mat &,
                            const arma::mat &,
                            const arma::mat &,
                            arma::uword, arma::uword,
                            arma::uword, arma::uword, arma::uword,
                            const arma::uword);

void sample_alpha_and_alpha_star_expa_crosssec(const OtherVals &,
                                               MCMCDraws &,
                                               const arma::mat &,
                                               const arma::uword &);

void sample_alpha_and_alpha_star_expa_longit(
    DatagenVals & datagenvals,
            const OtherVals &, MCMCDraws &,
            const arma::mat &, const arma::uword &);

void sample_gamma_expa(const OtherVals &,
                       MCMCDraws &,
                       const arma::uword &);

void sample_Sigma_and_lambda_expa(const OtherVals &, MCMCDraws&,
                                          const arma::mat &,
                                          const arma::uword &);

void sample_Sigma_and_zeta_expa(const OtherVals &,
                                MCMCDraws &,
                                const arma::mat &,
                                const arma::mat &,
                                const arma::uword &);

void sample_omega(const OtherVals &, MCMCDraws &,
                  const arma::uword &);

#endif
