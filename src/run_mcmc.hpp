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

#ifndef HEADERS_RUN_MCMC_HPP_
#define HEADERS_RUN_MCMC_HPP_

#include "struct_datagenvals.hpp"
#include "struct_mcmcdraws.hpp"
#include "struct_othervals.hpp"

#include <armadillo>
#include <nlohmann/json.hpp>

#include <string>

double log_sum_exp(double[], size_t);

arma::field<arma::vec> compute_V_quantities(const arma::mat&);

arma::uvec calc_vec_of_indices(const arma::uword &, const arma::uword &);

void prepare_thresholds(arma::field<arma::mat> &,
                        const arma::Col<arma::uword> &,
                        const arma::uword,
                        arma::uword,
                        double, double);

void initialize_alpha_and_alpha_star_cross_sec(const OtherVals &,
                                               const arma::field<arma::mat> &,
                                               const arma::mat &,
                                               const arma::mat &,
                                               const arma::mat &,
                                               arma::mat &,
                                               arma::mat &,
                                               arma::umat &);

void initialize_alpha_and_alpha_star_longit(const OtherVals &,
                                            const arma::field<arma::mat> &,
                                            const arma::mat &,
                                            const arma::mat &,
                                            const arma::mat &,
                                            const arma::mat &,
                                            arma::mat &,
                                            arma::mat &,
                                            arma::umat &);


void perform_transformations(const OtherVals &,
                             MCMCDraws &,
                             const arma::uword &);

// OtherVals set_up_othervals(const nlohmann::json &);

void set_up_othervals(OtherVals &,
                      const nlohmann::json &,
                      std::string, bool);

// void set_up_othervals_data_analysis(OtherVals &,
//                                     const nlohmann::json &);

// begin find_best_permutation function section

arma::uword factorial(arma::uword);

int indexofSmallestElement(double array[], arma::uword);

arma::umat load_umat(std::string, std::string);
arma::uvec get_best_perm(const arma::umat &,
                         arma::uword);

arma::uword find_index_num_of_best_permutation(arma::cube &, arma::mat &,
                                               OtherVals &);

arma::umat calc_mode_alphas(arma::umat &, OtherVals &);

arma::uword find_index_num_of_best_inverse_perm_cb(arma::umat &,
                                                   arma::umat &,
                                                   OtherVals &);

// end find_best_permutation function section

double calc_theta_jmc(arma::uword, const arma::Row<arma::uword> &,
                      const arma::vec &, const arma::rowvec &);

void calc_theta_j_quantities(arma::mat &,
                             const arma::vec &,
                             const arma::rowvec &,
                             const arma::uword M_j,
                             const OtherVals &);

arma::mat prepare_w(const arma::mat &, const OtherVals &,
                    const arma::Mat<arma::uword> &);


void update_class_counts(MCMCDraws &, const OtherVals &,
                         const arma::uvec &);

arma::urowvec calc_total_class_counts_for_single_draw(
            const OtherVals &,
            const arma::uvec &,
            const arma::uword &);

void run_mcmc(DatagenVals &, OtherVals &, MCMCDraws &,
              const arma::umat &,
              const arma::mat &,
              bool);

nlohmann::json read_json_file(std::string);

void run_replication(std::string,
                     int,
                     int,
                     int,
                     std::string,
                     std::string,
                     std::string,
                     std::string,
                     bool,
                     std::string);

void run_data_analysis_chain(int,
                             std::string,
                             std::string,
                             std::string,
                             int,
                             bool,
                             std::string);

#endif
