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

#ifndef HEADERS_INITIALIZATIONS_HPP_
#define HEADERS_INITIALIZATIONS_HPP_

#include <armadillo>
#include "run_mcmc.hpp"

#include "struct_mcmcdraws.hpp"
#include "struct_othervals.hpp"

void find_initial_alphas(std::string, std::string,
                         const int, const int, const int,
                         arma::uvec, const int);

void initialize_star_variable(arma::mat &,
                              const arma::field<arma::mat> &,
                              const arma::uvec &,
                              const arma::uword &,
                              const arma::umat &,
                              const arma::uword &, const arma::uword &);

void initialize_mcmc_variables(const OtherVals &,
                               MCMCDraws &,
                               const arma::umat &,
                               const arma::mat &,
                               int,
                               const DatagenVals &,
                               bool);

#endif
