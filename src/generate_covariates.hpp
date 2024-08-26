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

#ifndef HEADERS_COVARIATES_HPP_
#define HEADERS_COVARIATES_HPP_

#include <armadillo>

arma::vec generate_age_data(arma::uword, arma::uword, bool, double &);
arma::vec generate_sex_data(arma::uword, arma::uword);
arma::mat generate_covariates(arma::uword, arma::uword, double &);
arma::mat generate_covariates_scenario(arma::uword, arma::uword, int);

#endif
