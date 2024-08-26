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

#ifndef HEADERS_LATENT_STATE_RELATED_HPP_
#define HEADERS_LATENT_STATE_RELATED_HPP_

#include <armadillo>

#include <string>

arma::urowvec generate_d_b_k_alpha(arma::urowvec &,
                                   arma::uword &,
                                   arma::uword &);

arma::umat generate_design_matrix(arma::umat,
                                  arma::uword,
                                  arma::uvec,
                                  arma::uword,
                                  arma::uvec);

arma::Col<arma::uword> calculate_basis_vector(arma::Col<arma::uword>,
                                              arma::uword);

arma::uvec convert_alpha_to_class_numbers(
    const arma::umat &, const arma::uvec &);

arma::mat get_design_vectors_from_alpha(
                const arma::Mat<arma::uword>&,
                const arma::Mat<arma::uword>&,
                const arma::Col<arma::uword>&);

arma::urowvec convert_position_number_to_alpha(arma::uword,
                                               arma::uword,
                                               arma::uvec);

#endif
