/*
 * This file is part of "probitlcm" which is released under GPL v3.
 *
 * Copyright (c) 2024 Jesse Bowers <jessemb2@illinois.edu>.
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

#ifndef HEADERS_CROSSCOUNTS_HPP_
#define HEADERS_CROSSCOUNTS_HPP_

#include <armadillo>

arma::uvec getBasisVec(const arma::umat& x, const arma::uvec& nvalues, const arma::uvec& items);
arma::uvec aggregateCounts(const arma::uvec& xcol, arma::uword nvalues);
arma::uvec aggregateAllCounts(const arma::umat& x, const arma::uvec& nvalues);
int diffcounts(const arma::umat& x1, const arma::umat& x2, const arma::uvec& nvalues);
arma::uvec concatenate(const std::vector<arma::uvec>& vecs);

#endif
