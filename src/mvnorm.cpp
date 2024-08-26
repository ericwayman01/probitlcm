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

#include <armadillo>

#include <cmath>

/*
 * This function "log_det" is essentially the same as the log_det function from
 * "statslib" by Keith O'Hara; the relevant path in that software is:
 * include/stats_incl/misc/matrix_ops/log_det.hpp
 * That file was distributed under the Apache License, Version 2.0.
 * The copyrights for that release are as follows:
 * Copyright (C) 2011-2023 Keith O'Hara
 */
double log_det(const arma::mat & X) {
    arma::mat tmpmat = arma::chol(X, "lower");
    arma::vec my_vec = tmpmat.diag();
    return 2.0 * arma::accu(arma::log(my_vec));
}

/*
 * This function "log_mvnorm_pdf" was based on "dmvnorm" from
 * "statslib" by Keith O'Hara; the relevant path in that software is:
 * include/stats_incl/dens/dmvnorm.ipp
 * That file was distributed under the Apache License, Version 2.0.
 * The copyrights for that release are as follows:
 * Copyright (C) 2011-2023 Keith O'Hara
 */
double log_mvnorm_pdf(const arma::rowvec & X,
                      const arma::rowvec & mu_par,
                      const arma::mat & Rmat) {
    arma::uword K = Rmat.n_rows;
    const double cons_term = -0.5 * K * std::log(2.0 * arma::datum::pi);
    arma::rowvec X_term = X - mu_par;
    arma::vec X_term_trans = arma::trans(X_term);
    // arma::solve(Rmat, X_term) gives Rmat^{-1} X_term
    double quad_term = arma::as_scalar(
        arma::dot(X_term_trans, arma::solve(Rmat, X_term_trans)));
    double return_value = cons_term - 0.5 * (log_det(Rmat) + quad_term);
    return return_value;
}
