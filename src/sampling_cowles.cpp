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

#include "struct_mcmcdraws.hpp"
#include "struct_othervals.hpp"

#include "latent_state_related.hpp"
#include "sampling_generic.hpp"

#include <armadillo>

#include <cmath>


// these are block 1 for this model

double calculate_block1_i_part(const OtherVals & othervals,
                               const MCMCDraws& draws,
                               const arma::mat & dMat_beta,
                               const arma::Mat<arma::uword> & Ymat,
                               const arma::rowvec & kappa_j_prop,
                               const arma::rowvec & kappa_j_previous,
                               const arma::uword j) {
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    double lnipart = 0.0;
    double numer_left, numer_right, denom_left, denom_right;
    int y_n_j;
    double d_n_beta_j;
    for (arma::uword n = 1; n <= N; ++n) { // 1-based
        // 1-based
        d_n_beta_j = dMat_beta(n - 1, j - 1);
        // 1-based
        y_n_j = Ymat(n - 1, j - 1); // 1-based
        numer_left = arma::normcdf(kappa_j_prop(y_n_j + 1),
                                   d_n_beta_j, 1.0);
        numer_right = arma::normcdf(kappa_j_prop(y_n_j),
                                    d_n_beta_j, 1.0);
        denom_left = arma::normcdf(kappa_j_previous(y_n_j + 1),
                                   d_n_beta_j, 1.0);
        denom_right = arma::normcdf(kappa_j_previous(y_n_j),
                                    d_n_beta_j, 1.0);
        lnipart += std::log(numer_left - numer_right) -
            std::log(denom_left - denom_right);
    }
    return lnipart;    
}

void sample_ystar_and_kappa(OtherVals & othervals,
                            MCMCDraws & draws,
                            const arma::Mat<arma::uword> & Ymat,
                            const arma::uword & draw_number) {
    // get necessary parameter references
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & J = othervals.dimensions.at("J");
    const double & sigma_kappa_sq =
        othervals.fixed_constants.at("sigma_kappa_sq");
    double sigma_kappa = sqrt(sigma_kappa_sq);
    // perform necessary pre-calculations
    arma::mat dMat = get_design_vectors_from_alpha(
        draws.alpha_previous, othervals.design_matrix, othervals.basis_vector);
    arma::mat dMat_beta = dMat * draws.beta.slice(draw_number - 1);
    for (arma::uword j = 1; j <= J; ++j) { // 1-based
        const arma::rowvec & kappa_j_previous = draws.kappa(j - 1).row(
            draw_number - 1);
        arma::uword M_j = othervals.M_j_s(j - 1);
        arma::rowvec kappa_j_prop(M_j + 1);
        kappa_j_prop(0) = -1.0 * arma::datum::inf;
        kappa_j_prop(1) = 0.0;
        kappa_j_prop(M_j) = arma::datum::inf;
        if (M_j > 2) {
            // looping over the m's, this does two things: 
            //     (1) draw the m'th item of the proposed kappa_j block, and
            //     (2) calculate the "m part" of R_k^(r) for this particular m
            // for (1), this uses the method from burkardt2014normal,
            //          page 24
            for (arma::uword m = 2; m <= M_j - 1; ++m) { // 0-based inherently
                kappa_j_prop(m) = sample_normal_truncated(
                    kappa_j_previous(m), sigma_kappa,
                    kappa_j_prop(m - 1), kappa_j_previous(m + 1));
                if (kappa_j_prop(m) <= kappa_j_prop(m - 1)) {
                    kappa_j_prop(m) = kappa_j_prop(m);
                }
            }
            double lnmpart = 0.0;
            double numer_left, numer_right, denom_left, denom_right;
            for (arma::uword m = 2; m <= M_j - 1; ++m) {
                numer_left = arma::normcdf(
                    kappa_j_previous(m + 1), kappa_j_previous(m),
                    sigma_kappa);
                numer_right = arma::normcdf(
                    kappa_j_prop(m - 1), kappa_j_previous(m), sigma_kappa);
                denom_left = arma::normcdf(
                    kappa_j_prop(m + 1), kappa_j_prop(m), sigma_kappa);
                denom_right = arma::normcdf(
                    kappa_j_previous(m - 1), kappa_j_prop(m), sigma_kappa);
                lnmpart += std::log(numer_left - numer_right) -
                    std::log(denom_left - denom_right);
            }
            double log_i_part = calculate_block1_i_part(othervals,
                                                        draws, 
                                                        dMat_beta,
                                                        Ymat,
                                                        kappa_j_prop,
                                                        kappa_j_previous,
                                                        j);
            double log_u = std::log(arma::randu());
            if (log_u <= std::min(0.0, lnmpart + log_i_part)) {
                // accept the proposed values
                if (draw_number >= othervals.burnin) {
                    othervals.ystar_and_kappa_accept_count += 1;
                }
                // kappa_j_prop.print("kappa_j_prop:");
                draws.kappa(j - 1).row(draw_number) = kappa_j_prop;
                // prepare to calculate y_j_star
                double d_n_beta_j;
                arma::uword y_n_j;
                // calculate y_j_star
                #pragma omp parallel for
                for (arma::uword n = 1; n <= N; ++n) { // 1-based
                    d_n_beta_j = dMat_beta(n - 1, j -1);
                    y_n_j = Ymat(n - 1, j - 1); // 1-based
                    // 1-based
                    draws.Ymat_star_current(
                        n - 1,
                        j - 1) = sample_normal_truncated(
                            d_n_beta_j, 1.0,
                            kappa_j_prop(y_n_j),
                            kappa_j_prop(y_n_j + 1));
                }
            } else { // proposal rejected
                draws.kappa(j - 1).row(draw_number) = kappa_j_previous;
                draws.Ymat_star_current.col(j - 1) =
                    draws.Ymat_star_previous.col(j - 1);
            }
        } else { // for M_j = 2, we always accept. nothing to do for kappas.
                 // sample only Y-stars
            // note: this has nothing to do with the hyperparameter
            //       sigma_kappa_j
            // nothing to do for kappas
            draws.kappa(j - 1).row(draw_number) = kappa_j_previous;
            arma::rowvec kappa_j = draws.kappa(j - 1).row(draw_number);
            // sample y_j_star
            double d_n_beta_j;
            arma::uword y_n_j;
            // calculate y_j_star
            #pragma omp parallel for
            for (arma::uword n = 1; n <= N; ++n) { // 1-based
                d_n_beta_j = dMat_beta(n - 1, j -1);
                y_n_j = Ymat(n - 1, j - 1); // 1-based
                // 1-based
                draws.Ymat_star_current(
                    n - 1,
                    j - 1) = sample_normal_truncated(
                        d_n_beta_j, 1.0,
                        kappa_j(y_n_j), kappa_j(y_n_j + 1));
            }
        }
    }
}
