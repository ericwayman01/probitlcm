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
        
#include "struct_datagenvals.hpp"
#include "struct_mcmcdraws.hpp"
#include "struct_othervals.hpp"

#include "latent_state_related.hpp"
#include "log_of_norm_cdf.hpp"
#include "mvnorm.hpp"
#include "run_mcmc.hpp"
#include "sampling.hpp"
#include "sampling_generic.hpp"

#include <armadillo>

#include <vector>
#include <algorithm> // for std::min, std::max, std::max_element
#include <cmath> // for std::log, std::sqrt



// functions for sampling block 2 for this model

std::vector<double> calc_part_of_conds(
            const arma::Mat<arma::uword> & design_matrix,
            arma::vec beta_j, arma::uword h, arma::uword H) {
    std::vector<double> results;
    // for the outer loop
    arma::Row<arma::uword> d_u(H);
    arma::vec beta_j_less_h(H - 1);
    // for the inner loop
    arma::Row<arma::uword> d_v(H);
    double dotprod;
    beta_j_less_h = beta_j;
    beta_j_less_h(0) = 0;
    beta_j_less_h(h - 1) = 0;
    // 0-based
    arma::Row<arma::uword> d_u_copy;
    for (arma::uword i = 1; i < H; ++i) {
        d_u = design_matrix.row(i);
        d_u(0) = 0;
        d_u(h - 1) = 0;
        arma::uvec ids = arma::find(d_u == 1);
        for (auto const &j: ids) {
            d_v = design_matrix.row(j);
            d_v(0) = 0;
            if (d_v(h - 1) == 0) {
                dotprod = arma::dot(d_u - d_v, beta_j_less_h);
                results.push_back(-1.0 * dotprod);                    
            }
        }
    }
    return results;    
}

double calc_L(std::vector<double> & part_of_conds) {
    double overall_max = *std::max_element(part_of_conds.begin(),
                                           part_of_conds.end());
    return overall_max;
}

bool check_cond1(std::vector<double> & part_of_conds) {
    double overall_max = *std::max_element(part_of_conds.begin(),
                                           part_of_conds.end());
    if (overall_max <= 0) {
        return true;
    } else {
        return false;
    }
}

void sample_beta_and_delta(DatagenVals & datagenvals,
                           const OtherVals & othervals,
                           MCMCDraws & draws,
                           arma::uword draw_number) {
    const arma::uword & J = othervals.dimensions.at("J");
    const arma::uword & H = othervals.dimensions.at("H");
    double omega = draws.omega(draw_number - 1);
    arma::mat dMat = get_design_vectors_from_alpha(draws.alpha_previous,
                                                   othervals.design_matrix,
                                                   othervals.basis_vector);
    arma::mat dMat_prime_dMat = arma::trans(dMat) * dMat;
    arma::mat dMat_prime_y_star = arma::trans(dMat) * draws.Ymat_star_current;
    double sigma_beta_sq = othervals.fixed_constants.at("sigma_beta_sq");
    double sigma_beta = std::sqrt(sigma_beta_sq);
    // mat to write over, which we assign to draws struct at the end
    arma::Mat<arma::uword> delta_draw = draws.delta.slice(draw_number - 1);
    // enter loops
    // note: the vector we use for beta is really the previous block beta,
    //       whose elements are, one-by-one, either replaced or left unchanged.
    //       so we set our (new) beta_j to previous_beta_j, and work through
    //       element-by-element.
    // Note: Parallelized with OpenMP
    #pragma omp parallel for
    for (arma::uword j = 1; j <= J; ++j) {
        arma::vec beta_j = draws.beta.slice(draw_number - 1).col(j - 1);
        arma::Col<arma::uword> delta_j = draws.delta.slice(
            draw_number - 1).col(j - 1);
        for (arma::uword h = 1; h <= H; ++h) {
            arma::uword delta_hj_prime;
            double beta_hj_prime;
            std::vector<double> part_of_conds = calc_part_of_conds(
                othervals.design_matrix, beta_j, h, H);
            double L = calc_L(part_of_conds);
            if (h == 1) {
                L = -1 * arma::datum::inf;
            }
            // Calculate c_2 and then c_1
            double c_2_sq_denom = dMat_prime_dMat(h - 1,
                                                  h - 1) + (1 / sigma_beta_sq);
            double c_2_sq = 1 / c_2_sq_denom;
            double c_2 = std::sqrt(c_2_sq);
            arma::vec beta_j_tmp = beta_j;
            beta_j_tmp(h - 1) = 0;
            arma::vec c_1_tmp = dMat_prime_y_star.col(j - 1) -
                (dMat_prime_dMat * beta_j_tmp);
            double c_1 = c_2_sq * c_1_tmp(h - 1);
            double c_1_sq = c_1 * c_1;
            if (check_cond1(part_of_conds)) { // it is not a point mass
                double log_numer_part_1 = std::log(omega) -
                    std::log(arma::normcdf(-1.0 * L / sigma_beta, 0.0, 1.0)) +
                    std::log(c_2) - std::log(sigma_beta);                
                double log_numer_part_2 = (c_1_sq / (2 * c_2_sq)) +
                    log_of_norm_cdf(-(L - c_1) / c_2);
                double log_numer = log_numer_part_1 + log_numer_part_2;
                double s_array[2];
                // double log_denom = log_numer + std::log(1 - omega);
                s_array[0] = log_numer;
                s_array[1] = std::log(1 - omega);
                double log_denom = log_sum_exp(s_array, 2);
                double omega_tilde = std::exp(log_numer - log_denom);
                delta_hj_prime = sample_bernoulli_distn(omega_tilde);
            } else { // cond1 not satisfied, so it's a point-mass
                     // distribution at 1
                delta_hj_prime = 1;
            }
            if (h == 1) {
                delta_hj_prime = 1;
            }

            // delta_hj_prime = datagenvals.delta(h - 1, j - 1);
            // now that delta_hj has been sampled, sample beta_hj 
            if (delta_hj_prime == 0) {
                beta_hj_prime = 0;
            } else {
                if (h == 1) {
                    beta_hj_prime = arma::randn(arma::distr_param(c_1, c_2));
                    if (std::isnan(beta_hj_prime)) {
                        std::cout << "c_1 = "
                                  << std::to_string(c_1) << std::endl;
                        std::cout << "c_2 = "
                                  << std::to_string(c_2) << std::endl;
                        dMat_prime_y_star.print("dMat_prime_y_star:");
                    }
                } else {
                    beta_hj_prime = sample_normal_truncated(c_1, c_2, L,
                                                            arma::datum::inf);
                    if (std::isnan(beta_hj_prime)) {
                        std::cout << "c_1 = "
                                  << std::to_string(c_1) << std::endl;
                        std::cout << "c_2 = "
                                  << std::to_string(c_2) << std::endl;
                        std::cout << "L = "
                                  << std::to_string(L) << std::endl;
                        dMat_prime_y_star.print("dMat_prime_y_star:");
                    }
                }
            }
            // assign the sampled values
            beta_j.at(h - 1) = beta_hj_prime;
            // std::cout << "beta_hj_prime = " << beta_hj_prime << std::endl;
            delta_j.at(h - 1) = delta_hj_prime;
        }
        draws.beta.slice(draw_number).col(j - 1) = beta_j;
        draws.delta.slice(draw_number).col(j - 1) = delta_j;
    }
}

// helper functions for sampling block 3 for this model

SigmaRelatVals perform_sigma_calculations(const OtherVals & othervals,
                                          const arma::mat & Sigma) {
    const arma::uword & K = othervals.dimensions.at("K");
    SigmaRelatVals sigmarelatvals;
    sigmarelatvals.Sigma_neg_k_neg_k_s = arma::field<arma::mat>(K);
    sigmarelatvals.Sigma_k_neg_k_s = arma::field<arma::mat>(K);
    sigmarelatvals.Sigma_neg_k_k_s = arma::field<arma::mat>(K);
    sigmarelatvals.Sigma_k_k_s = arma::vec(K);
    // do sigma calculations
    // initialize Sigma temporary variables
    arma::mat Sigma_neg_k_neg_k(K - 1, K - 1);
    arma::mat Sigma_k_neg_k(1, K - 1);
    arma::mat Sigma_neg_k_k(K - 1, 1);
    double Sigma_k_k;
    // do the sigma_k calculations, one for each k
    arma::uvec vec_of_indices;
    for (arma::uword k = 1; k <= K; ++k) { // recall always K >= 2
        arma::uvec uvec_k = {k - 1};
        vec_of_indices = calc_vec_of_indices(k, K);
        Sigma_neg_k_neg_k = Sigma.submat(vec_of_indices,
                                         vec_of_indices);
        Sigma_k_neg_k = Sigma.submat(uvec_k, vec_of_indices);
        Sigma_neg_k_k = Sigma.submat(vec_of_indices, uvec_k);
        Sigma_k_k = Sigma(k - 1, k - 1);
        // assign results to components of fields
        sigmarelatvals.Sigma_neg_k_neg_k_s(k - 1) = Sigma_neg_k_neg_k;
        sigmarelatvals.Sigma_k_neg_k_s(k - 1) = Sigma_k_neg_k;
        sigmarelatvals.Sigma_neg_k_k_s(k - 1) = Sigma_neg_k_k;
        sigmarelatvals.Sigma_k_k_s(k - 1) = arma::as_scalar(Sigma_k_k);
    }
    return sigmarelatvals;
}

// cross-sec
arma::uword draw_alpha_nk(const OtherVals & othervals,
                          MCMCDraws & draws,
                          std::vector<double> cond_mean_and_sqrt,
                          const arma::urowvec & alpha_n,
                          double n,
                          double k,
                          const arma::uword draw_number) {
    const arma::uword & J = othervals.dimensions.at("J");
    arma::uword L_k = othervals.L_k_s(k - 1);
    arma::urowvec art_alpha_n = alpha_n;
    arma::vec prob_art_alpha_nk_values(L_k); // to store the values
    arma::vec log_multinoulli_cdf_values(L_k);
    double s_array[L_k];
    // std::vector<double> s_vector;
    // aliases
    const arma::mat & beta = draws.beta.slice(draw_number);
    const arma::rowvec & Ymat_star_n = draws.Ymat_star_current.row(n - 1);
    const arma::rowvec& gamma_expa_k = draws.gamma_expa(
        k - 1).row(draw_number - 1);
    // do work
    for (arma::uword i = 0; i <= L_k - 1; ++i) { // 0-based (inherently)
        art_alpha_n(k - 1) = i;
        arma::mat art_d_n = get_design_vectors_from_alpha(
            art_alpha_n, othervals.design_matrix, othervals.basis_vector);
        arma::vec art_d_n_beta_trans = arma::trans(art_d_n * beta);
        arma::vec y_n_star_trans = arma::trans(Ymat_star_n);        
        arma::vec single_log_norm_pdf = log_normpdf(
            y_n_star_trans, art_d_n_beta_trans, arma::ones(J));
        double log_g_negative_1 = arma::sum(single_log_norm_pdf);
        double log_cdf1 = log_of_norm_cdf((gamma_expa_k(i + 1) - 
                                           cond_mean_and_sqrt[0]) /
                                          cond_mean_and_sqrt[1]);
        double log_cdf2 = log_of_norm_cdf((gamma_expa_k(i) - 
                                           cond_mean_and_sqrt[0]) /
                                          cond_mean_and_sqrt[1]);

        double log_of_cdf_diff = calc_log_of_normal_cdf_difference(log_cdf1,
                                                                   log_cdf2);
        s_array[i] = log_g_negative_1 + log_of_cdf_diff;
    }
    double log_c = -1.0 * log_sum_exp(s_array, L_k);    
    for (arma::uword i = 0; i < L_k; ++i) {
        log_multinoulli_cdf_values(i) = log_c + log_sum_exp(s_array, i + 1);
    }
    double log_u = std::log(arma::randu());
    arma::uword i = 0;
    while (log_multinoulli_cdf_values(i) < log_u) {
        ++i;
    }
    // deal with the rare numerical issue that occurs even on the log scale
    //     (i.e. all cdf values are 0)
    if (i > L_k - 1) {
        i = L_k - 1;
    }
    arma::uword alpha_nk = i;
    return alpha_nk;
}

std::vector<double> calc_cond_mean_and_sqrt_crosssec(
            SigmaRelatVals sigmarelatvals,
            const arma::rowvec & data_n,
            const arma::mat & alpha_star_expa_draw,
            const arma::mat & alpha_star_slope_expa,
            const arma::uvec & inner_vec_of_indices,
            arma::uword n, arma::uword k, arma::uword N) {
    // initialize variable
    std::vector<double> cond_mean_and_sqrt;
    // do calculations
    arma::vec alpha_star_slope_expa_k = alpha_star_slope_expa.col(
        k - 1);
    arma::mat alpha_star_slope_expa_neg_k = alpha_star_slope_expa.cols(
        inner_vec_of_indices);
    arma::vec b = arma::inv(sigmarelatvals.Sigma_neg_k_neg_k_s(k - 1)) *
        sigmarelatvals.Sigma_neg_k_k_s(k - 1);
    double a_part_one = arma::as_scalar(data_n * alpha_star_slope_expa_k);
    double a_part_two = arma::as_scalar(data_n *
                                        alpha_star_slope_expa_neg_k * b);
    double a = a_part_one - a_part_two;
    double Sigma_kknegk_part_two = arma::as_scalar(
        sigmarelatvals.Sigma_k_neg_k_s(k - 1) *
        arma::inv(sigmarelatvals.Sigma_neg_k_neg_k_s(k - 1)) *
        sigmarelatvals.Sigma_neg_k_k_s(k - 1));
    double Sigma_kknegk = sigmarelatvals.Sigma_k_k_s(k - 1) -
        Sigma_kknegk_part_two;
    arma::rowvec myrow = alpha_star_expa_draw.row(n - 1);
    arma::rowvec myrow_two = myrow.cols(inner_vec_of_indices);
    double my_mean = a + arma::as_scalar(myrow_two * b);
    // store results in variable
    cond_mean_and_sqrt.push_back(my_mean);
    cond_mean_and_sqrt.push_back(sqrt(Sigma_kknegk));
    return cond_mean_and_sqrt;
}

void sample_alpha_and_alpha_star_expa_crosssec(
            const OtherVals & othervals, MCMCDraws & draws,
            const arma::mat & Xmat, const arma::uword & draw_number) {
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & K = othervals.dimensions.at("K");
    arma::mat Sigma = draws.Sigma.slice(draw_number - 1);
    arma::mat lambda_expa = draws.lambda_expa.slice(draw_number - 1);
    arma::mat alpha_star_expa_draw = draws.alpha_star_expa_previous;
    arma::Mat<arma::uword> alpha_draw = draws.alpha_previous;
    SigmaRelatVals sigmarelatvals = perform_sigma_calculations(
        othervals, Sigma);
    #pragma omp parallel for
    for (arma::uword n = 1; n <= N; ++n) {
        arma::rowvec X_n = Xmat.row(n - 1);
        for (arma::uword k = 1; k <= K; ++k) {
            const arma::rowvec & gamma_expa_k = draws.gamma_expa(
                k - 1).row(draw_number - 1);
            arma::uvec inner_vec_of_indices = calc_vec_of_indices(k, K);
            std::vector<double> cond_mean_and_sqrt = calc_cond_mean_and_sqrt_crosssec(
                    sigmarelatvals, X_n, alpha_star_expa_draw,
                    lambda_expa, inner_vec_of_indices,
                    n, k, N);
            arma::uword alpha_nk = draw_alpha_nk(
                othervals,
                draws,
                cond_mean_and_sqrt,
                alpha_draw.row(n - 1),
                n,
                k,
                draw_number);
            alpha_draw(n - 1, k - 1) = alpha_nk;
            // alpha_draw(n - 1, k - 1) = datagenvals.alpha(n - 1, k - 1);
            alpha_star_expa_draw(n - 1, k - 1) = sample_normal_truncated(
                cond_mean_and_sqrt[0], cond_mean_and_sqrt[1],
                gamma_expa_k(alpha_nk), gamma_expa_k(alpha_nk + 1));
        }
    }
    draws.alpha_current = alpha_draw;
    draws.alpha_star_expa_current = alpha_star_expa_draw;
}

// block 4 for this model
void sample_gamma_expa(const OtherVals & othervals,
                       MCMCDraws & draws,
                       const arma::uword & draw_number) {
    const arma::uword & K = othervals.dimensions.at("K");
    for (arma::uword k = 0; k < K; ++k) { // 0-based
        double gamma_expa_kl_lb, gamma_expa_kl_ub;
        double max_alpha_star_expas, min_alpha_star_expas;
        arma::uword L_k = othervals.L_k_s(k);
        arma::mat gamma_expa_k_draws_mat = draws.gamma_expa(k);
        arma::rowvec gamma_expa_k = gamma_expa_k_draws_mat.row(draw_number - 1);
        const arma::Col<arma::uword> & alpha_k = draws.alpha_current.col(k);
        const arma::vec & alpha_star_expa_k = draws.alpha_star_expa_current.col(
            k);
        arma::uvec ids;
        if (L_k > 2) {
            for (arma::uword l = 2; l < L_k; ++l) { // 0-based
                // calculate lower bound of truncated exponential
                ids = arma::find(alpha_k == l - 1);
                if (ids.n_elem == 0) {
                    gamma_expa_kl_lb = gamma_expa_k(l - 1);
                } else {
                    max_alpha_star_expas = arma::max(
                        alpha_star_expa_k.elem(ids));
                    gamma_expa_kl_lb = std::max(max_alpha_star_expas,
                                                gamma_expa_k(l - 1));
                }
                // calculate upper bound of truncated exponential
                ids = arma::find(alpha_k == l);
                if (ids.n_elem == 0) {
                    gamma_expa_kl_ub = gamma_expa_k(l+1);
                } else {
                    min_alpha_star_expas = arma::min(
                        alpha_star_expa_k.elem(ids));
                    gamma_expa_kl_ub = std::min(min_alpha_star_expas,
                                                gamma_expa_k(l+1));
                }
                // now sample gamma_{kl}
                if (l < L_k - 1) {
                    gamma_expa_k(l) = arma::randu(
                        arma::distr_param(gamma_expa_kl_lb, gamma_expa_kl_ub));
                } else {
                    gamma_expa_k(l) = sample_exponential_truncated(
                        0.001, gamma_expa_kl_lb, gamma_expa_kl_ub);
                }
            }
            draws.gamma_expa(k).row(draw_number) = gamma_expa_k;
        } else {
            draws.gamma_expa(k).row(draw_number) = draws.gamma_expa(k).row(
                draw_number - 1); // will just be (-inf, constant, inf)
        }
    }
}

// block 5
void sample_Sigma_and_lambda_expa(const OtherVals & othervals,
                                  MCMCDraws & draws,
                                  const arma::mat & Xmat,
                                  const arma::uword & draw_number) {
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & K = othervals.dimensions.at("K");
    const arma::uword & D = othervals.dimensions.at("D");
    const arma::mat & alpha_star_expa = draws.alpha_star_expa_current;
    // Sample Sigma
    arma::mat I_D(D, D);
    I_D.eye();
    arma::mat X_trans_X = arma::trans(Xmat) * Xmat;
    arma::mat X_trans_X_plus_I_D_inv = arma::inv(X_trans_X + I_D);
    arma::mat L_2_hat = X_trans_X_plus_I_D_inv *
        arma::trans(Xmat) * alpha_star_expa;
    arma::mat S_part_1 = alpha_star_expa - Xmat * L_2_hat;
    arma::mat S = arma::trans(S_part_1) * S_part_1 +
        arma::trans(L_2_hat) * L_2_hat;
    arma::mat I_K(K, K);
    I_K.eye();
    arma::mat my_matrix = arma::symmatu(S + I_K);
    arma::mat Sigma = arma::iwishrnd(my_matrix, N + K + 1);
    // set Sigma result
    draws.Sigma.slice(draw_number) = Sigma;
    // Sample lambda_expa
    arma::mat lambda_expa = sample_matrix_variate_normal(
        L_2_hat, X_trans_X_plus_I_D_inv, Sigma, D, K);
    draws.lambda_expa.slice(draw_number) = lambda_expa;
}

// block 6 for this model
void sample_omega(const OtherVals & othervals, MCMCDraws & draws,
                  const arma::uword & draw_number) {
    const arma::uword & H = othervals.dimensions.at("H");
    const arma::uword & J = othervals.dimensions.at("J");
    const double & omega_0 = othervals.fixed_constants.at("omega_0");
    const double & omega_1 = othervals.fixed_constants.at("omega_1");
    // do real work
    double a, b;
    double sum_delta = arma::accu(draws.delta.slice(draw_number));
    a = sum_delta + omega_0;
    b = H * J - sum_delta + omega_1;
    draws.omega(draw_number) = sample_beta_distn(a, b);
}
