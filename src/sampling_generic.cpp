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

#include "log_of_norm_cdf.hpp"
#include "sampling_generic.hpp"
#include "qnorm.hpp"

#include <armadillo>

#include <cmath> // for std::log1p

// math helpers

// https://github.com/scipy/scipy/issues/13923
double calc_log_of_normal_cdf_difference(double log_cdf2, double log_cdf1) {
    // cdf2 > cdf1
    // since cdf2 > cdf1,
    // log(cdf2 - cdf1) = log(cdf2) + log(1 - exp(log(cdf1) - log(cdf2)))
    // So log(cdf2 - cdf1) = log(cdf2) + std::log1p(-1.0 * log_cdf2 - log_cdf1);
    // we return log(cdf2 - cdf1)
    return log_cdf2 + std::log1p(-1.0 * std::exp(log_cdf1 - log_cdf2));
}

// from https://stackoverflow.com/questions/45943750/calculating-log-sum-exp-function-in-c
// formula: https://nhigham.com/2021/01/05/what-is-the-log-sum-exp-function/
double log_sum_exp(double arr[], size_t count) {
   if (count > 0) {
      double maxVal = arr[0];
      double sum = 0;
      for (unsigned int i = 1; i < count; i++) {
         if (arr[i] > maxVal){
            maxVal = arr[i];
         }
      }
      for (unsigned int i = 0; i < count; i++) {
          sum += std::exp(arr[i] - maxVal);
      }
      return std::log(sum) + maxVal;
   }
   else {
      return 0.0;
   }
}


// https://github.com/scipy/scipy/issues/13923
double sample_normal_truncated(const double mu_bar, const double sigma_bar,
                               const double a, const double b) {
    double u =  arma::randu();
    double q2 = arma::normcdf(b, mu_bar, sigma_bar);
    double q1 = arma::normcdf(a, mu_bar, sigma_bar);
    double q = q1 + u * (q2 - q1);
    
    if (q <= 0.0) {
        q = 0.0000000000001;
    } else if (q >= 1.0) {
        q = 0.9999999999999;
    }
    double quant = qnorm(q, 0, 1);
    return mu_bar + sigma_bar * quant;
}

double sample_exponential_truncated(const double rate_param,
                                    const double lower_val,
                                    const double upper_val) {
    double u =  arma::randu();
    double part1 = u * (std::exp(-rate_param * lower_val) -
                        std::exp(-rate_param * upper_val));
    return (-1.0 / rate_param) * std::log(
        std::exp(-rate_param * lower_val) - part1);
}

// M is p by n, Sigma is p by p, Psi is n by n
arma::mat sample_matrix_variate_normal(arma::mat M,
                                       arma::mat Sigma, arma::mat Psi,
                                       arma::uword p, arma::uword n) {
    arma::vec vec_of_M_prime = arma::reshape(arma::trans(M), n * p, 1);
    arma::mat Sigma_otimes_Psi = arma::kron(Sigma, Psi);
    arma::vec vec_of_X_prime = arma::mvnrnd(vec_of_M_prime, Sigma_otimes_Psi);
    arma::mat X_trans = reshape(vec_of_X_prime, n, p);
    arma::mat X = trans(X_trans);
    return X;
}

arma::mat sample_matrix_variate_normal_indep_rows(arma::mat M,
                                                  arma::mat Psi,
                                                  arma::uword p,
                                                  arma::uword n) {
    arma::mat X(p, n);
    arma::vec X_row_prime;
    for (arma::uword i = 0; i < p; ++i) {
        X_row_prime = arma::mvnrnd(trans(M.row(i)), Psi);
        X.row(i) = arma::trans(X_row_prime);
    }
    return X;
}

double sample_beta_distn(const double shape, const double scale) {    
    double x = arma::randg(arma::distr_param(shape, 1.0));
    double y = arma::randg(arma::distr_param(scale, 1.0));
    return x / (x + y);
}

arma::uword sample_bernoulli_distn(const double p) {
    double u = arma::randu();
    if (u < p) {
        return 1;
    } else {
        return 0;
    }
}

// more specific sampling situations

// ordinal "sampling" of Ymat from latent variable values Zmat
void sample_ordinal_values(arma::umat & Ymat,
                           const arma::mat & Zmat,
                           const arma::field<arma::mat> & gamma,
                           const arma::Col<arma::uword> & L_k_s,
                           const arma::uword N,
                           const arma::uword T,
                           const arma::uword K) {
    for (arma::uword k = 1; k <= K; ++k) {
        arma::uword L_k = L_k_s(k - 1);
        const arma::rowvec & gamma_k = gamma(k - 1).row(0);
        arma::Col<arma::uword> Ymat_k(N);
        for (arma::uword l = 0; l <= L_k - 1; ++l) {
            // elems1
            arma::uvec elems1_full(N);
            elems1_full.zeros();
            arma::uvec elems1 = arma::find(
                Zmat.col(k - 1) > gamma_k(l));
            elems1_full.elem(elems1).fill(1);
            // elems2
            arma::uvec elems2_full(N);
            elems2_full.zeros();
            arma::uvec elems2 = arma::find(
                Zmat.col(k - 1) <= gamma_k(l + 1));
            elems2_full.elem(elems2).fill(1);
            arma::uvec elems_full = elems1_full % elems2_full;
            arma::uvec elems = arma::find(elems_full == 1);
            Ymat_k.elem(elems).fill(l);
        }
        Ymat.col(k - 1) = Ymat_k;
    }
}

// ordinal "sampling" of Ymat from latent variable values Zmat
void sample_ordinal_values_newer(arma::umat & Ymat,
                                 const arma::mat & Zmat,
                                 const arma::field<arma::mat> & gamma,
                                 const arma::uvec & L_k_s,
                                 const arma::uword K,
                                 const arma::uword draw_number) {
    arma::uword nrows = Ymat.n_rows;
    for (arma::uword k = 1; k <= K; ++k) {
        arma::uword L_k = L_k_s(k - 1);
        const arma::rowvec & gamma_k = gamma(k - 1).row(draw_number);
        arma::uvec Ymat_k(nrows);
        for (arma::uword l = 0; l <= L_k - 1; ++l) {
            // elems1
            arma::uvec elems1_full(nrows);
            elems1_full.zeros();
            arma::uvec elems1 = arma::find(
                Zmat.col(k - 1) > gamma_k(l));
            elems1_full.elem(elems1).fill(1);
            // elems2
            arma::uvec elems2_full(nrows);
            elems2_full.zeros();
            arma::uvec elems2 = arma::find(
                Zmat.col(k - 1) <= gamma_k(l + 1));
            elems2_full.elem(elems2).fill(1);
            arma::uvec elems_full = elems1_full % elems2_full;
            arma::uvec elems = arma::find(elems_full == 1);
            Ymat_k.elem(elems).fill(l);
        }
        Ymat.col(k - 1) = Ymat_k;
    }
}
