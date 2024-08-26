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

#include "sampling_generic.hpp"

#include <armadillo>

arma::vec generate_age_data(arma::uword T, arma::uword N, bool rescale,
                            double & rescale_factor) {
    // Ages 25-49.999: 55%
    // Ages 50-59.999: 20%
    // Ages >=60: 25%

    arma::Col<arma::uword> age_bins(N);
    arma::vec age_pmf_values = {0.55, 0.20, 0.25};
    arma::vec age_cdf_values = arma::cumsum(age_pmf_values);
    double u;
    for (arma::uword n = 0; n < N; ++n) {
        u = arma::randu();
        arma::uword i = 0;
        while (age_cdf_values(i) <= u) {
            ++i;
        }
        age_bins(n) = i;
    }
    arma::Col<arma::uword> full_age_vec(N);
    arma::Col<arma::uword> initial_age_vec(N);
    double mu_bar, sigma_sq_bar, a, b;
    sigma_sq_bar = 3.0;
    for (arma::uword n = 0; n < N; ++n) {
        if (age_bins(n) == 0) {
            mu_bar = 37.5; a = 25.0; b = 50.0;
        } else if (age_bins(n) == 1) {
            mu_bar = 55.0; a = 50.0; b = 60.0;
        } else { //age_bins(n) == 2
            mu_bar = 65.0; a = 60.0; b = 70.0;
        }
        initial_age_vec(n) = (arma::uword) sample_normal_truncated(mu_bar,
                                                                   sigma_sq_bar,
                                                                   a, b);
    }
    full_age_vec = initial_age_vec;
    arma::vec full_age_vec_double = arma::conv_to<arma::vec>::from(
        full_age_vec);
    // log transform
    // full_age_vec_double = arma::log(full_age_vec_double);
    // Perform rescaling
    if (rescale) {
        rescale_factor = arma::stddev(full_age_vec_double);
        std::cout << "assigned rescale_factor" << std::endl;
        full_age_vec_double = (
            full_age_vec_double - arma::mean(full_age_vec_double)) /
            arma::stddev(full_age_vec_double);
    } else {
        full_age_vec_double = full_age_vec_double -
            arma::mean(full_age_vec_double);
        rescale_factor = 1.0;
    }
    return full_age_vec_double;
}

arma::vec generate_sex_data(arma::uword T, arma::uword N) {
    // Assume 60% female
    double prob_param = 0.6;
    // do bernoulli generation
    arma::Col<arma::uword> initial_sex_vec(N);
    arma::Col<arma::uword> full_sex_vec(N);
    // std::bernoulli_distribution my_bern_dist(prob_param);
    for (arma::uword n = 1; n <= N; ++n) {
        // sex_vec(n - 1) = my_bern_dist(random_engine);
        initial_sex_vec(n - 1) = sample_bernoulli_distn(prob_param);
    }
    full_sex_vec = initial_sex_vec;
    arma::vec full_sex_vec_double = arma::conv_to<arma::vec>::from(
        full_sex_vec);
    return full_sex_vec_double;
}


arma::mat generate_covariates(arma::uword T, arma::uword N,
                              double & age_rescale_factor) {
    // Age, sex (skip doing BMI for now)
    arma::uword D = 2;
    // Age and assigned sex will end up being floating point.
    arma::mat covariates(N, D);
    covariates.col(0) = generate_age_data(T, N, true, age_rescale_factor);
    covariates.col(1) = generate_sex_data(T, N);
    return covariates;
}

arma::mat generate_covariates_scenario(arma::uword T, arma::uword N,
                                       int seed_value_scenario) {
    arma::arma_rng::set_seed(seed_value_scenario);
    // Age, sex (skip doing BMI for now)
    arma::uword D = 2;
    // Age and assigned sex will end up being floating point.
    arma::mat covariates(N, D);
    double unused_age_rescale_factor;
    covariates.col(0) = generate_age_data(T, N, false,
                                          unused_age_rescale_factor);
    covariates.col(1) = generate_sex_data(T, N);
    return covariates;
}
