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

#include "struct_datagenvals.hpp"
#include "struct_othervals.hpp"

#include "datagen.hpp"
#include "file_io.hpp"
#include "generate_covariates.hpp"
#include "latent_state_related.hpp"
#include "run_mcmc.hpp"
#include "sampling_generic.hpp"
      
void do_datagen_and_write_values(const OtherVals & othervals,
                                 DatagenVals & datagenvals,
                                 int scenario_number_zb) {
    std::cout << "begin datagen" << std::endl;
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & K = othervals.dimensions.at("K");
    const arma::uword & J = othervals.dimensions.at("J");
    const arma::uword & H_K = othervals.dimensions.at("H_K");
    datagenvals.omega = 0.5;
    datagenvals.Rmat.load(othervals.scenario_datagen_params_path + "/" +
                          "datagen_Rmat.txt");
    datagenvals.Sigma = datagenvals.Rmat;
    datagenvals.lambda.load(
        othervals.scenario_datagen_params_path + "/" + "datagen_lambda.txt");
    datagenvals.lambda = datagenvals.lambda;
    datagenvals.lambda_expa = datagenvals.lambda;
    datagenvals.beta.load(othervals.scenario_datagen_params_path + "/" +
                          "datagen_beta.txt");
    datagenvals.delta.load(othervals.scenario_datagen_params_path + "/" +
                           "datagen_delta.txt");
    // gammas
    datagenvals.gamma = arma::field<arma::mat>(K);
    std::string fpath;
    for (arma::uword k = 1; k <= K; ++k) {
        arma::uword L_k = othervals.L_k_s(k - 1);
        // initialize full structures
        datagenvals.gamma(k - 1) = arma::mat(1, L_k + 1);
        fpath = othervals.scenario_datagen_params_path + "/" +
            "datagen_gamma_" +
            std::to_string(k) + ".txt";
        arma::rowvec gamma_k;
        gamma_k.load(fpath);
        datagenvals.gamma(k - 1).row(0) = gamma_k;
    }
    // // load kappas
    // kappas
    datagenvals.kappa = arma::field<arma::mat>(J);
    for (arma::uword j = 1; j <= J; ++j) {
        arma::uword M_j = othervals.M_j_s(j - 1);
        // initialize full structures
        datagenvals.kappa(j - 1) = arma::mat(1, M_j + 1);
        std::string fname = pad_string_with_zeros(3, "datagen_kappa_j_",
                                                  std::to_string(j));
        fpath = othervals.scenario_datagen_params_path + "/" + fname;
        arma::rowvec kappa_j;
        kappa_j.load(fpath);
        datagenvals.kappa(j - 1).row(0) = kappa_j;
    }
    arma::mat tmpmat = generate_covariates(T, N,
                                           datagenvals.age_rescale_factor);
    arma::mat intercept(N, 1);
    intercept.ones();
    datagenvals.Xmat = arma::join_rows(intercept, tmpmat);
    // generate alpha stars
    datagenvals.alpha_star = arma::mat(N, K);
    datagenvals.alpha_star_expa = arma::mat(N, K);
    datagenvals.alpha = arma::umat(N, K);
    // initialize alpha_star and write it to disk
    initialize_alpha_and_alpha_star_cross_sec(othervals,
                                              datagenvals.gamma,
                                              datagenvals.Xmat,
                                              datagenvals.lambda,
                                              datagenvals.Rmat,
                                              datagenvals.alpha_star,
                                              datagenvals.alpha_star_expa,
                                              datagenvals.alpha);
    // generate Ymat_star
    arma::mat d_of_alpha = get_design_vectors_from_alpha(
        datagenvals.alpha, othervals.design_matrix,
        othervals.basis_vector);
    datagenvals.Ymat_star = arma::mat(N, J);
    arma::mat d_times_beta = d_of_alpha * datagenvals.beta;
    for (arma::uword n = 1; n <= N; ++n) {
        for (arma::uword j = 1; j <= J; ++j) {
            datagenvals.Ymat_star(
                n - 1,
                j - 1) = arma::randn(
                    arma::distr_param(d_times_beta(n - 1,
                                                   j - 1), 1.0));
        }
    }
    datagenvals.Ymat = arma::Mat<arma::uword>(N, J);
    sample_ordinal_values(datagenvals.Ymat,
                          datagenvals.Ymat_star,
                          datagenvals.kappa,
                          othervals.M_j_s,
                          N,
                          T,
                          J);
    // //// initialize and then calculate values for theta_j matrices
    datagenvals.theta_j_mats = arma::field<arma::mat>(J);
    for (arma::uword j = 1; j <= J; ++j) {
        datagenvals.theta_j_mats(j - 1) = arma::mat(H_K,
                                                    othervals.M_j_s(j - 1));
        // initialize it to all zeros for now; these will be overwritten
        //     in the following for loop
        std::string fname = pad_string_with_zeros(3, "datagen_theta_j_mat_",
                                                  std::to_string(j));
        fpath = othervals.scenario_datagen_params_path + "/" + fname;
        arma::mat theta_j;
        theta_j.load(fpath);
        datagenvals.theta_j_mats(j - 1) = theta_j;
    }
    // write datagen values
    write_datagen_values(othervals, datagenvals);
    std::cout << "end datagen" << std::endl;
}
