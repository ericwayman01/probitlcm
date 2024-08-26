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

#include "initializations.hpp"
#include "latent_state_related.hpp"
#include "sampling_generic.hpp"

#include <armadillo>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // type conversion

namespace py = pybind11;

void find_initial_alphas(std::string process_dir, std::string path_to_Ymat,
                         const int T, const int N, const int K,
                         arma::uvec L_k_s, const int seed_value) {
    std::vector<int> L_k_s_vec = arma::conv_to< std::vector<int> >::from(L_k_s);
    py::function find_initial_alphas_py =
        py::reinterpret_borrow<py::function>(
            py::module::import(
                "probitlcm.exec_numpy").attr("find_initial_alphas")
            );
    find_initial_alphas_py(process_dir, path_to_Ymat,
                           T, N, K, L_k_s_vec, seed_value);
}


void initialize_star_variable(arma::mat & Ymat_star_current,
                              const arma::field<arma::mat> & thresholds_field,
                              const arma::uvec & levels_vec,
                              const arma::uword & num_of_dims,
                              const arma::umat & Ymat,
                              const arma::uword & T, const arma::uword & N) {
    // nomenclature used here is for the Ymat and Ymat_star case
    arma::uword y_nj;
    arma::uword idx;
    arma::uword M_j;
    arma::rowvec kappa_j;
    for (arma::uword j = 1; j <= num_of_dims; ++j) {
        M_j = levels_vec(j - 1);
        kappa_j = thresholds_field(j - 1).row(0);
        for (arma::uword n = 1; n <= N; ++n) {
            idx = n - 1;
            y_nj = Ymat(idx, j - 1);
            if (y_nj == 0) {
                Ymat_star_current(idx, j - 1) = -0.3;
            } else if (y_nj == M_j - 1) {
                Ymat_star_current(idx, j - 1) = kappa_j(M_j - 1) + 0.3;
            } else {
                Ymat_star_current(idx, j - 1) =
                    (kappa_j(y_nj + 1) - kappa_j(y_nj)) / 2.0;
            }
        }
    }
}

void initialize_mcmc_variables(const OtherVals & othervals,
                               MCMCDraws & draws,
                               const arma::umat & Ymat,
                               const arma::mat & Xmat,
                               int seed_value,
                               const DatagenVals & datagenvals,
                               bool is_simulation) {
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & J = othervals.dimensions.at("J");
    const arma::uword & K = othervals.dimensions.at("K");
    const arma::uword & H = othervals.dimensions.at("H");
    const arma::uword & D = othervals.dimensions.at("D");
    const arma::uword & H_K = othervals.dimensions.at("H_K");
    const arma::uword & total_chain_length = othervals.total_chain_length;
    // initializations for original model
    std::string path_to_Ymat;
    if (is_simulation) {
        path_to_Ymat = othervals.replic_path + "/" + "datagen_Ymat.txt";
    } else { // is a data analysis
        path_to_Ymat = othervals.dataset_dir + "/" + "responses.txt";
    }
    // set up Ymat_pred_chunk (draws from posterior predictive distn)
    if (!is_simulation) {
        draws.Ymat_pred_chunk = arma::ucube(N, J, othervals.stream_max_val);
    }
    find_initial_alphas(othervals.replic_path, path_to_Ymat,
                        T, N, K,
                        othervals.L_k_s, seed_value);
    draws.alpha_current = arma::umat(N, K);
    draws.alpha_previous = arma::umat(N, K);
    draws.alpha_current.load(
        othervals.replic_path + "/" + "initial_alphas.txt");
    // draws.alpha_current = datagenvals.alpha;

    arma::mat Ymat_doubles = arma::conv_to<arma::mat>::from(Ymat);
    arma::mat initial_design_vectors = get_design_vectors_from_alpha(
        draws.alpha_current, othervals.design_matrix, othervals.basis_vector);
    
    draws.beta = arma::cube(H, J, total_chain_length);
    // estimate beta
    draws.beta.slice(0) = arma::solve(initial_design_vectors, Ymat_doubles);
    // make beta non-negative
    draws.beta.slice(0).transform( [](double val) {
        return (val < 0.0) ? 0.0 : val; } );
    // draws.beta.slice(0) = datagenvals.beta;

    draws.delta = arma::Cube<arma::uword>(H, J, total_chain_length);
    draws.delta.slice(0).ones(); // initialize delta to all ones
    
    draws.kappa = arma::field<arma::mat>(J);
    prepare_thresholds(draws.kappa,
                       othervals.M_j_s,
                       J,
                       total_chain_length,
                       1.6, 2.4);

    draws.Ymat_star_current = arma::mat(N, J);
    initialize_star_variable(draws.Ymat_star_current,
                             draws.kappa,
                             othervals.M_j_s,
                             J,
                             Ymat,
                             T, N);
    // draws.Ymat_star_current = datagenvals.Ymat_star;
    draws.Ymat_star_previous = arma::mat(N, J);

    arma::mat alpha_doubles = arma::conv_to<arma::mat>::from(
        draws.alpha_current);
    draws.lambda = arma::cube(D, K, total_chain_length);
    draws.lambda.slice(0) = arma::solve(Xmat, alpha_doubles);
    // draws.lambda.slice(0) = datagenvals.lambda;
    // lambda_expa
    draws.lambda_expa = arma::cube(D, K, total_chain_length);
    draws.lambda_expa.slice(0) = draws.lambda.slice(0);
    // draws.lambda_expa.slice(0) = datagenvals.lambda;

    // initialize gamma_k's
    draws.gamma = arma::field<arma::mat>(K);
    draws.gamma_expa = arma::field<arma::mat>(K);
    prepare_thresholds(draws.gamma,
                       othervals.L_k_s,
                       K,
                       total_chain_length,
                       1.6, 2.4);
    for (arma::uword k = 1; k <= K; ++k) {
        arma::uword L_k = othervals.L_k_s(k - 1);
        draws.gamma_expa(k - 1) = arma::mat(total_chain_length, L_k + 1);
        draws.gamma_expa(k - 1).row(0) = draws.gamma(k - 1).row(0);
    }

    draws.alpha_star_current = arma::mat(N, K);
    initialize_star_variable(draws.alpha_star_current,
                             draws.gamma,
                             othervals.L_k_s,
                             K,
                             draws.alpha_current,
                             T, N);
    // draws.alpha_star_current = datagenvals.alpha_star;
    
    draws.alpha_star_previous = arma::mat(N, K);
    draws.alpha_star_expa_current = arma::mat(N, K);
    draws.alpha_star_expa_current = draws.alpha_star_current;
    draws.alpha_star_expa_previous = arma::mat(N, K);

    draws.Rmat = arma::cube(K, K, total_chain_length);
    arma::mat I_K = arma::eye(K, K);
    draws.Rmat.slice(0) = I_K;
    // draws.Rmat.slice(0) = datagenvals.Rmat;
    
    draws.Sigma = arma::cube(K, K, total_chain_length);
    draws.Sigma.slice(0) = draws.Rmat.slice(0);
    
    draws.omega = arma::vec(total_chain_length);
    draws.omega(0) = 0.5;
    // draws.omega(0) = datagenvals.omega;

    // misc stuff
    // initialize derived parameters
    draws.theta_j_mats_sums = arma::field<arma::mat>(J);
    for (arma::uword j = 1; j <= J; ++j) {
        draws.theta_j_mats_sums(j - 1) = arma::mat(
            H_K, othervals.M_j_s(j - 1));
        draws.theta_j_mats_sums(j - 1).zeros();
    }
    draws.class_counts = arma::umat(N, H_K);
    // initialize total_class_counts_per_draw
    draws.total_class_counts_per_draw = arma::umat(total_chain_length, H_K);
    draws.total_class_counts_per_draw.zeros();
    arma::uvec current_r_alpha_class_numbers = convert_alpha_to_class_numbers(
            draws.alpha_current, othervals.basis_vector);
    // set total class counts to all zeros at the start
    draws.total_class_counts_per_draw.row(
        0) = calc_total_class_counts_for_single_draw(
            othervals, current_r_alpha_class_numbers, 0);
}
