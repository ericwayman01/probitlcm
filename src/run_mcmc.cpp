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

#include "datagen.hpp"
#include "file_io.hpp"
#include "generate_covariates.hpp"
#include "initializations.hpp"
#include "latent_state_related.hpp"
#include "log_of_norm_cdf.hpp"
#include "run_mcmc.hpp"
#include "sampling.hpp"
#include "sampling_cowles.hpp"
#include "sampling_generic.hpp"

#include <armadillo>
#include <nlohmann/json.hpp>

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <cmath> // for std::pow

arma::field<arma::vec> compute_V_quantities(const arma::mat& Sigma) {
    arma::field<arma::vec> diagonals(2);
    arma::vec Sigma_diag = Sigma.diag();
    arma::vec V_sqrt_vec = arma::sqrt(Sigma_diag);
    diagonals(0) = V_sqrt_vec;
    arma::vec V_negative_sqrt_vec = 1 / V_sqrt_vec;
    diagonals(1) = V_negative_sqrt_vec;
    return diagonals;
}

arma::uvec calc_vec_of_indices(const arma::uword & k, const arma::uword & K) {
    arma::uvec vec_of_indices;
    if (K == 2) {
        if (k == 1) {
            vec_of_indices = {1};
        } else { // k = 2
            vec_of_indices = {0};
        }
    } else { // K > 2
        if (k == 1) {
            vec_of_indices = arma::linspace<arma::uvec>(1, K - 1, K - 1);
        } else if (k == K) {
            vec_of_indices = arma::linspace<arma::uvec>(0, K - 2, K - 1);
        } else { // 1 < k < K
            arma::uvec vec_of_indices_part_1 = arma::linspace<arma::uvec>(
                0, k - 2, k - 1);
            arma::uvec vec_of_indices_part_2 = arma::linspace<arma::uvec>(
                k, K - 1, K - k);
            vec_of_indices = join_cols(vec_of_indices_part_1,
                                       vec_of_indices_part_2);
        }
    }
    return vec_of_indices;
}

void prepare_thresholds(arma::field<arma::mat> & gamma,
                        const arma::Col<arma::uword> & L_k_s,
                        const arma::uword K,
                        arma::uword n_rows,
                        double unif_lower_bound, double unif_upper_bound) {
    for (arma::uword k = 1; k <= K; ++k) {
        arma::uword L_k = L_k_s(k - 1);
        // initialize full structures
        gamma(k - 1) = arma::mat(n_rows, L_k + 1);
        // initialize specific values for each
        arma::rowvec gamma_k(L_k + 1);
        gamma_k(0) = -1.0 * arma::datum::inf;
        gamma_k(1) = 0.0;
        gamma_k(L_k) = arma::datum::inf;
        // if L_k > 2, initialize other values
        if (L_k > 2) {
            double distr_param1, distr_param2;
            distr_param1 = unif_lower_bound / L_k;
            distr_param2 = unif_upper_bound / L_k;
            arma::rowvec myranvec = arma::randu<arma::rowvec>(
                L_k - 2, arma::distr_param(distr_param1, distr_param2));
            for (arma::uword l = 2; l <= L_k - 1; ++l) { // 0-based
                gamma_k(l) = gamma_k(l - 1) + myranvec(l - 2);
            }
        }
        gamma(k - 1).row(0) = gamma_k;
    }
}

// used in both datagen and initializations
void initialize_alpha_and_alpha_star_cross_sec(
        const OtherVals & othervals,
        const arma::field<arma::mat> & gamma,
        const arma::mat & Xmat,
        const arma::mat & lambda,
        const arma::mat & Rmat,
        arma::mat & alpha_star,
        arma::mat & alpha_star_expa,
        arma::umat & alpha) {
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & K = othervals.dimensions.at("K");
    // build alpha_star and alpha
    arma::mat alpha_star_mean(N, K);
    // build alpha star
    alpha_star_mean = Xmat * lambda;
    for (arma::uword n = 1; n <= N; ++n) {
        alpha_star.row(n - 1) = arma::trans(
            arma::mvnrnd(arma::trans(alpha_star_mean.row(n - 1)),
                         Rmat));
    }
    alpha_star_expa = alpha_star;
    // build alpha
    alpha = arma::umat(N, K);
    sample_ordinal_values(alpha,
                          alpha_star,
                          gamma,
                          othervals.L_k_s,
                          N,
                          T, // T is 1 here
                          K);
}

void perform_transformations(const OtherVals & othervals,
                             MCMCDraws & draws,
                             const arma::uword & draw_number) {
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & K = othervals.dimensions.at("K");
    arma::field<arma::vec> diagonals = compute_V_quantities(
        draws.Sigma.slice(draw_number));
    arma::vec V_sqrt_vec = diagonals(0);
    arma::vec V_negative_sqrt_vec = diagonals(1);
    arma::mat V_sqrt =  arma::diagmat(V_sqrt_vec);
    arma::mat V_negative_sqrt = arma::diagmat(V_negative_sqrt_vec);
    // calculate complete data model quantities
    //// deal with gammas
    for (arma::uword k = 1; k <= K; ++k) {
        draws.gamma(k - 1).row(draw_number) = 
            draws.gamma_expa(k - 1).row(draw_number) *
            V_negative_sqrt_vec(k - 1);
    }
    //// regular variables
    draws.alpha_star_current = draws.alpha_star_expa_current * V_negative_sqrt;
    draws.Rmat.slice(draw_number) = V_negative_sqrt *
        draws.Sigma.slice(draw_number) * V_negative_sqrt;
    draws.lambda.slice(draw_number) = draws.lambda_expa.slice(draw_number) *
            V_negative_sqrt;
}

void set_up_othervals(OtherVals & othervals,
                      const nlohmann::json & json_object,
                      std::string scenario_or_chain_path, bool is_simulation) {
    // chain length related
    othervals.chain_length_after_burnin = json_object.at(
        "chain_length_after_burnin");
    othervals.burnin = json_object.at("burnin");
    othervals.total_chain_length = othervals.chain_length_after_burnin +
        othervals.burnin;
    // vectors
    std::vector<arma::uword> tmpvec;
    tmpvec = json_object.at("L_k_s").get<std::vector<arma::uword>>();
    othervals.L_k_s = arma::uvec(tmpvec);
    // Dimensionalities
    std::map<std::string, arma::uword> dimensions;
    if (is_simulation) {
        othervals.dimensions["N"] = json_object.at("N");
        othervals.dimensions["J"] = json_object.at("J");
        if (json_object.at("covariates") == "age_assignedsex") {
            othervals.dimensions["D"] = 3;
        } else if (json_object.at("covariates") == "age_assignedsex_medical") {
            othervals.dimensions["D"] = 4;
        }
    }
    othervals.dimensions["K"] = json_object.at("K");
    const arma::uword & T = othervals.dimensions.at("T");
    // Order
    othervals.order = json_object.at("order");
    if (othervals.order < othervals.dimensions["K"]) {
        std::string filename = scenario_or_chain_path + "/" +
            "pos_to_remove.txt";
        arma::uvec pos_to_remove;
        pos_to_remove.load(filename);
        othervals.pos_to_remove = pos_to_remove;
    }
    arma::umat alphas;
    std::string fpath = othervals.scenario_path + "/" + "alphas_table.txt";
    alphas.load(fpath);
    othervals.design_matrix = generate_design_matrix(alphas,
                                                     othervals.dimensions["K"],
                                                     othervals.L_k_s,
                                                     othervals.order,
                                                     othervals.pos_to_remove);
    othervals.basis_vector = calculate_basis_vector(othervals.L_k_s,
                                                    othervals.dimensions["K"]);
    othervals.dimensions["H"] = othervals.design_matrix.n_cols;
    othervals.dimensions["H_K"] = arma::prod(othervals.L_k_s);
    // params
    // fixed constants
    std::vector<std::string> fixed_constants_list =
        {"omega_0", "omega_1", "sigma_kappa_sq", "sigma_beta_sq"};
    for (auto & x : fixed_constants_list) {
        othervals.fixed_constants[x] = json_object.at(x);
    }
    if (is_simulation) {
        // misc stuff
        othervals.rho = json_object.at("rho");
        othervals.q = json_object.at("q");
        const arma::uword & q = othervals.q;
        //set up M_j_s from per_effect_M_j_s
        othervals.M_j_s = arma::Col<arma::uword>(othervals.dimensions["J"]);
        // clear tmpvec so it can be used again here
        tmpvec.clear();
        tmpvec = json_object.at(
            "per_effect_M_j_s").get<std::vector<arma::uword>>();
        arma::uword z = 0;
        arma::uword l = 0;
        for (auto x : tmpvec) {
            while (l < q) {
                othervals.M_j_s(z * q + l) = x;
                ++l;
            }
            l = 0;
            ++z;
        }
        // simulation related
        othervals.beta_element_main_effect = 2.0;
        othervals.beta_element_interaction_effect = 1.0;
        // done
    }
    // set up Ymat_pred_chunk (draws from posterior predictive distn)
    if (!is_simulation) {
        othervals.stream_number = 1;
        othervals.stream_slice_ctr = 1;
        othervals.stream_max_val = 100;
        othervals.thinning_interval = 10;
    }
    // hyperparameter tuning-related
    othervals.ystar_and_kappa_accept_count = 0;
}


// begin find_best_permutation function section

arma::uword factorial(arma::uword n) {
    if (n < 0 ) {
        return 0;
    }
    return !n ? 1 : n * factorial(n - 1);
}

// from https://stackoverflow.com/questions/13558081/return-index-of-smallest-element-in-array
int indexofSmallestElement(double array[], arma::uword size) {
    int index = 0;
    int size_int = (int) size;
    for(int i = 1; i < size_int; i++) {
        if(array[i] < array[index])
            index = i;
    }
    return index;
}

arma::umat load_umat(std::string scenario_path, std::string fname) {
    std::string filename = scenario_path + "/" + fname;
    arma::umat my_mat;
    my_mat.load(filename);
    return my_mat;
}

arma::uvec get_best_perm(const arma::umat & my_perms_mat,
                         arma::uword index_of_best_perm) {
    return arma::trans(my_perms_mat.row(index_of_best_perm));
}

arma::umat calc_mode_alphas(arma::umat & class_counts, OtherVals & othervals) {
    // prepare aliases
    arma::uword T = othervals.dimensions["T"];
    arma::uword N = othervals.dimensions["N"];
    arma::uword K = othervals.dimensions["K"];    
    arma::uvec L_k_s = othervals.L_k_s;
    // do work
    arma::umat mode_alphas(N, K);
    arma::uword idx;
    arma::uword mode_class;
    for (arma::uword n = 1; n <= N; ++n) {
        idx = n - 1;
        mode_class = arma::index_max(class_counts.row(idx));
        mode_alphas.row(idx) = convert_position_number_to_alpha(mode_class,
                                                                K,
                                                                L_k_s);
    }
    return mode_alphas;
}

// based on pos nums. "cb" means "class-based."
arma::uword find_index_num_of_best_inverse_perm_cb(arma::umat & class_counts,
                                                   arma::umat & datagen_alpha,
                                                   OtherVals & othervals) {
    // prepare variables
    arma::uword T = othervals.dimensions["T"];
    arma::uword N = othervals.dimensions["N"];
    arma::uvec L_k_s = othervals.L_k_s;
    // find mode_alphas matrix
    arma::umat mode_alphas;
    std::string fpath = othervals.replic_path + "/" + "mode_alphas" + ".txt";
    mode_alphas.load(fpath, arma::arma_ascii);
    // consider all dimension-induced permutations of that mode_alphas matrix
    arma::umat table_inverse_perms_of_dims;
    fpath = othervals.scenario_path + "/" + "table_inverse_perms_of_dims.txt";
    table_inverse_perms_of_dims.load(fpath);
    arma::uword num_of_inverse_perms = table_inverse_perms_of_dims.n_rows;
    arma::urowvec inverse_perm;
    arma::umat potential_alpha;
    arma::uvec metric(num_of_inverse_perms);
    metric.zeros();
    for (arma::uword i = 0; i < num_of_inverse_perms; ++i) {
        inverse_perm = table_inverse_perms_of_dims.row(i);
        potential_alpha = mode_alphas.cols(inverse_perm);
        unsigned int metric_value = 0;
        for (arma::uword j = 0; j < N; ++j) {
            if (arma::all(potential_alpha.row(j) == datagen_alpha.row(j))) {
                metric_value += 1;
            }
        }
        metric(i) = metric_value;
    }
    arma::uword index_of_best_inverse_perm = metric.index_max();
    // write index of best inverse permutation to file
    std::ofstream myfstream(othervals.replic_path + "/" +
                            "index_of_best_inverse_permutation.txt");
    myfstream << std::to_string(index_of_best_inverse_perm);
    return index_of_best_inverse_perm;
}

// based on pos nums
arma::uword find_index_num_of_best_permutation(arma::cube & X_hat,
                                               arma::mat & X,
                                               OtherVals & othervals) {
    // arma::uword K = othervals.dimensions["K"];
    arma::uword H = othervals.dimensions["H"];
    // create uvec of indices post-burnin
    arma::uword burnin = othervals.burnin;
    arma::uword total_chain_length = othervals.total_chain_length;
    arma::uword chain_length_after_burnin = total_chain_length - burnin;
    arma::uvec indices_post_burnin = arma::linspace<arma::uvec>(
        burnin, // burnin + 1 - 1
        total_chain_length - 1,
        chain_length_after_burnin);
    // load table of inverse perms of pos nums
    arma::umat table_inverse_perms_of_pos_nums = load_umat(
        othervals.scenario_path, "table_inverse_perms_of_pos_nums.txt");
    arma::Row<arma::uword> inverse_perm_of_pos_nums(H);
    // arma::uword number_of_permutations = factorial(K);
    // get number of permutations from loaded table
    arma::uword number_of_permutations =
        table_inverse_perms_of_pos_nums.n_rows;
    double mad_values[number_of_permutations];
    const arma::cube & X_hat_post_burnin = X_hat.slices(indices_post_burnin);
    arma::uvec indices;
    // loop over table of inverse perms of pos nums and do mad calculations
    for (arma::uword i = 0; i < number_of_permutations; ++i) {
        inverse_perm_of_pos_nums = table_inverse_perms_of_pos_nums.row(i);
        indices = arma::trans(inverse_perm_of_pos_nums);
        // calculate norm
        double norm_value = 0.0;
        for (arma::uword s = 0; s < chain_length_after_burnin; ++s) {
            norm_value += arma::norm(
                arma::vectorise(X_hat_post_burnin.slice(s).rows(indices) - X),
                1);
        }
        mad_values[i] = norm_value;
    }
    int index_of_best_inverse_permutation = indexofSmallestElement(
        mad_values, number_of_permutations);
    // write index of best inverse permutation to file
    std::ofstream myfstream(othervals.replic_path + "/" +
                            "index_of_best_inverse_permutation.txt");
    myfstream << std::to_string(index_of_best_inverse_permutation);
    return arma::uword(index_of_best_inverse_permutation);
}

// end find_best_permutation function section

double calc_theta_jmc(arma::uword m, const arma::Row<arma::uword> & design_vec,
                      const arma::vec & beta_j, const arma::rowvec & kappa_j) {
    return arma::normcdf(kappa_j(m) - arma::dot(design_vec, beta_j)) -
        arma::normcdf(kappa_j(m - 1) - arma::dot(design_vec, beta_j));
}

void calc_theta_j_quantities(arma::mat & theta_j_mat_sums,
                             const arma::vec & beta_j,
                             const arma::rowvec & kappa_j,
                             const arma::uword M_j,
                             const OtherVals & othervals) {
    const arma::uword & H_K = othervals.dimensions.at("H_K");
    double theta_jmc;
    for (arma::uword c = 1; c <= H_K; ++c) {
        for (arma::uword m = 1; m <= M_j; ++m) {
            theta_jmc = calc_theta_jmc(
                m, othervals.design_matrix.row(c - 1),
                beta_j, kappa_j);
            theta_j_mat_sums(c - 1, m - 1) = theta_j_mat_sums(c - 1, m - 1) +
                theta_jmc;
        }
    }
}

void update_class_counts(MCMCDraws & draws, const OtherVals & othervals,
                         const arma::uvec & current_r_alpha_class_numbers) {
    arma::uword T = othervals.dimensions.at("T");
    arma::uword N = othervals.dimensions.at("N");
    arma::uword tmpnum;
    // Do calculations
    for (arma::uword n = 1; n <= N; ++n) {
        tmpnum = current_r_alpha_class_numbers(n - 1); // 1-based
        draws.class_counts(n - 1, tmpnum) =
            draws.class_counts(n - 1, tmpnum) + 1;
        // 1-based
    }
}

arma::urowvec calc_total_class_counts_for_single_draw(
            const OtherVals & othervals,
            const arma::uvec & current_r_alpha_class_numbers,
            const arma::uword & draw_number) {
    const arma::uword & H_K = othervals.dimensions.at("H_K");
    arma::vec myedges = arma::linspace(-0.5, H_K - 1.5, H_K);
    arma::vec current_class_numbers_vec = arma::conv_to<arma::vec>::from(
        current_r_alpha_class_numbers);
    arma::uvec myhist = arma::hist(current_class_numbers_vec, myedges);
    return arma::trans(myhist);
}

void run_mcmc(DatagenVals & datagenvals,
              OtherVals & othervals, MCMCDraws & draws,
              const arma::umat & Ymat,
              const arma::mat & Xmat,
              bool is_simulation) {
    std::cout << "beginning run_mcmc" << std::endl;
    arma::uword draw_number = 1;
    // note: draw_number is zero-based (the first draw is the initialization
    //           draw). it runs from 0 through othervals.total_chain_length - 1
    //           inclusive.
    std::ofstream log_file(othervals.replic_path + "/" +
                           "log_loops.txt");
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & J = othervals.dimensions.at("J");
    const arma::uword & K = othervals.dimensions.at("K");
    const arma::uword & D = othervals.dimensions.at("D");
    arma::uvec current_r_alpha_class_numbers;
    // begin setup for posterior predictive check code
    arma::mat alpha_star_ppc(N, K);
    arma::umat alpha_ppc(N, K);
    arma::mat design_mat_ppc;
    arma::mat ystar_ppc(N, J);
    arma::mat eye_ppc(J, J);
    eye_ppc.eye();
    // end setup for posterior predictive check code
    // execute sampling
    while (draw_number < othervals.total_chain_length) {
        // set "previous" values for large matrices
        draws.Ymat_star_previous = draws.Ymat_star_current;
        draws.alpha_previous = draws.alpha_current;
        draws.alpha_star_expa_previous = draws.alpha_star_expa_current;
        // perform draws
        //// b1
        sample_ystar_and_kappa(othervals, draws, Ymat, draw_number);
        // draws.Ymat_star_current = datagenvals.Ymat_star;
        // const arma::uword & J = othervals.dimensions.at("J");
        // for (arma::uword j = 1; j <= J; ++j) {
        //     draws.kappa(j - 1).row(draw_number) = datagenvals.kappa(
        //         j - 1).row(0);
        // }
        sample_beta_and_delta(datagenvals, othervals, draws, draw_number);
        // draws.delta.slice(draw_number) = datagenvals.delta;
        // draws.beta.slice(draw_number) = datagenvals.beta;
        sample_alpha_and_alpha_star_expa_crosssec(othervals, draws, Xmat,
                                                  draw_number);
        sample_gamma_expa(othervals, draws, draw_number);
        // const arma::uword & K = othervals.dimensions.at("K");
        // for (arma::uword k = 1; k <= K; ++k) {
        //     draws.gamma_expa(k - 1).row(draw_number) = datagenvals.gamma(
        //         k - 1).row(0);
        //     draws.gamma(k - 1).row(draw_number) = datagenvals.gamma(
        //         k - 1).row(0);
        //     draws.gamma_expa(k - 1).row(draw_number).print("gamma:");
        // }
        sample_Sigma_and_lambda_expa(othervals, draws,
                                     Xmat, draw_number); 
        sample_omega(othervals, draws, draw_number); // b7
        // draws.omega(draw_number) = datagenvals.omega;
        // transform parameters of expanded model to parameters of
        //     complete data model
        perform_transformations(othervals, draws, draw_number);
        arma::uword d_n_plus_one = draw_number + 1;
        if (!is_simulation) { // draw y values from predictive distn
            if ((d_n_plus_one % othervals.thinning_interval == 0)
                && (draw_number > othervals.burnin)) {
                // once we are dealing with the "right" chunk,
                //     sample the slice and update the counters
                alpha_star_ppc = sample_matrix_variate_normal_indep_rows(
                    Xmat * draws.lambda.slice(draw_number),
                    draws.Rmat.slice(draw_number),
                    N, K);                    
                sample_ordinal_values_newer(
                    alpha_ppc, alpha_star_ppc, draws.gamma,
                    othervals.L_k_s, K, draw_number);
                design_mat_ppc = get_design_vectors_from_alpha(
                    alpha_ppc, othervals.design_matrix,
                    othervals.basis_vector);
                ystar_ppc = sample_matrix_variate_normal_indep_rows(
                    design_mat_ppc * draws.beta.slice(draw_number),
                    eye_ppc,
                    N, J);
                sample_ordinal_values_newer(
                    draws.Ymat_pred_chunk.slice(othervals.stream_slice_ctr - 1),
                    ystar_ppc, draws.kappa,
                    othervals.M_j_s, J, draw_number);                
                // if stream_slice_ctr has reached largest possible value,
                //     save chunk and zero out cube so we can begin
                //     the next chunk "cleanly"
                if (othervals.stream_slice_ctr == othervals.stream_max_val) {
                    std::string fname = pad_string_with_zeros(
                        3, "Ymat_pred_chunk_",
                        std::to_string(othervals.stream_number));
                    std::string fpath = othervals.replic_path + "/" + fname;
                    draws.Ymat_pred_chunk.save(fpath, arma::arma_ascii);
                    draws.Ymat_pred_chunk.zeros(); // "clear" the cube
                    // reset the appropriate counters
                    othervals.stream_slice_ctr = 1;
                    othervals.stream_number += 1;
                } else {
                    othervals.stream_slice_ctr += 1;
                }
            }
        }
        if (draw_number > othervals.burnin) {
            for (arma::uword j = 1; j <= J; ++j) {
                calc_theta_j_quantities(draws.theta_j_mats_sums(j - 1),
                                        draws.beta.slice(draw_number).col(
                                            j - 1),
                                        draws.kappa(j - 1).row(draw_number),
                                        othervals.M_j_s(j - 1),
                                        othervals);
            }
        }
        // update alpha_counts
        current_r_alpha_class_numbers = convert_alpha_to_class_numbers(
            draws.alpha_current, othervals.basis_vector);
        if (draw_number > othervals.burnin) {
            update_class_counts(draws, othervals,
                                current_r_alpha_class_numbers);
        }
        // note: does not include draw zero (the initialization)
        draws.total_class_counts_per_draw.row(
            draw_number) = calc_total_class_counts_for_single_draw(
                othervals, current_r_alpha_class_numbers, draw_number);
        // perform logging and increment draw number
        log_file << "finished loop " << std::to_string(draw_number)
                 << std::endl;
        draw_number += 1;
    }
    std::cout << "Chain done." << std::endl;
}

nlohmann::json read_json_file(std::string path_to_json_file) {
    std::ifstream fstream_for_json(path_to_json_file);
    nlohmann::json json_object = nlohmann::json::parse(fstream_for_json);
    return json_object;
}

void run_replication(std::string jsonfilename_stem,
                     int scenario_number_zb,
                     int replicnum,
                     int number_of_replics,
                     std::string other_json_files_path,
                     std::string scenario_path,
                     std::string scenario_datagen_params_path,
                     std::string replic_path,
                     bool hyperparam_tuning,
                     std::string tuning_path) {
    // read in the json file
    std::cout << "tuning true/false: "
              << std::to_string(hyperparam_tuning) << std::endl;
    std::string fulljsonfilename = scenario_path + "/" + jsonfilename_stem +
        ".json";
    nlohmann::json json_object = read_json_file(fulljsonfilename);
    // merge in other json files
    fulljsonfilename = other_json_files_path + "/" + "01_fixed_vals.json";
    nlohmann::json json_object_two = read_json_file(fulljsonfilename);
    json_object.update(json_object_two);
    if (hyperparam_tuning) {
        fulljsonfilename = other_json_files_path + "/" +
            "02_current_tuning_hyperparam_vals.json";
        json_object_two = read_json_file(fulljsonfilename);
        json_object.update(json_object_two);
    } else {
        fulljsonfilename = other_json_files_path + "/" +
            "02_list_decided_hyperparam_vals.json";
        json_object_two = read_json_file(fulljsonfilename);
        nlohmann::json json_object_three = json_object_two.at(
            "sigma_kappa_sq");
        std::vector<int> myvec1 = json_object_three.at(0);
        std::vector<double> myvec2 = json_object_three.at(1);
        // find and assign proper value for sample size
        int i = 0;
        for (auto & element : myvec1) {
            if (element == json_object.at("N")) {
                json_object["sigma_kappa_sq"] = myvec2[i];
            }
            ++i;
        }
    }
    // }
    // // build othervals from json mostly
    OtherVals othervals;
    if (hyperparam_tuning) {
        json_object["chain_length_after_burnin"] = 1000;
        json_object["burnin"] = 0;
        json_object["total_chain_length"] = 1000;
    }
    // set up path related
    othervals.jsonfilename_stem = jsonfilename_stem;
    othervals.replicnum = replicnum;
    othervals.scenario_path = scenario_path;
    othervals.scenario_datagen_params_path = scenario_datagen_params_path;
    othervals.replic_path = replic_path;
    othervals.tuning_path = tuning_path;
    // set up other values
    othervals.dimensions["T"] = json_object.at("T");
    set_up_othervals(othervals, json_object, scenario_path, true);
    // declare datagenvals and draws
    DatagenVals datagenvals;
    MCMCDraws draws;
    // set random seed
    int replicnum_zb = replicnum - 1;
    int seed_value = number_of_replics * scenario_number_zb + replicnum_zb;
    arma::arma_rng::set_seed(seed_value);
    // set up: latent state formulation related quantities
    // do things to run model
    do_datagen_and_write_values(othervals, datagenvals,
                                scenario_number_zb);
    // reset seed
    arma::arma_rng::set_seed(seed_value);
    initialize_mcmc_variables(othervals, draws,
                              datagenvals.Ymat, datagenvals.Xmat,
                              seed_value, datagenvals, true);
    std::cout << "run_mcmc started" << std::endl;
    run_mcmc(datagenvals, othervals, draws, datagenvals.Ymat, datagenvals.Xmat,
             true);
    write_mcmc_output(draws, datagenvals, othervals, true,
                      hyperparam_tuning, scenario_number_zb);
}

void run_data_analysis_chain(int setup_num,
                             std::string dataset_dir,
                             std::string data_analysis_path,
                             std::string chain_results_path,
                             int chainnum,
                             bool hyperparam_tuning,
                             std::string tuning_path) {
    //// read in data files
    arma::umat Ymat;
    arma::field<std::string> empty_header(0);
    std::string fpath = dataset_dir + "/" + "responses.csv";
    Ymat.load(arma::csv_name(fpath, empty_header,
                             arma::csv_opts::strict));
    arma::mat Xmat;
    fpath = dataset_dir + "/" + "covariates.csv";
    Xmat.load(arma::csv_name(fpath, empty_header,
                             arma::csv_opts::strict));
    //// save off copies of these data files in arma_ascii format
    ////     so later they can be loaded using the same functions
    ////     used to load all other results (this will help avoid potential
    ////     breakage in the future due to potential differences in the way CSV
    ////     files are loaded in Python)
    fpath = dataset_dir + "/" + "responses.txt";
    Ymat.save(fpath, arma::arma_ascii);
    fpath = dataset_dir + "/" + "covariates.txt";
    Xmat.save(fpath, arma::arma_ascii);
    // read in the json file
    std::string setup_num_string = std::to_string(setup_num);
    setup_num_string.insert(setup_num_string.begin(),
                            4 - setup_num_string.size(), '0');    
    std::string fulljsonfilename = data_analysis_path + "/" +
        "setup_" + setup_num_string + ".json";
    std::ifstream fstream_for_json(fulljsonfilename);
    nlohmann::json json_object = nlohmann::json::parse(fstream_for_json);
    // merge in other json files
    fulljsonfilename = dataset_dir + "/" + "01_fixed_vals.json";
    nlohmann::json json_object_two = read_json_file(fulljsonfilename);
    json_object.update(json_object_two);
    if (hyperparam_tuning) {
        fulljsonfilename = dataset_dir + "/" +
            "02_current_tuning_hyperparam_vals.json";
        json_object_two = read_json_file(fulljsonfilename);
        json_object.update(json_object_two);
    } else {
        fulljsonfilename = dataset_dir + "/" +
            "02_decided_hyperparam_vals.json";
        json_object_two = read_json_file(fulljsonfilename);
        double hello = json_object_two.at(
            "sigma_kappa_sq");
        json_object["sigma_kappa_sq"] = hello;
    }
    // // build othervals from json mostly
    OtherVals othervals;
    if (hyperparam_tuning) {
        json_object["chain_length_after_burnin"] = 1000;
        json_object["burnin"] = 0;
        json_object["total_chain_length"] = 1000;
    }
    // set up path related
    othervals.jsonfilename_stem = setup_num_string;
    othervals.scenario_path = data_analysis_path;
    othervals.replic_path = chain_results_path;
    othervals.dataset_dir = dataset_dir; // only in use with data analysis code
    othervals.tuning_path = tuning_path;
    //// find N, J, and D
    othervals.dimensions["T"] = json_object.at("T");
    othervals.dimensions["N"] = Ymat.n_rows;
    othervals.dimensions["J"] = Ymat.n_cols;
    othervals.dimensions["D"] = Xmat.n_cols;
    // find M_j_s and set it in othervals
    const arma::uword J = othervals.dimensions.at("J");
    arma::uvec M_j_s(J);
    for (arma::uword j = 1; j <= J; ++j) {
        // in the following line, we must add 1 since the responses are 0-based
        M_j_s(j - 1) = arma::max(Ymat.col(j - 1)) + 1;
    }
    othervals.M_j_s = M_j_s;
    // set up other values
    set_up_othervals(othervals, json_object, data_analysis_path, false);
    // continue
    // declare draws
    MCMCDraws draws;
    // set random seed equal to zero-based chain number
    int seed_value = chainnum - 1;
    arma::arma_rng::set_seed(seed_value);
    DatagenVals datagenvals; // create empty datagenvals struct
    std::cout << "start initialize_mcmc_variables" << std::endl;
    initialize_mcmc_variables(othervals, draws,
                              Ymat, Xmat,
                              seed_value, datagenvals, false);
    std::cout << "start run_mcmc" << std::endl;
    run_mcmc(datagenvals, othervals, draws, Ymat, Xmat, false);
    std::cout << "start write_mcmc_output" << std::endl;
    write_mcmc_output(draws, datagenvals, othervals, false, hyperparam_tuning,
                      setup_num);
}
