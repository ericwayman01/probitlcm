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

#include "file_io.hpp"
#include "latent_state_related.hpp"
#include "run_mcmc.hpp"

#include <string>
#include <fstream>
#include <sstream> // for std::stringstream
#include <iomanip> // for std::setprecision

#include <armadillo>


// helper function

std::string pad_string_with_zeros(int padded_length, std::string fnamestem,
                                  std::string s) {
    s.insert(s.begin(), padded_length - s.size(), '0');
    std::string fname = fnamestem + s + ".txt";
    return fname;
}

// major functions

void write_datagen_values(const OtherVals & othervals,
                          const DatagenVals & datagenvals) {
    const arma::uword & T = othervals.dimensions.at("T");
    std::string filename;
    std::string filename_index;
    filename = othervals.replic_path + "/" + "datagen_Xmat.txt";
    datagenvals.Xmat.save(filename, arma::arma_ascii);
    filename = othervals.replic_path + "/" + "datagen_alpha_star.txt";
    datagenvals.alpha_star.save(filename, arma::arma_ascii);    
    filename = othervals.replic_path + "/" + "datagen_alpha.txt";
    datagenvals.alpha.save(filename, arma::arma_ascii);
    filename = othervals.replic_path + "/" + "datagen_Ymat_star.txt";
    datagenvals.Ymat_star.save(filename, arma::arma_ascii);    
    filename = othervals.replic_path + "/" + "datagen_Ymat.txt";
    datagenvals.Ymat.save(filename, arma::arma_ascii);

    filename = othervals.replic_path + "/" + "design_matrix.txt";
    othervals.design_matrix.save(filename, arma::arma_ascii);
    arma::uvec class_numbers = convert_alpha_to_class_numbers(
        datagenvals.alpha, othervals.basis_vector);
    filename = othervals.replic_path + "/" +
        "datagen_class_numbers.txt";
    class_numbers.save(filename, arma::arma_ascii);
}

void write_mcmc_output(MCMCDraws & draws, DatagenVals & datagenvals,
                       OtherVals & othervals,
                       bool is_simulation,
                       bool hyperparam_tuning,
                       int scenario_number_zb) {
    // note that in this function, datagenvals is only used if
    //     perform_inverse_permutations is true.
    // write out results
    arma::uword chain_length_after_burnin =
        othervals.chain_length_after_burnin;
    const arma::uword & T = othervals.dimensions.at("T");
    const arma::uword & N = othervals.dimensions.at("N");
    const arma::uword & J = othervals.dimensions.at("J");
    const arma::uword & K = othervals.dimensions.at("K");
    const arma::uword & H_K = othervals.dimensions.at("H_K");
    std::string filename;
    //// write out kappas
    for (arma::uword j = 1; j <= J; ++j) { 
        std::string fname = pad_string_with_zeros(3, "draws_kappa_",
                                                  std::to_string(j));
        filename = othervals.replic_path + "/" + fname;
        // write out the matrix of kappa draws
        draws.kappa(j - 1).save(filename, arma::arma_ascii);
        // write out the rowvec of the average of the above matrix
        fname = pad_string_with_zeros(3, "average_kappa_", std::to_string(j));
        filename = othervals.replic_path + "/" + fname;
        arma::rowvec average_kappa_j = arma::sum(draws.kappa(j - 1), 0);
        average_kappa_j = average_kappa_j / chain_length_after_burnin;
        average_kappa_j.save(filename, arma::arma_ascii);
    }

    filename = othervals.replic_path + "/" + "draws_omega" + ".txt";
    draws.omega.save(filename, arma::arma_ascii);

    // before permuting, calculate our performance metric
    arma::umat mode_alphas = calc_mode_alphas(draws.class_counts, othervals);
    filename = othervals.replic_path + "/" + "mode_alphas.txt";
    mode_alphas.save(filename, arma::arma_ascii);

    arma::uvec mode_class_nums = convert_alpha_to_class_numbers(
        mode_alphas, othervals.basis_vector);
    arma::field<arma::mat> theta_j_mats_avgs(J);
    for (arma::uword j = 1; j <= J; ++j) {
        arma::mat theta_j_mat_sum = draws.theta_j_mats_sums(j - 1);
        theta_j_mats_avgs(j - 1) = theta_j_mat_sum / chain_length_after_burnin;
    }
    if (is_simulation) {
        double total_metric_value = 0;
        arma::uword mode_class;
        for (arma::uword n = 1; n <= N; ++n) {
            mode_class = mode_class_nums(n - 1);
            for (arma::uword j = 1; j <= J; ++j) {
                total_metric_value += arma::norm(
                    theta_j_mats_avgs(j - 1) .row(mode_class) -
                    datagenvals.theta_j_mats(j - 1).row(mode_class),
                    1);
            }
        }
        // save as matrix for safe loading later
        double metric_value = total_metric_value / (J * N * H_K);
        filename = othervals.replic_path + "/" + "metric_value" + ".txt";
        arma::mat metric_value_mat(1, 1);
        metric_value_mat(0, 0) = metric_value;
        metric_value_mat.save(filename, arma::arma_ascii);
    }
    
    arma::uvec best_inverse_perm_dims; // also used after this if statement
    if (is_simulation) {
        // perform inverse permutation before writing draws to disk
        // find index number of best permutation
        arma::uword index_of_best_inverse_perm =
            find_index_num_of_best_inverse_perm_cb(
                draws.class_counts, datagenvals.alpha, othervals);
        std::cout << "index_of_best_inverse_perm = "
                  << std::to_string(index_of_best_inverse_perm) << std::endl;
        // arma::uword index_of_best_perm = find_index_num_of_best_permutation(
        //     draws.beta, datagenvals.beta, othervals);
        // use this index number to load the inverse permutation
        //     of the pos nums
        arma::umat table_inverse_perms_of_pos_nums = load_umat(
            othervals.scenario_path, "table_inverse_perms_of_pos_nums.txt");
        arma::uvec best_inverse_perm_pos_nums = get_best_perm(
            table_inverse_perms_of_pos_nums,
            index_of_best_inverse_perm);            
        // do the same, i.e. load the best inverse permutation of the class
        //     nums (since the index of the best permutation is the same)
        arma::umat table_inverse_perms_of_class_nums = load_umat(
            othervals.scenario_path, "table_inverse_perms_of_class_nums.txt");
        arma::uvec best_inverse_perm_class_nums = get_best_perm(
            table_inverse_perms_of_class_nums, index_of_best_inverse_perm);
        // perform inverse permutations that use pos nums
        draws.beta.each_slice([best_inverse_perm_pos_nums](arma::mat& X)
            {X = X.rows(best_inverse_perm_pos_nums);});
        draws.delta.each_slice([best_inverse_perm_pos_nums](
                                   arma::Mat<arma::uword> & X)
            {X = X.rows(best_inverse_perm_pos_nums);});
        // perform inverse permutation that uses class nums
        for (arma::uword j = 1; j <= J; ++j) {
            // draws.theta_j_mats_sums(j - 1) = draws.theta_j_mats_sums(
            //     j - 1).rows(best_inverse_perm_class_nums);
            theta_j_mats_avgs(j - 1) = theta_j_mats_avgs(j - 1).rows(
                best_inverse_perm_class_nums);
        }
        draws.class_counts = draws.class_counts.cols(
            best_inverse_perm_class_nums);
        // use the index number to load the inverse permutation
        //     of the dims
        arma::umat table_inverse_perms_of_dims = load_umat(
            othervals.scenario_path, "table_inverse_perms_of_dims.txt");
        best_inverse_perm_dims = get_best_perm(
            table_inverse_perms_of_dims, index_of_best_inverse_perm);
        // now perform inverse permutations that use dims
        // do permutation
        draws.lambda.each_slice([best_inverse_perm_dims](arma::mat & X)
            {X = X.cols(best_inverse_perm_dims);});
        // rescale second row
        draws.lambda.each_slice(
            [age_rescale_factor = datagenvals.age_rescale_factor](
                arma::mat & X)
                {X.row(1) = X.row(1) / age_rescale_factor;}
            );
        draws.Rmat.each_slice([best_inverse_perm_dims](arma::mat & X)
            {X = X.cols(best_inverse_perm_dims);});
        draws.Rmat.each_slice([best_inverse_perm_dims](arma::mat & X)
            {X = X.rows(best_inverse_perm_dims);});
        draws.Sigma.each_slice([best_inverse_perm_dims](arma::mat & X)
            {X = X.cols(best_inverse_perm_dims);});
        draws.Sigma.each_slice([best_inverse_perm_dims](arma::mat & X)
            {X = X.rows(best_inverse_perm_dims);});
    }

    // write gammas
    // note that gammas will be saved as they are permuted, due to
    //     how they are stored as opposed to most other variables
    arma::uword index = 0; // 0-based
    std::string k_string;
    // find new index for gamma vector based on inverse permutation
    for (arma::uword k = 1; k <= K; ++k) {
        if (is_simulation) {
            index = best_inverse_perm_dims(k - 1);
        }
        // index = k - 1;
        k_string = std::to_string(k);
        filename = othervals.replic_path + "/" +
            "draws_gamma_" + k_string + ".txt";
        // write out the matrix of kappa draws, using the proper index
        //     (from the inverse permutation)
        draws.gamma(index).save(filename, arma::arma_ascii);
        // write out the rowvec of the average of the above matrix
        filename = othervals.replic_path + "/" +
            "average_gamma_" + k_string + ".txt";
        arma::rowvec average_gamma_k = arma::sum(
            draws.gamma(index).rows(othervals.burnin,
                                    othervals.total_chain_length - 1), 0);
        average_gamma_k = average_gamma_k / chain_length_after_burnin;
        average_gamma_k.save(filename, arma::arma_ascii);
    }

    write_output_cube(othervals, draws.delta, "delta");    
    write_output_cube_average(othervals, draws.delta, "delta",
                              chain_length_after_burnin);
    write_output_cube(othervals, draws.beta, "beta");
    write_output_cube_average(othervals, draws.beta, "beta",
                              chain_length_after_burnin);

    write_output_cube(othervals, draws.Rmat, "Rmat");
    write_output_cube_average(othervals, draws.Rmat, "Rmat",
                              chain_length_after_burnin);

    write_output_cube(othervals, draws.lambda, "lambda");
    write_output_cube_average(othervals, draws.lambda, "lambda",
                              chain_length_after_burnin);
    
    // draws.alpha_counts.cols(best_permutation);    
    // filename = othervals.replic_path + "/" + "alpha_counts" + ".txt";
    // draws.alpha_counts.save(filename, arma::arma_ascii);

    // //// write out theta_j_mat_avg values
    for (arma::uword j = 1; j <= J; ++j) {
        // // theta_j_mat_avg
        // arma::mat theta_j_mat_sum = draws.theta_j_mats_sums(j - 1);
        // arma::mat theta_j_mat_avg = theta_j_mat_sum / chain_length_after_burnin;
        std::string fname = pad_string_with_zeros(3, "theta_j_mat_avg_",
                                                  std::to_string(j));
        filename = othervals.replic_path + "/" + fname;
        // theta_j_mat_avg.save(filename, arma::arma_ascii);
        theta_j_mats_avgs(j - 1).save(filename, arma::arma_ascii);
    }
    write_output_cube(othervals, draws.Sigma, "Sigma");
    write_output_cube_average(othervals, draws.Sigma, "Sigma",
                              chain_length_after_burnin);
    // write out alpha_counts
    filename = othervals.replic_path + "/" + "class_counts" + ".txt";
    draws.class_counts.save(filename, arma::arma_ascii);
    filename = othervals.replic_path + "/" + "total_class_counts_per_draw" +
        ".txt";
    draws.total_class_counts_per_draw.save(filename, arma::arma_ascii);    

    // write hyperparam_tuning results
    if (hyperparam_tuning) { // later, perhaps change variable name to reflect
                             // difference between scenario number zb and
                             // setup number, which is not zero-based
        std::string s;
        std::string s_prefix;
        if (is_simulation) {
            s = std::to_string(scenario_number_zb + 1);
            s_prefix = "scenario_";
        } else {
            s = std::to_string(scenario_number_zb);
            s_prefix = "setup_";
        }
        s.insert(s.begin(), 3 - s.size(), '0');
        // scale count properly
        int total_samples = 0;
        for (arma::uword j = 1; j <= J; ++j) { // 1-based
            if (othervals.M_j_s(j - 1) > 2) {
                total_samples += 1;
            }
        }
        if (total_samples == 0) {
            std::string info_message =
                "No kappas were ever drawn from the proposal density, "
                "so nothing was recorded.";
            std::cout << info_message << std::endl;
        } else {
            double count_float = static_cast<double>(
                othervals.ystar_and_kappa_accept_count) / total_samples;
            double acceptance_rate = count_float /
                othervals.chain_length_after_burnin;
            std::string fname = s_prefix + s + "_" +
                "sigma_kappa_sq" + "_" + "results.txt";
            std::string fpath = othervals.tuning_path + "/" + fname;
            std::ofstream tuning_file;
            tuning_file.open(fpath, std::ios_base::app);
            std::string content_string = std::to_string(
                othervals.fixed_constants["sigma_kappa_sq"]);
            // deal with stringstream
            std::stringstream stream;
            stream << std::fixed << std::setprecision(2) << acceptance_rate;
            std::string acceptance_rate_str = stream.str();
            content_string = content_string + "," + acceptance_rate_str;
            tuning_file << content_string << std::endl;
            tuning_file.close();
        }
    }
    filename = othervals.replic_path + "/" + "done" + ".txt";
    
    std::ofstream donefile(filename);
    donefile << "replication finished." << std::endl;
}
