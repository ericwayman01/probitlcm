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

#ifndef HEADERS_STRUCT_OTHERVALS_HPP_
#define HEADERS_STRUCT_OTHERVALS_HPP_

#include <armadillo>

struct OtherVals {
    // path related
    std::string jsonfilename_stem;
    int replicnum;
    std::string scenario_path;
    std::string scenario_datagen_params_path;
    std::string replic_path;
    std::string dataset_dir; // only in use for data analysis code
    std::string tuning_path; // only in use for simulations
                             //     with hyperparameter tuning
    // chain length related
    arma::uword burnin;
    arma::uword chain_length_after_burnin;
    arma::uword total_chain_length;
    // dimensions and levels related
    arma::Col<arma::uword> M_j_s;
    arma::Col<arma::uword> L_k_s;
    std::map<std::string, arma::uword> dimensions;
    arma::uword order;
    arma::uvec pos_to_remove; // in use iff order < K
    // fixed constants
    std::map<std::string, double> fixed_constants;
    // misc stuff
    arma::uword q;
    double rho;
    // latent state formulation related
    arma::umat design_matrix;
    arma::uvec basis_vector;
    // simulation related
    double beta_element_main_effect;
    double beta_element_interaction_effect;
    // for posterior predictive checks
    // these are all 1-based
    arma::uword stream_number;
    arma::uword stream_slice_ctr;
    arma::uword stream_max_val;
    arma::uword thinning_interval;
    // ystar_and_kappa draw accept counts
    arma::uword ystar_and_kappa_accept_count;
};

#endif
