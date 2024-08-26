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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <carma>
#include <armadillo>

#include <string>

#include "generate_covariates.hpp"
#include "latent_state_related.hpp"
#include "run_mcmc.hpp"
#include "crosscounts.hpp"

arma::mat load_arma_mat_np(std::string path) {
    arma::mat mymat;
    mymat.load(path);
    return mymat;
}

arma::umat load_arma_umat_np(std::string path) {
    arma::umat mymat;
    mymat.load(path);
    return mymat;
}

arma::cube load_arma_cube_np(std::string path) {
    arma::cube mycube;
    mycube.load(path);
    return mycube;
}

arma::ucube load_arma_ucube_np(std::string path) {
    arma::ucube mycube;
    mycube.load(path);
    return mycube;
}


void save_arma_uvec_np(arma::uvec my_uvec, std::string fpath) {
    my_uvec.save(fpath, arma::arma_ascii);
}

void save_arma_umat_np(arma::umat my_umat, std::string fpath) {
    my_umat.save(fpath, arma::arma_ascii);
}

void save_arma_mat_np(arma::mat my_mat, std::string fpath) {
    my_mat.save(fpath, arma::arma_ascii);
}

py::array_t<double> generate_covariates_scenario_np(
            arma::uword T, arma::uword N, int seed_value_scenario) {
    arma::mat covariates = generate_covariates_scenario(T, N,
                                                        seed_value_scenario);
    return carma::mat_to_arr(covariates, true);
}

py::array_t<unsigned int> generate_design_matrix_np(
            arma::umat alphas,
            arma::uword K,
            arma::uvec L_k_s,
            arma::uword ord,
            arma::uvec pos_to_remove) {
    arma::umat mytestdawg = generate_design_matrix(alphas,
                                                   K, L_k_s, ord,
                                                   pos_to_remove);
    return carma::mat_to_arr(mytestdawg, true);
}

py::array_t<double> get_design_vectors_from_alpha_np(arma::umat alpha,
                                                     arma::umat design_matrix,
                                                     arma::uvec basis_vector) {
    arma::mat mymat = get_design_vectors_from_alpha(alpha,
                                                    design_matrix,
                                                    basis_vector);
    return carma::mat_to_arr(mymat, true);
}

// py::array_t<unsigned int> convert_alpha_to_class_numbers_np(
//                 arma::umat alpha,
//                 arma::uvec basis_vector) {
//     arma::uvec mymat = convert_alpha_to_class_numbers(alpha, basis_vector);
//     return carma::mat_to_arr(mymat, true);
// }

py::array_t<unsigned int> aggregateAllCounts_np(arma::umat x,
                                                arma::uvec nvalues) {
    arma::umat myvec = aggregateAllCounts(x, nvalues);
    return carma::mat_to_arr(myvec, true);
}

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        C++ Extension Module
        --------------------

        .. currentmodule:: probitlcm._core

        .. autosummary::
           :toctree: _generate

           run_replication
           run_data_analysis_chain
           load_arma_mat_np
           load_arma_umat_np
           load_arma_cube_np
           load_arma_ucube_np
           save_arma_uvec_np
           save_arma_umat_np
           save_arma_mat_np
           generate_design_matrix_np
           generate_covariates_scenario_np
           calculate_basis_vector
           get_design_vectors_from_alpha_np
           convert_alpha_to_class_numbers
           aggregateAllCounts_np
    )pbdoc";
    
    m.def("run_replication", &run_replication, R"pbdoc(
        Run a replication

        Some other explanation about the run_replication function.
    )pbdoc");

    m.def("run_data_analysis_chain", &run_data_analysis_chain, R"pbdoc(
        Run a single chain for a data analysis

        Some other explanation about the run_data_analysis_chain.
    )pbdoc");
   
    m.def("load_arma_mat_np", &load_arma_mat_np);
    m.def("load_arma_umat_np", &load_arma_umat_np);
    m.def("load_arma_cube_np", &load_arma_cube_np);
    m.def("load_arma_ucube_np", &load_arma_ucube_np);

    m.def("save_arma_uvec_np", &save_arma_uvec_np);
    m.def("save_arma_umat_np", &save_arma_umat_np);
    m.def("save_arma_mat_np", &save_arma_mat_np);

    m.def("generate_design_matrix_np", &generate_design_matrix_np, R"pbdoc(
        Wrapper function for generate_design_matrix.

        More explanation.
    )pbdoc");
    m.def("generate_covariates_scenario_np", &generate_covariates_scenario_np);

    m.def("calculate_basis_vector", &calculate_basis_vector);
    m.def("get_design_vectors_from_alpha_np",
          &get_design_vectors_from_alpha_np);
    m.def("convert_alpha_to_class_numbers",
          &convert_alpha_to_class_numbers);
    m.def("aggregateAllCounts_np",
          &aggregateAllCounts_np);
}
