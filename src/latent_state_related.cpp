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

// Functions related to converting to and from position numbers and alphas,
//      etc.

arma::urowvec generate_d_b_k_alpha(const arma::urowvec & alpha,
                                   const arma::uword & L_k,
                                   const arma::uword & k) {
    arma::urowvec d_b_k_alpha(L_k);
    arma::uword alpha_k = alpha(k - 1);
    arma::uword l = 0;
    while (l < L_k) {
        if (alpha_k >= l) {
            d_b_k_alpha(l) = 1;
        } else {
            d_b_k_alpha(l) = 0;
        }
        ++l;
    }
    return d_b_k_alpha;
}

arma::umat generate_design_matrix(arma::umat alphas,
                                  arma::uword K,
                                  arma::uvec L_k_s,
                                  arma::uword ord,
                                  arma::uvec pos_to_remove) {
    arma::uword nrows = alphas.n_rows;
    arma::uword H_K = arma::prod(L_k_s);
    arma::umat design_matrix(H_K, H_K);
    // do stuff
    for (arma::uword i = 1; i <= nrows; ++i) {
        arma::field<arma::urowvec> d_b_k_alphas(K);
        arma::urowvec alpha = alphas.row(i - 1);
        for (arma::uword k = 1; k <= K; ++k) {
            d_b_k_alphas(k - 1) = generate_d_b_k_alpha(alpha, L_k_s(k - 1), k);
        }
        arma::urowvec building_alpha = arma::kron(d_b_k_alphas(K - 1 - 1),
                                                  d_b_k_alphas(K - 1));
        if (K > 2) {
            for (arma::uword k = K - 2; k >= 1; --k) {
                building_alpha = arma::kron(d_b_k_alphas(k - 1),
                                            building_alpha);
            }
        }
        design_matrix.row(i - 1) = building_alpha;
    }
    if (ord < K) {
        design_matrix.shed_cols(pos_to_remove);
    }
    return design_matrix;
}

arma::uvec calculate_basis_vector(arma::uvec Lk_s, arma::uword K) {
    arma::uvec basis_vector(K);
    basis_vector(K - 1) = 1;
    arma::uword z = 0;
    if (K > 1) {
        while(z <= K - 2) {
            basis_vector(z) = arma::prod(Lk_s.subvec(z+1, K - 1));
            ++z;
        }
    }
    return basis_vector;
}

arma::uvec convert_alpha_to_class_numbers(const arma::umat & alpha,
                                          const arma::uvec & basis_vector) {
    return alpha * basis_vector;
}

arma::mat get_design_vectors_from_alpha(
                const arma::umat & alpha,
                const arma::umat & design_matrix,
                const arma::uvec & basis_vector) {
    arma::uvec class_numbers = convert_alpha_to_class_numbers(
        alpha, basis_vector);
    arma::uvec class_numbers_2 = arma::conv_to<arma::uvec>::from(
        class_numbers);
    arma::umat d = design_matrix.rows(class_numbers_2);
    arma::mat d_double = arma::conv_to<arma::mat>::from(d);
    return d_double;
}

arma::urowvec convert_position_number_to_alpha(arma::uword position_number,
                                               arma::uword K,
                                               arma::uvec L_k_s) {
    arma::uword z = K - 1;
    arma::uword N = position_number;
    arma::urowvec alpha(K);
    while (z > 0) {
        alpha(z) = N % L_k_s(z);
        N = N / L_k_s(z);
        --z;
    }
    alpha(0) = N;
    return alpha;
}
