# This file is part of "probitlcm" which is released under GPL v3.
#
# Copyright (c) 2022-2024 Eric Alan Wayman <ewayman2@illinois.edu>.
#
# This program is FLO (free/libre/open) software: you can redistribute
# it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## standard library
import pathlib
import itertools
import json

## other modules
import numpy as np

from probitlcm import report_helpers

# translation of
# arma::Col<arma::uword> calculate_basis_vector(arma::Col<arma::uword>,
#                                               arma::uword);
# from latent_state_related.cpp
def calculate_basis_vector(L_k_s, K):
    basis_vector = np.empty(K, dtype=np.uint)
    basis_vector[K - 1] = 1
    z = 0
    if K > 1:
        while (z <= K - 2):
            basis_vector[z] = np.prod(L_k_s[z+1:K])
            z = z + 1
    return basis_vector

# translation (but for single vector only) of
# arma::Col<arma::uword> convert_alpha_to_class_numbers(
#                 arma::Mat<arma::uword>,
#                 arma::Col<arma::uword>);
# from latent_state_related.cpp
def convert_alpha_to_class_number(alpha, basis_vector):
    class_number = int(np.dot(alpha, basis_vector))
    return class_number
    
def create_alphas_table(scenario_dict):
    K = scenario_dict["K"]
    L_k_s = scenario_dict["L_k_s"]
    C = np.prod(L_k_s)
    # create table of possible alphas
    alpha = np.empty((C, K), dtype=np.uint)
    for c in range(C):
        alpha[c, :] = np.unravel_index(c, L_k_s)
    return alpha

def find_pos_to_keep(ord, L_k_s, K):
    my_dict = dict()
    for k in range(1, K + 1):
        my_dict[k] = list(range(L_k_s[k - 1]))
    seq_K = list(range(1, K + 1))
    # calculate which dimensions to eliminate from position vector
    # this can be calculated for any dimension permutation and one
    #     would get the same thing, so calculate it for the identity
    #     permutation
    full_levels_vec_set = set()
    for my_ord in range(1, ord + 1):
        for comb in itertools.combinations(seq_K, my_ord):
            building_blocks = list()
            for j in seq_K:
                if j in comb:
                    building_blocks.append(my_dict[j])
                else:
                    building_blocks.append([0])
            full_levels_vec_set.update(itertools.product(*building_blocks))
    pos_to_keep_set = set()
    basis_vector = calculate_basis_vector(L_k_s, K)
    for x in full_levels_vec_set:
        pos_to_keep_set.add(convert_alpha_to_class_number(x, basis_vector))
    return sorted(pos_to_keep_set)

def build_pos_related_structures(order, L_k_s, K):
    pos_to_keep = list()
    if order < K:
        pos_to_keep = find_pos_to_keep(order, L_k_s, K)
        # note that this list is actually the inverse permutation of the
        #      resulting vector
    else: # order == K 
        pos_to_keep = list(range(0, np.prod(L_k_s)))
        # everything
    # create dict which gives mapping of numbers
    pos_mapping_dict = dict()
    i = 0
    for x in pos_to_keep:
        pos_mapping_dict[x] = i
        i = i + 1
    return (pos_mapping_dict, pos_to_keep)


def pos_to_remove_helper(order, L_k_s, K, H_K):
    pos_mapping_dict, pos_to_keep = build_pos_related_structures(order,
                                                                 L_k_s, K)
    # find pos_to_remove and, if T > 1, pos_to_remove_trans
    seq_H_K = set(range(0, H_K))
    pos_to_remove = sorted(seq_H_K.difference(set(pos_to_keep)))
    H = len(pos_to_keep)
    # part four point five
    return pos_to_remove
    
def find_pos_to_remove_and_effects_tables(scenario_dict, fixed_other_vals_dict,
                                          mydict, alpha):
    T = fixed_other_vals_dict["T"]
    K = scenario_dict["K"]
    L_k_s = scenario_dict["L_k_s"]
    order = fixed_other_vals_dict["order"]
    H_K = np.prod(L_k_s)
    results_dict = dict()
    pos_to_remove = list()
    pos_to_remove_trans = list()
    if order < K:
        pos_to_remove = pos_to_remove_helper(order, L_k_s, K, H_K)
        pos_to_remove_np = np.array(pos_to_remove, dtype=np.uint)
        results_dict["pos_to_remove"] = pos_to_remove_np
    # create effects table (same as alphas_table if ord == K)
    effects_table = np.copy(alpha)
    if order < K:
        effects_table = np.delete(effects_table, pos_to_remove, 0)
    results_dict["effects_table"] = effects_table
    return results_dict


def calculate_M_j_s(scenario_dict):
    per_effect_M_j_s = scenario_dict["per_effect_M_j_s"]
    z = 0
    l = 0
    q = 5
    total_num_of_M_j_s = q * len(per_effect_M_j_s)
    M_j_s = np.empty(total_num_of_M_j_s, dtype='uint32')
    for x in per_effect_M_j_s:
        while l < q:
            M_j_s[z * q + l] = x
            l = l + 1
        l = 0
        z = z + 1
    M_j_s = np.reshape(M_j_s, (total_num_of_M_j_s, 1))
    return M_j_s


## begin permutation related 

def build_table_of_perms_of_dims(perm_mat_row_one_lists, K):
    list_of_perms_of_dims = list()
    perm_mat_row_two_lists = list()
    for x in perm_mat_row_one_lists:
        perm_mat_row_two_lists.append(list(itertools.permutations(x)))
    for x in itertools.product(*perm_mat_row_two_lists):
        my_perm = np.empty(K)
        for i, s in enumerate(perm_mat_row_one_lists):
            my_perm[s] = x[i]
        list_of_perms_of_dims.append(my_perm)
    table_of_perms_of_dims = np.array(list_of_perms_of_dims, dtype='uint')
    return table_of_perms_of_dims

def find_perm_of_class_numbers(alpha, perm_of_dim):
    B = alpha[:, perm_of_dim]
    perm_of_class_numbers = np.lexsort(np.rot90(B))
    return perm_of_class_numbers

def find_pos_perms_and_inverse_perms(results_dict,
        order, table_perms_of_dims, H, K, table_perms_of_class_nums,
        pos_to_remove, pos_mapping_dict, alpha,
        meas = True):
    # create table of perms of pos nums and inverse perms of pos nums
    num_of_perms = np.shape(table_perms_of_dims)[0]
    table_perms_of_pos_nums = np.empty((num_of_perms, H), dtype='uint')
    table_inverse_perms_of_pos_nums = np.empty_like(table_perms_of_pos_nums)
    for i, perm_of_dim in enumerate(table_perms_of_dims):
        perm_of_class_nums = table_perms_of_class_nums[i, :]
        perm_of_pos_nums = np.copy(perm_of_class_nums)
        if order < K: # perform necessary deletions and renumberings
            perm_of_pos_nums = np.delete(perm_of_pos_nums, pos_to_remove, 0)
            perm_of_pos_nums = np.vectorize(pos_mapping_dict.get)(
                perm_of_pos_nums)
            inverse_perm_of_pos_nums = np.argsort(perm_of_pos_nums)
            table_perms_of_pos_nums[i, :] = perm_of_pos_nums
            table_inverse_perms_of_pos_nums[i, :] = inverse_perm_of_pos_nums
        else: # simply make a copy
            inverse_perm_of_pos_nums = np.argsort(perm_of_pos_nums)
            table_perms_of_pos_nums[i, :] = perm_of_pos_nums
            table_inverse_perms_of_pos_nums[i, :] = inverse_perm_of_pos_nums
        # assign perm and inverse perm of pos nums to table
    # save pos_to_remove for use in cpp function generate_design_matrix
    # note: only save it if it's a non-empty vector, i.e. if order < K
    # save table of perms and inverse perms of class nums and pos nums
    key_name = ""
    if meas is True:
        key_name = "table_perms_of_pos_nums"
    else: 
        key_name = "table_perms_of_pos_nums_trans"
    results_dict[key_name] = table_perms_of_pos_nums
    key_name = ""
    if meas is True:
        key_name = "table_inverse_perms_of_pos_nums"
    else:
        key_name = "table_inverse_perms_of_pos_nums_trans"
    results_dict[key_name] = table_inverse_perms_of_pos_nums
    return results_dict

def create_perm_and_inverse_perm_tables(mydict, alpha, scenario_dict,
                                        fixed_other_vals_dict):
    results_dict = dict()
    lists_of_col_nums = list(mydict.values())
    K = scenario_dict["K"]
    T = fixed_other_vals_dict["T"]
    table_perms_of_dims = build_table_of_perms_of_dims(lists_of_col_nums, K)
    num_of_perms = np.shape(table_perms_of_dims)[0]
    ### first find and save tables of perms and inverse perms of dims
    table_inverse_perms_of_dims = np.empty_like(table_perms_of_dims)
    for i, perm in enumerate(table_perms_of_dims):
        # find inverse permutation of each perm
        table_inverse_perms_of_dims[i, :] = np.argsort(perm)
    results_dict["table_perms_of_dims"] = table_perms_of_dims
    results_dict["table_inverse_perms_of_dims"] = table_inverse_perms_of_dims
    # create the tables of perms and inverse perms of class nums
    L_k_s = scenario_dict["L_k_s"]
    K = scenario_dict["K"]
    H_K = np.prod(L_k_s)
    table_perms_of_class_nums = np.empty((num_of_perms, H_K), dtype='uint')
    table_inverse_perms_of_class_nums = np.empty_like(table_perms_of_class_nums)
    for i, perm_of_dim in enumerate(table_perms_of_dims):
        perm_of_class_nums = find_perm_of_class_numbers(alpha, perm_of_dim)
        inverse_perm_of_class_nums = np.argsort(perm_of_class_nums)
        table_perms_of_class_nums[i, :] = perm_of_class_nums
        table_inverse_perms_of_class_nums[i, :] = inverse_perm_of_class_nums
    results_dict["table_perms_of_class_nums"] = table_perms_of_class_nums
    results_dict["table_inverse_perms_of_class_nums"] = \
        table_inverse_perms_of_class_nums
    ### now, deal with position number-related tasks
    #### now do the general finding and save tables of perms and inverse perms
    #### of class numbers
    ## first do order
    order = fixed_other_vals_dict["order"]
    pos_mapping_dict, pos_to_keep = build_pos_related_structures(
        order, L_k_s, K)
    # find pos_to_remove
    seq_H_K = set(range(0, H_K))
    pos_to_remove = sorted(seq_H_K.difference(set(pos_to_keep)))
    H = len(pos_to_keep)
    # create table of perms of pos nums and inverse perms of pos nums
    results_dict = find_pos_perms_and_inverse_perms(results_dict,
        order, table_perms_of_dims, H, K, table_perms_of_class_nums,
        pos_to_remove, pos_mapping_dict,
        alpha, True)
    return results_dict

## end permutation related 
