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

import numpy as np
import pathlib
from probitlcm import _core

# for find_initial_alphas
from sklearn.decomposition import NMF
from sklearn.cluster import KMeans
import itertools

### for initialization for mcmc

def find_initial_alphas(process_dir, path_to_Ymat, T, N, K, L_k_s, seed_value):
    Ymat = _core.load_arma_umat_np(path_to_Ymat)
    Ymat = Ymat + 1
    model = NMF(n_components = K, init='random', solver ='mu',
                beta_loss='kullback-leibler',
                random_state=seed_value,
                max_iter=10000)
    W = model.fit_transform(Ymat)
    # _core.save_arma_mat_np(W, str(process_dir.joinpath("W_KL.txt")))
    H = model.components_
    # _core.save_arma_mat_np(H, str(process_dir.joinpath("H_KL.txt")))
    reshaped_data = list()
    for k in range(K):
        reshaped_data.append(W[:, k].reshape(-1, 1))
    unique_levels = list(set(L_k_s))
    # remap the L_k_s
    idx_unique_levels = list(range(len(unique_levels)))
    num_of_unique_levels = len(unique_levels)
    idx_vector = list(range(K))
    mapping_dict = dict()
    inverse_mapping_dict = dict()
    for x in idx_unique_levels:
        mapping_dict[unique_levels[x]] = x
        inverse_mapping_dict[x] = unique_levels[x]
    L_k_s_adj = [mapping_dict[x] for x in L_k_s]
    # perform K-means clustering
    performances = np.empty((K, num_of_unique_levels)) # K is n_components
    my_labels = list() # K entries, each of which is a list of len
                       #     num_of_unique_levels
    for k in range(K):
        my_labels_sublist = list()
        for idx, l in enumerate(unique_levels):
            kmeans = KMeans(n_clusters=l, init="random", n_init=10,
                            max_iter=10000, random_state=seed_value,
                            algorithm="elkan").fit(
                reshaped_data[k])
            # sort clusters in ascending order
            my_idx = np.argsort(kmeans.cluster_centers_.sum(axis=1))
            lut = np.zeros_like(my_idx)
            lut[my_idx] = np.arange(l)
            performances[k, idx] = kmeans.inertia_
            my_labels_sublist.append(lut[kmeans.labels_])
        my_labels.append(my_labels_sublist)
    permuted_dims = list(set(itertools.permutations(L_k_s_adj, K)))
    my_sums = np.empty((len(permuted_dims)))
    for idx, x in enumerate(permuted_dims):
        my_sums[idx] = np.sum(performances[idx_vector, x])
    best_fit_permuted_dims_remapped = permuted_dims[my_sums.argmin()]
    alphas = np.empty((N, K), dtype=np.uint)
    for k in range(K):
        alphas[:, k] = my_labels[k][best_fit_permuted_dims_remapped[k]]
    L_k_s_adj_np = np.array(L_k_s_adj)
    current_levels_adj = np.array(best_fit_permuted_dims_remapped)
    perms = itertools.permutations(range(current_levels_adj.size))
    desired_perm = None
    for x in perms:
        if np.array_equal(np.take(current_levels_adj, x), L_k_s_adj_np):
            desired_perm = x
            break
    # apply inverse_perm to alphas columns
    alphas = alphas[:, desired_perm]
    # save it
    process_dir = pathlib.Path(process_dir)
    _core.save_arma_umat_np(alphas,
                            str(process_dir.joinpath("initial_alphas.txt")))
