/*
 * This file is part of "probitlcm" which is released under GPL v3.
 *
 * Copyright (c) 2024 Jesse Bowers <jessemb2@illinois.edu>.
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

#include "crosscounts.hpp"

#include <armadillo>

#include <vector>

// For this subset of items, convert those columns of x to an joint variable (basis vec), currently ony support 2 items exactly
// nvalues its a vector of length J with the number of factor levels of each item
arma::uvec getBasisVec(const arma::umat& x, const arma::uvec& nvalues,
                       const arma::uvec& items) {
    arma::uvec keyvec (nvalues.n_elem, arma::fill::zeros);
    keyvec(items(0)) = 1;
    keyvec(items(1)) = nvalues(items(0));
    
    arma::uvec mytmp = x * keyvec;    
    return x * keyvec;
}

// aggregate the jointvariable into a vector of counts
arma::uvec aggregateCounts(const arma::uvec& xcol, arma::uword nvalues) {
    arma::uvec counts (nvalues, arma::fill::zeros);
  
    for (arma::uword i=0; i < xcol.n_elem; i++) {
        counts(xcol(i)) += 1;
    }
  
    return counts;
}


// get counts for every pair of items
arma::uvec aggregateAllCounts(const arma::umat& x, const arma::uvec& nvalues) {
    // std::cout << "aggregateAllCounts" << "start" << std::endl;
    unsigned int i;
    unsigned int nitems = nvalues.n_elem;
    arma::uvec sizes (nitems * (nitems-1) / 2);
    i = 0;
    for (arma::uword item0=0; item0<nitems; item0++) {
        for (arma::uword item1=item0+1; item1<nitems; item1++) {
            sizes(i) = nvalues(item0) * nvalues(item1);
            i += 1;
        }
    }
  
    arma::uword total_size = arma::accu(sizes);
    arma::uvec counts_out (total_size);
  
    unsigned int last_index = 0;
    i = 0;
    for (arma::uword item0=0; item0<nitems; item0++) {
        for (arma::uword item1=item0+1; item1<nitems; item1++) {
            counts_out.subvec(last_index,
                              last_index+sizes(i) - 1) = aggregateCounts(
                getBasisVec(x, nvalues, arma::uvec{item0, item1})
                , nvalues(item0) * nvalues(item1)
                );
            last_index += sizes(i);
            i += 1;
        }
    }
  
  
    return counts_out;
}
