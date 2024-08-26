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

#include <cmath>

double norm_cdf(double x) {
    return (1 + std::erf(x / std::sqrt(2))) / 2;
}

double Q_denom(double x) {
    double a = 0.339;
    double b = 5.510;
    double result = (1 - a) * x + a * std::sqrt(std::pow(x, 2) + b);
    return result;
}

double norm_pdf(double x) {
    return (1.0 / std::sqrt(2.0 * arma::datum::pi)) *
        std::exp(-1.0 * std::pow(x, 2) / 2.0);
}

// this is exact
double log_of_norm_pdf(double x) {
    return (-1.0 * std::log(2 * arma::datum::pi) / 2.0) -
        (std::pow(x, 2) / 2.0);
}

// this is an approximation
double log_of_norm_cdf(double x) {
    double thresh = 3.0;
    double lower_thresh = -1.0 * thresh;
    double out = 0.0;
    if (x >= lower_thresh && x <= thresh) {
        out = std::log(norm_cdf(x));
    } else if (x > thresh) {
        out = std::log1p(-1.0 * norm_pdf(x) / Q_denom(x));
    } else if (x < lower_thresh) {
        out = log_of_norm_pdf(-1.0 * x) - std::log(Q_denom(-1.0 * x));
    }
    return out;
}
