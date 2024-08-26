/*
 * This file, originally downloaded from
 * http://www.wilmott.com/messageview.cfm?catid=10&threadid=38771
 * was posted there without attribution. The content is mostly
 * from the R project's src/nmath/qnorm.c file.
 * That file was licensed under the GNU General Public License v2.
 * The copyrights for the file are as follows:
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2000--2023 The R Core Team
 *  Copyright (C) 1998       Ross Ihaka
 *  based on AS 241 (C) 1988 Royal Statistical Society
 */

/*
 *     Compute the quantile function for the normal distribution.
 *
 *     For small to moderate probabilities, algorithm referenced
 *     below is used to obtain an initial approximation which is
 *     polished with a final Newton step.
 *
 *     For very large arguments, an algorithm of Wichura is used.
 *
 *  REFERENCE
 *
 *     Beasley, J. D. and S. G. Springer (1977).
 *     Algorithm AS 111: The percentage points of the normal distribution,
 *     Applied Statistics, 26, 118-121.
 *
 *      Wichura, M.J. (1988).
 *      Algorithm AS 241: The Percentage Points of the Normal Distribution.
 *      Applied Statistics, 37, 477-484.
 */
#include <stdio.h>
#include <math.h>
#include "qnorm.hpp"

#define R_D_Cval(p) (lower_tail ? (1 - (p)) : (p))  /*  1 - p */

#define R_DT_CIv(p) (log_p ? (lower_tail ? -expm1(p) : exp(p)) \
                             : R_D_Cval(p))
#define R_D_Lval(p)   (lower_tail ? (p) : (1 - (p)))
#define ML_NEGINF   ((-1.0) / 0.0)
#define R_D__0  (log_p ? ML_NEGINF : 0.)        /* 0 */
#define R_D__1  (log_p ? 0. : 1.)           /* 1 */
#define ML_POSINF       (1.0 / 0.0)

#define R_DT_0  (lower_tail ? R_D__0 : R_D__1)      /* 0 */
#define R_DT_1  (lower_tail ? R_D__1 : R_D__0)      /* 1 */
#define R_Q_P01_check(p)            \
         if ((log_p  && p > 0) ||            \
                        (!log_p && (p < 0 || p > 1)) )      \
    return 0;
#define R_DT_qIv(p) (log_p ? (lower_tail ? exp(p) : - expm1(p)) \
                             : R_D_Lval(p))

#define DBL_EPSILON 0.0000001

double expm1(double x)
{
    double y, a = fabs(x);

    if (a < DBL_EPSILON) return x;
    if (a > 0.697) return exp(x) - 1;  /* negligible cancellation */

    if (a > 1e-8)
    y = exp(x) - 1;
    else /* Taylor expansion, more accurate in this range */
    y = (x / 2 + 1) * x;

    /* Newton step for solving   log(1 + y) = x   for y : */
    /* WARNING: does not work for y ~ -1: bug in 1.5.0 */
    y -= (1 + y) * (log1p (y) - x);
    return y;
}

                                                          
double qnorm(double p, double mu, double sigma)
{
    double p_, q, r, val;
     int lower_tail = 1;
     int log_p = 0;
    if (p == R_DT_0)     return ML_NEGINF;
    if (p == R_DT_1)     return ML_POSINF;
    R_Q_P01_check(p);

    if(sigma  < 0)     return 0;
    if(sigma == 0)     return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
     val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

     /* r = min(p, 1-p) < 0.075 */
     if (q > 0)
         r = R_DT_CIv(p);/* 1-p */
     else
         r = p_;/* = R_DT_Iv(p) ^=  p */

     r = sqrt(- ((log_p &&
               ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
              p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

     if(q < 0.0)
         val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

