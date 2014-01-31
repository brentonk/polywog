#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

// Reimplementation of the R function quantile(..., type = 7), the default
// type, as described in ?quantile

// [[Rcpp::export]]
NumericVector quantile7(NumericVector x, NumericVector probs) {
    // Check validity of probs
    if (is_true(any(probs < 0)) || is_true(any(probs > 1)))
        stop("'probs' must all be between 0 and 1");

    // Set up output vector and intermediate values
    int n = x.size();
    int len_probs = probs.size();
    int j;
    NumericVector y = clone(x);
    NumericVector res(len_probs);
    double lo, hi, p, m, gamma;

    // Compute the quantile corresponding to each element of probs
    for (int i = 0; i < len_probs; i++) {
        p = probs[i];
        m = 1 - p;
        j = floor(n*p + m);
        gamma = n*p + m - j;

        // Grab the j'th and (j+1)'th order statistics
        std::nth_element(y.begin(), y.begin() + j - 1, y.end());
        lo = y[j-1];
        std::nth_element(y.begin(), y.begin() + j, y.end());
        hi = y[j];

        res[i] = (1 - gamma) * lo + gamma * hi;
    }

    return res;
}
