#include <Rcpp.h>

using namespace Rcpp;

// Compute the polynomial expansion of a vector, according to a 'poly_terms'
// object as returned by polym2()
NumericVector rawToPoly(NumericVector x, IntegerMatrix poly_terms) {
    int n_poly = poly_terms.nrow();
    int n_variables = poly_terms.ncol();

    if (x.size() != n_variables)
        stop("'x' must be the same length as the number of columns in 'poly_terms'");

    NumericVector ans(n_poly + 1, 1.0); // Initialize each as 1 since we're
                                        // multiplying within the loop
    IntegerVector powers;

    // Outside loop: Each row of the polyTerms matrix
    for (int i = 0; i < n_poly; i++) {
        // Inside loop: Each element of the row, which is the power of the
        // corresponding variable of the x vector
        powers = poly_terms.row(i);
        for (int j = 0; j < n_variables; j++) {
            if (powers[j] > 0)
                ans[i+1] *= pow(x[j], powers[j]);
        }
    }

    return ans;
}

