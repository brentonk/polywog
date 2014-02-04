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

// Extract a specified column (by zero-based index) from a sparse matrix
NumericVector columnFromSparse(S4 X, int j) {
    // Set up the output as an appropriately-sized vector of 0s
    IntegerVector dim = X.slot("Dim");
    int nrow = dim[0];
    NumericVector res(nrow, 0.0);

    // The values of a "dgCMatrix" object X are stored as follows:
    //   o X@p contains the cumulative number of non-zero values by column
    //   o The indices of the non-zero values of the j'th column are contained
    //     in X@i[X@p[j]:(X@p[j+1]-1)]
    //   o The corresponding values are stored in X@x[X@p[j]:(X@p[j+1]-1)]
    IntegerVector p = X.slot("p");
    int ind_start = p[j];
    int ind_end = p[j+1];

    // Loop through the indices corresponding to the j'th column and fill in
    // the result vector appropriately
    int ind;
    IntegerVector indices = X.slot("i");
    NumericVector values = X.slot("x");
    for (int i = ind_start; i < ind_end; i++) {
        ind = indices[i];
        res[ind] = values[i];
    }

    return res;
}

