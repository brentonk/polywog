// Utility functions for working with sparse matrices

#include <Rcpp.h>

// Extract a specified column (by zero-based index) from a sparse matrix
Rcpp::NumericVector columnFromSparse(Rcpp::S4 X, int j) {
    // Set up the output as an appropriately-sized vector of 0s
    Rcpp::IntegerVector dim = X.slot("Dim");
    int nrow = dim[0];
    Rcpp::NumericVector res(nrow, 0.0);

    // The values of a "dgCMatrix" object X are stored as follows:
    //   o X@p contains the cumulative number of non-zero values by column
    //   o The indices of the non-zero values of the j'th column are contained
    //     in X@i[X@p[j]:(X@p[j+1]-1)]
    //   o The corresponding values are stored in X@x[X@p[j]:(X@p[j+1]-1)]
    Rcpp::IntegerVector p = X.slot("p");
    int ind_start = p[j];
    int ind_end = p[j+1];

    // Loop through the indices corresponding to the j'th column and fill in
    // the result vector appropriately
    int ind;
    Rcpp::IntegerVector indices = X.slot("i");
    Rcpp::NumericVector values = X.slot("x");
    for (int i = ind_start; i < ind_end; i++) {
        ind = indices[i];
        res[ind] = values[i];
    }

    return res;
}
