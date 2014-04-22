#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix computeExpandMatrix(NumericMatrix X,
                                  IntegerMatrix poly_terms)
{
    int n_obs = X.nrow();
    int n_raw_terms = poly_terms.ncol();
    int n_poly_terms = poly_terms.nrow();
    NumericMatrix ans(n_obs, n_poly_terms);

    // Outer loop: each polynomial term (i.e., each column of the output
    // matrix, or row of poly_terms)
    for (int i = 0; i < n_poly_terms; i++) {
        NumericVector thisColumn(n_obs, 1.0);

        // Inner loop: the power of each raw term represented in the given
        // polynomial term
        for (int j = 0; j < n_raw_terms; j++) {
            int thisPower = poly_terms(i, j);
            for (int d = 0; d < thisPower; d++)
                thisColumn = thisColumn * X(_, j);
        }

        ans(_, i) = thisColumn;
    }

    return ans;
}
