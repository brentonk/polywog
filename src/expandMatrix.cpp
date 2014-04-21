#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix expandMatrix(NumericMatrix X,
                           IntegerMatrix poly_terms,
                           bool intercept = false)
{
    int n_obs = X.nrow();
    int n_raw_terms = poly_terms.ncol();
    int n_poly_terms = poly_terms.nrow();
    int n_col_out = intercept ? n_poly_terms + 1 : n_poly_terms;
    NumericMatrix ans(n_obs, n_col_out);

    // Outer loop: each polynomial term (i.e., each column of the output
    // matrix, or row of poly_terms)
    for (int i = 0; i < n_col_out; i++) {
        NumericVector thisColumn(n_obs, 1.0);

        // Inner loop: the power of each raw term represented in the given
        // polynomial term
        for (int j = 0; j < n_raw_terms; j++) {
            // Don't touch the raw terms when creating the intercept
            if (intercept && i == 0)
                break;

            int thisPower = poly_terms(intercept ? i - 1 : i, j);
            for (int d = 0; d < thisPower; d++)
                thisColumn = thisColumn * X(_, j);
        }

        ans(_, i) = thisColumn;
    }

    return ans;
}
