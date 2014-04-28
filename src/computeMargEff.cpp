#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector computeMargEff(NumericMatrix X,
                             IntegerMatrix poly_terms,
                             NumericVector coef,
                             LogicalVector coef_is_zero,
                             int xvar_col)
{
    // NOTE: Make sure the R wrapper decrements xvar_col by 1 when passing it,
    // in line with C++-style indexing

    int n_coef = poly_terms.nrow();
    int n_terms = poly_terms.ncol();
    int n_obs = X.nrow();
    NumericVector ans(n_obs, 0.0);
    NumericMatrix::Column x = X(_, xvar_col);

    // Outer loop: each estimated coefficient
    for (int i = 0; i < n_coef; i++) {
        int d_xvar = poly_terms(i, xvar_col);

        // Move on if the corresponding term of the model isn't a function of
        // the variable in question, or if the coefficient is 0
        if (coef_is_zero[i + 1] || d_xvar == 0)
            continue;

        // The portion of the effect corresponding to this variable
        NumericVector effect = coef[i + 1] * d_xvar * pow(x, d_xvar - 1);

        // Inner loop: multiply by any interactions
        for (int j = 0; j < n_terms; j++) {
            if (j == xvar_col) {
                continue;
            } else if (int d_j = poly_terms(i, j)) {
                effect = effect * pow(X(_, j), d_j);
            }
        }

        ans += effect;
    }

    return ans;
}
