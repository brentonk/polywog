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

// [[Rcpp::export]]
List newPredictPolywogC(NumericMatrix X,
                        IntegerMatrix poly_terms,
                        List coef,
                        bool avg,
                        bool interval,
                        double level)
{
    // Storage for final output
    //
    // avg -> just a single average and its confidence interval
    // interval -> prediction and confidence interval for each observation
    // neither -> just a prediction for each observation
    int n_obs = X.nrow();
    NumericVector predicted(avg ? 1 : n_obs, 0.0);
    NumericVector low(avg ? 1 :
                      interval ? n_obs :
                      0);
    NumericVector high(avg ? 1 :
                       interval ? n_obs :
                       0);

    // Extract the main coefficients and the bootstrap matrix
    NumericVector coef_main = coef["main"];
    int n_boot = 0;
    S4 boot_matrix;
    if (avg || interval) {
        boot_matrix = as<S4>(coef["boot"]);
        IntegerVector dim_boot = boot_matrix.slot("Dim");
        n_boot = dim_boot[1];  // Bootstrap iterations = number of columns
    }
    NumericVector fit_boot(n_boot, 0.0);

    // Set up for the confidence interval calculations
    Environment stats("package:stats");
    Function quantile = stats["quantile"];
    NumericVector probs(avg || interval ? 2 : 0);
    if (avg || interval) {
        probs[0] = 0.5 - (level/2);
        probs[1] = 1 - probs[0];
    }

    // Loop over each row (observation) in X
    for (int i = 0; i < n_obs; i++) {
        // Calculate predicted values from the main model fit
        NumericVector x_row = X.row(i);
        NumericVector x_row_poly = rawToPoly(x_row, poly_terms);
        if (avg) {
            predicted[0] += sum(x_row_poly * coef_main) / n_obs;
        } else {
            predicted[i] = sum(x_row_poly * coef_main);
        }

        // Loop through the bootstrap coefficients
        if (avg || interval) {
            for (int j = 0; j < n_boot; j++) {
                NumericVector coef_boot = columnFromSparse(boot_matrix, j);

                if (avg) {
                    fit_boot[j] += sum(x_row_poly * coef_boot) / n_obs;
                } else {
                    fit_boot[j] = sum(x_row_poly * coef_boot);
                }
            }

            // The confidence interval is per-observation in the 'interval'
            // case, so compute it inside the loop
            if (!avg) {
                NumericVector conf_int = quantile(fit_boot, probs);
                low[i] = conf_int[0];
                high[i] = conf_int[1];
            }
        }
    }

    // In the 'average' case, the confidence interval is for the mean, rather
    // than for each individual prediction, so compute it outside the loop
    // over bootstrap coefficients
    if (avg) {
        NumericVector conf_int = quantile(fit_boot, probs);
        low[0] = conf_int[0];
        high[0] = conf_int[1];
    }

    if (avg || interval) {
        return List::create(_["pred"] = predicted,
                            _["lwr"] = low,
                            _["upr"] = high);
    } else {
        return List::create(_["pred"] = predicted);
    }
}
