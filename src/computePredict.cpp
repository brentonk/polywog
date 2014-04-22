#include <Rcpp.h>
#include "logit.h"
#include "polynomials.h"
#include "sparse.h"

using namespace Rcpp;

// Wrapper for R quantile()
NumericVector computeConfInt(NumericVector x,
                             double level)
{
    Environment stats("package:stats");
    Function quantile = stats["quantile"];
    NumericVector probs = NumericVector::create(0.5 - (level/2),
                                                0.5 + (level/2));
    return quantile(x, probs);
}

// [[Rcpp::export]]
List computePredict(NumericMatrix X,
                    IntegerMatrix poly_terms,
                    List coef,
                    bool forPredVals,
                    bool interval,
                    bool bag,
                    double level,
                    bool transform = false)
{
    // Storage for final output
    //
    // interval -> prediction and confidence interval for each observation
    // forPredVals -> just a single average (and its confidence interval if
    //                requested)
    // neither -> just a prediction for each observation
    int n_obs = X.nrow();
    int n_pred, n_interval;
    if (forPredVals) {
        n_pred = 1;
        n_interval = interval ? 1 : 0;
    } else {
        n_pred = n_obs;
        n_interval = interval ? n_obs : 0;
    }
    NumericVector predicted(n_pred);
    NumericVector low(n_interval), high(n_interval);

    // Extract the main coefficients and the bootstrap matrix
    NumericVector coef_main = coef["main"];
    int n_boot = 0;
    S4 boot_matrix;
    if (interval || bag) {
        boot_matrix = as<S4>(coef["boot"]);
        IntegerVector dim_boot = boot_matrix.slot("Dim");
        n_boot = dim_boot[1];  // Bootstrap iterations = number of columns
    }
    NumericVector fit_boot(interval ? n_boot : 0);

    // Loop over each row (observation) in X
    for (int i = 0; i < n_obs; i++) {
        // Calculate predicted values from the main model fit
        NumericVector x_row = X.row(i);
        NumericVector x_row_poly = rawToPoly(x_row, poly_terms);
        double xb_row = sum(x_row_poly * coef_main);
        if (transform)
            xb_row = transformLogit(xb_row);

        // Store result from main fit (unless bagging)
        if (forPredVals) {
            predicted[0] += xb_row / n_obs;
        } else if (!bag) {
            predicted[i] = xb_row;
        }

        // Loop through the bootstrap coefficients (don't need to wrap this
        // inside a conditional because n_boot is 0 if the bootstrap matrix
        // isn't being accessed)
        for (int j = 0; j < n_boot; j++) {
            NumericVector coef_boot = columnFromSparse(boot_matrix, j);
            double xb_boot = sum(x_row_poly * coef_boot);
            if (transform)
                xb_boot = transformLogit(xb_boot);

            if (forPredVals) {
                // Under predVals, to calculate a confidence interval, we are
                // averaging across *observations* for each *bootstrap
                // coefficient* (and then taking the order statistics of the
                // resulting averages)
                fit_boot[j] += xb_boot / n_obs;
            } else {
                // Under bootstrap aggregation, to calculate the final
                // predicted value for each *observation*, we are averaging
                // across the predicted value given according to each
                // *bootstrap coefficient*
                if (bag)
                    predicted[i] += xb_boot / n_boot;

                // In the normal case (i.e., not predVals), to calculate the
                // confidence interval for each observation, we are using the
                // order statistics of the set of predicted values given each
                // bootstrap coefficient
                if (interval)
                    fit_boot[j] = xb_boot;
            }
        }

        // The confidence interval is per-observation in the normal case, so
        // compute it inside the loop over observations after the inner loop
        // over bootstrap coefficients is finished
        if (!forPredVals && interval) {
            NumericVector conf_int = computeConfInt(fit_boot, level);
            low[i] = conf_int[0];
            high[i] = conf_int[1];
        }
    }

    // Under predVals, the confidence interval is for the overall mean, rather
    // than for each individual prediction, so compute it outside the loop
    // over observations
    if (forPredVals && interval) {
        NumericVector conf_int = computeConfInt(fit_boot, level);
        low[0] = conf_int[0];
        high[0] = conf_int[1];
    }

    if (interval) {
        return List::create(_["fit"] = predicted,
                            _["lwr"] = low,
                            _["upr"] = high);
    } else {
        return List::create(_["fit"] = predicted);
    }
}
