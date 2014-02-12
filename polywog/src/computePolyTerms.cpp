#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List computePolyTerms(int degree,
                      int k_expand,
                      int k_lin)
{
    // The return object is a list of matrices, one for each degree up to the
    // number specified (plus an additional one if there are any terms
    // included linearly)
    List ans(k_lin > 0 ? degree + 1: degree);

    int k = 1;
    for (int d = 1; d <= degree; d++) {
        int last_k = k;
        k *= k_expand;

        IntegerMatrix poly_terms(k, k_expand + k_lin);
        IntegerMatrix last_poly_terms;
        if (d > 1)
            last_poly_terms = as<IntegerMatrix>(ans[d-2]);

        // What we want to produce:
        //
        //   d = 1 -> an identity matrix (augmented to the right with k_lin
        //            columns of zeros)
        //            
        //   d > 1 -> stack the matrix for the previous degree k_expand times,
        //            each time adding a column of ones to the corresponding
        //            term
        for (int i = 0; i < k; i++) {
            if (d == 1) {
                // Add one along the diagonal
                poly_terms(i, i) = 1;
            } else {
                // How many times have we run through the previous matrix?
                // That's the index of the term we'll be adding to
                int term_to_add = i / last_k;
                int row_of_last = i % last_k;

                // For degree > 1, stack the previous matrix k_expand times,
                // each time incrementing the corresponding term by 1
                poly_terms(i, _) = last_poly_terms(row_of_last, _);
                poly_terms(i, term_to_add) += 1;
            }
        }

        ans[d-1] = poly_terms;
    }

    // Compute the matrix for the terms included linearly, if there are any.
    // This is an identity matrix of size k_lin, augmented with k_expand
    // columns of zeros to the left
    if (k_lin > 0) {
        IntegerMatrix linear_terms(k_lin, k_expand + k_lin);
        for (int i = 0; i < k_lin; i++)
            linear_terms(i, k_expand + i) = 1;

        ans[degree] = linear_terms;
    }

    return ans;
}
