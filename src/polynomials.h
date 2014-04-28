#ifndef POLYWOG_POLYNOMIALS_H
#define POLYWOG_POLYNOMIALS_H

#include <Rcpp.h>

Rcpp::NumericVector rawToPoly(Rcpp::NumericVector x,
                              Rcpp::IntegerMatrix poly_terms,
                              bool intercept = true);

#endif  // POLYWOG_POLYNOMIALS_H
