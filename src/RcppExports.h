#include <Rcpp.h>
#include <cmath>

#ifndef INCLUDED_UTILS_H
#define INCLUDED_UTILS_H

using namespace Rcpp;
using namespace std;

const double MINUS_LOG_SQRT_2_PI = -0.9189385;

inline double log_ncr(int n, int r) {
	double result = 0.0;
	for (int i = r + 1; i <= n; i++)
		result += log((double)i);
	for (int i = 1; i <= n - r; i++)
		result -= log((double)i);
	return result;
}

inline int random_integer(int exc_max)
{
	return (int)(unif_rand() * (double)exc_max) % exc_max;
}

inline double sq(double x) { return x * x; }

inline double log_likelihood_normal(
	double mean,
	double sd,
	double value
){ return MINUS_LOG_SQRT_2_PI - log(sd) - 0.5 * sq((value - mean)/sd); }

Rcpp::IntegerVector sample_int(int n, int r);

RcppExport SEXP R_es(SEXP S, SEXP p, SEXP r);

RcppExport SEXP R_gset(
	SEXP N,
	SEXP S,
	SEXP p,
	SEXP r,
	SEXP min_its,
	SEXP max_its,
	SEXP signif,
	SEXP log_dismiss
);

#endif
