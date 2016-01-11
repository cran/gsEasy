#include "RcppExports.h"

using namespace Rcpp;
using namespace std;

Rcpp::IntegerVector sample_int(int n, int r) {
	Rcpp::IntegerVector result(r);
	Rcpp::LogicalVector still_in(n);
	for (int i = 0; i < n; i++)
		still_in[i] = true;

	for (int i = 0; i < r; i++) {
		do {
			result[i] = random_integer(n);
		}
		while (!still_in[result[i]]);
		still_in[result[i]] = false;
	}

	return result;
}

double es(IntegerVector S, double p, NumericVector r) {
	double N_r = 0.0;
	int N = r.length();
	int n = S.length();
	for (int i = 0; i < n; i++) {
		N_r += pow(r[S[i]], p);
	}
	double max_hit = 0.0;
	double p_hit = 0.0;

	for (int i = 0; i < n; i++) {
		p_hit += pow(r[S[i]], p)/N_r;
		double p_miss = ((double)(S[i] - i))/((double)N-(double)n);
		if ((p_hit - p_miss) > max_hit) max_hit = p_hit - p_miss;
	}

	return max_hit;
}

double gset(
	int N,
	IntegerVector S,
	double p,
	NumericVector r,
	int min_its,
	int max_its,
	double signif,
	double log_dismiss
) {
	Rcpp::RNGScope scope;

	int n = S.length();
	double max_hit = es(S, p, r);

	int as_sim = 0;
	for (int i = 1; i <= (max_its-1); i++) {
		//generate sample
		Rcpp::IntegerVector samp = sample_int(N, n);
		sort(samp.begin(), samp.end());
		
		//calculate es
		if (es(samp, p, r) >= max_hit) {
			as_sim += 1;
		}

		if (i >= min_its) {
			double sd = sqrt((double)i * signif * (1.0 - signif));
			if (log_dismiss > R::pnorm((double)as_sim, (double)i * signif, sd, 0, 1))
				return (double)(as_sim + 1) / (double)(i + 1);
		}
	}

	return (double)(as_sim + 1) / (double)max_its;
}

RcppExport SEXP R_es(SEXP S, SEXP p, SEXP r) {
	return wrap<double>(es(as<IntegerVector>(S), as<double>(p), as<NumericVector>(r)));
}

RcppExport SEXP R_gset(
	SEXP N,
	SEXP S,
	SEXP p,
	SEXP r,
	SEXP min_its,
	SEXP max_its,
	SEXP signif,
	SEXP log_dismiss
) {
	return(wrap<double>(gset(
		as<int>(N),
		as<IntegerVector>(S),
		as<double>(p),
		as<NumericVector>(r),
		as<int>(min_its),
		as<int>(max_its),
		as<double>(signif),
		as<double>(log_dismiss)
	)));
}
