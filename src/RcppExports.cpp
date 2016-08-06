#include "RcppExports.h"

using namespace Rcpp;
using namespace std;

IntegerVector sample_int(int n, int r) {
	IntegerVector result(r);
	LogicalVector still_in(n);
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

double es(IntegerVector S, NumericVector r) {
	double N_r = 0.0;
	int N = r.length();
	int n = S.length();
	for (int i = 0; i < n; i++) {
		N_r += r[S[i]];
	}
	double max_hit = 0.0;
	double p_hit = 0.0;

	for (int i = 0; i < n; i++) {
		p_hit += r[S[i]]/N_r;
		double p_miss = ((double)(S[i] - i))/((double)N-(double)n);
		if ((p_hit - p_miss) > max_hit) max_hit = p_hit - p_miss;
	}

	return max_hit;
}

double gset(
	int N,
	IntegerVector S,
	NumericVector r,
	int min_its,
	int max_its,
	double signif,
	double log_dismiss
) {
	int n = S.length();
	double max_hit = es(S, r);
	int as_sim = 0;
	int samples = 0;

	do {
		IntegerVector samp = sample_int(N, n);
		sort(samp.begin(), samp.end());
		samples++;
		as_sim += (int)(es(samp, r) >= max_hit);
	}
	while (samples < min_its || ((R::pnorm((double)as_sim, (double)samples*signif, sqrt((double)samples*signif*(1.0-signif)),false,true) > log_dismiss) && (samples < max_its)));

	return (double)(as_sim + 1) / (double)(samples + 1);
}

RcppExport SEXP R_es(SEXP S, SEXP r) {
	Rcpp::RObject __result;
	__result = wrap<double>(es(as<IntegerVector>(S), as<NumericVector>(r)));
	return __result;
}

RcppExport SEXP R_gset(
	SEXP N,
	SEXP S,
	SEXP r,
	SEXP min_its,
	SEXP max_its,
	SEXP signif,
	SEXP log_dismiss
) {
	RNGScope scope;
	Rcpp::RObject __result;

	__result = wrap<double>(gset(
		as<int>(N),
		as<IntegerVector>(S),
		as<NumericVector>(r),
		as<int>(min_its),
		as<int>(max_its),
		as<double>(signif),
		as<double>(log_dismiss)
	));
	return __result;
}
