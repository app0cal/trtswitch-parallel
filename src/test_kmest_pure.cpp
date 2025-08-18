#include <Rcpp.h>
#include "pure_kmest.h"
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame test_kmest_pure(NumericVector time, IntegerVector event,
                          std::string conftype = "log‐log",
                          double conflev = 0.95,
                          bool keep_censor = false) {
  std::vector<double> t = as<std::vector<double>>(time);
  std::vector<int>    e = as<std::vector<int>>(event);
  std::vector<KMPoint> result = kmest_pure(t, e, conftype, conflev, keep_censor);

  int n = result.size();
  NumericVector r_time(n), r_nrisk(n), r_nevent(n), r_ncensor(n),
                r_surv(n), r_stddev(n), r_lower(n), r_upper(n);

  for (int i = 0; i < n; ++i) {
    r_time[i]   = result[i].time;
    r_nrisk[i]  = result[i].nrisk;
    r_nevent[i] = result[i].nevent;
    r_ncensor[i]= result[i].ncensor;
    r_surv[i]   = result[i].survival;
    r_stddev[i] = result[i].std_err;   // renamed field
    r_lower[i]  = result[i].lower;
    r_upper[i]  = result[i].upper;
  }

  return DataFrame::create(
    _["time"]    = r_time,
    _["nrisk"]   = r_nrisk,
    _["nevent"]  = r_nevent,
    _["ncensor"] = r_ncensor,
    _["survival"]= r_surv,
    _["stddev"]  = r_stddev,
    _["lower"]   = r_lower,
    _["upper"]   = r_upper
  );
}