#include <Rcpp.h>
#include "survival_analysis.h"
#include "survival_analysis_pure.h"
#include "utilities.h"
#include "utilities_pure.h"
using namespace Rcpp;

// [[Rcpp::export]]
// ... keep your includes and using namespace Rcpp;

// [[Rcpp::export]]
List debug_aft_llik_score_info(DataFrame data1,
                               StringVector covariates_aft,
                               std::string dist,
                               double alpha = 0.05) {
  // 1) Fit the reference model
  NumericVector init(1, NA_REAL);
  List fit_ref = liferegcpp(
    data1, "", "", "pps", "", "event",
    covariates_aft, "", "", "", dist, init,
    /*robust=*/false, /*plci=*/false, alpha, 50, 1.0e-9
  );

  DataFrame parest_ref = as<DataFrame>(fit_ref["parest"]);
  CharacterVector par_names = parest_ref["par"];
  NumericVector   par_beta  = parest_ref["beta"];
  double llik_ref          = as<double>(fit_ref["loglik"]);

  // 2) Build pure-CPP parameters
  int n = data1.nrow();
  std::vector<double> tstart(n), tstop(n), weight(n, 1.0), offset(n, 0.0);
  std::vector<int> status(n), stratum(n, 1); // single stratum

  NumericVector timev  = data1["pps"];
  IntegerVector eventv = data1["event"];
  for (int i = 0; i < n; ++i) {
    tstart[i] = timev[i];
    if (eventv[i] == 1) {
      tstop[i]  = tstart[i];
      status[i] = 1;
    } else {
      tstop[i]  = NA_REAL;   // right-censored
      status[i] = 0;
    }
  }

  // IMPORTANT: add intercept column (col 0 = 1.0)
  int nvar_raw = covariates_aft.size();   // swtrt + base2 columns you passed
  int nvar     = 1 + nvar_raw;            // +1 for intercept
  std::vector<std::vector<double>> z(n, std::vector<double>(nvar, 1.0));
  for (int j = 0; j < nvar_raw; ++j) {
    Rcpp::String nm = covariates_aft[j];
    NumericVector col = data1[nm];
    for (int i = 0; i < n; ++i) z[i][j+1] = col[i];  // shift by 1
  }

  aftparams_pure param;
  param.dist    = dist;
  param.strata  = stratum;
  param.tstart  = tstart;
  param.tstop   = tstop;
  param.status  = status;
  param.weight  = weight;
  param.offset  = offset;
  param.z       = z;
  param.nstrata = 1;

  // 3) Build beta vector in the SAME order expected by pure code:
  //    [ (Intercept), covariates_aft..., Log(scale) ]
  int p = nvar + (dist == "exponential" ? 0 : param.nstrata);
  std::vector<double> beta_ref_full(p, std::numeric_limits<double>::quiet_NaN());

  auto get_by_name = [&](const std::string& key)->double {
    for (int i = 0; i < par_names.size(); ++i) {
      if (key == std::string(par_names[i])) return par_beta[i];
    }
    return std::numeric_limits<double>::quiet_NaN();
  };

  // Intercept
  beta_ref_full[0] = get_by_name("(Intercept)");

  // Covariates in the same order you constructed z
  for (int j = 0; j < nvar_raw; ++j) {
    beta_ref_full[1 + j] = get_by_name(std::string(covariates_aft[j]));
  }

  // Log(scale) (single stratum)
  if (dist != "exponential") {
    double logsig = get_by_name("Log(scale)");
    if (!R_finite(logsig)) logsig = get_by_name("Log(scale 1)");
    beta_ref_full[nvar] = logsig;
  }

  // 4) Compare llik/score/info at the SAME parameter vector
  double llik_pure = f_llik_1_cpp(p, beta_ref_full, &param);
  std::vector<double> score_pure = f_score_1_cpp(p, beta_ref_full, &param);
  std::vector<std::vector<double>> info_pure = f_info_1_cpp(p, beta_ref_full, &param);

  return List::create(
    _["llik_ref"]   = llik_ref,
    _["llik_pure"]  = llik_pure,
    _["beta_ref"]   = beta_ref_full,
    _["score_pure"] = score_pure,
    _["info_pure"]  = info_pure,
    _["fit_ref"]    = fit_ref
  );
}
