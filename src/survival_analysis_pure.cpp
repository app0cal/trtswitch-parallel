#include <Rcpp.h>
#include <R_ext/Applic.h>
//#include "utilities.h"
//#include "survival_analysis.h"

using namespace Rcpp;

//#if defined(_OPENMP)
//# include <omp.h>
//#endif

//custom includes (basically just utilities and survival analysis in pure cpp)
#include "survival_analysis_pure.h"
#include "utilities_pure.h"
//using namespace Rcpp;

// helper functions for converting Rcpp types to std::vector
// a few are also just rewritten rcpp functions to use std::vector or stl structs instead
/*
template<typename STL, typename RCPPTYPE>
STL to_std(const RCPPTYPE& rvec) {
  return Rcpp::as<STL>(rvec);
}
*/

// -----------------------------------------------------------------------------
// Pure-C++ equivalent of Rcpp::bygroup(data, factor):
//   - factor: length-n vector of group labels (e.g. rep or stratum)
//   - returns GroupIndex with index[] & lookup[]
// -----------------------------------------------------------------------------
group_index group_by(const std::vector<std::string>& factor) {
  int n = static_cast<int>(factor.size());
  group_index res;
  res.index.resize(n);

  // 1) assign each unique string a small integer ID
  std::unordered_map<std::string,int> map_id;
  map_id.reserve(n);
  int next_id = 1;

  for (int i = 0; i < n; ++i) {
    const auto &lvl = factor[i];
    auto it = map_id.find(lvl);
    if (it == map_id.end()) {
      map_id[lvl]    = next_id;
      res.index[i]   = next_id;
      ++next_id;
    }
    else {
      res.index[i]   = it->second;
    }
  }

  // 2) build lookup: one vector of row-indices per group
  res.lookup.resize(next_id);  // we’ll ignore slot 0
  for (int i = 0; i < n; ++i) {
    int g = res.index[i];
    res.lookup[g].push_back(i);
  }

  return res;
}

//this has variable functoin will depend on each function that uses it, so it might not be the same as it 
// would for liferegrcpp helper function 
/* hard coding this was a very bad idea do NOT do this
bool has_variable(const trial_data& d, const std::string& var_name) {
  if (var_name == "time")     return !d.pps.empty();
  if (var_name == "event")     return true;
  if (var_name == "swtrt")     return true;
  if (var_name == "aft")       return true;  // or check individual cov names elsewhere

  for (const auto& name : d.aft_names) {
    if (var_name == name) return true;
  }
  return false;
}
*/

int find_aft_index(const trial_data& d, const std::string& col_name) {
  for (size_t i = 0; i < d.aft_names.size(); ++i) {
    if (d.aft_names[i] == col_name) {
      return static_cast<int>(i);
    }
  }
  return -1; // not found
}

bool has_col_aft(const trial_data& d, const std::string& col_name) {
  if (col_name == "pps")     return !d.pps.empty();
  if (col_name == "event")    return !d.event.empty() && d.event.size() == d.pps.size();
  if (col_name == "swtrt")    return !d.swtrt.empty() && d.swtrt.size() == d.pps.size();
  
  //col name check
  int y = find_aft_index(d, col_name);
  if (y != -1 && d.aft[y].size() == d.pps.size()) {
    return true;
  }
  return false;
}

void safety_check_col_aft(const trial_data& d) {
  if (!has_col_aft(d, "pps"))     throw std::runtime_error("liferegcpp: Missing 'pps' column in trial_data");
  if (!has_col_aft(d, "event"))   throw std::runtime_error("liferegcpp: Missing 'event' column in trial_data");
  if (!has_col_aft(d, "swtrt"))   throw std::runtime_error("liferegcpp: Missing 'swtrt' column in trial_data");
  if (d.aft_names.size() != d.aft.size()) {
    throw std::runtime_error("liferegcpp: Mismatch between aft_names and aft data size");
  }
  for (size_t j = 0; j < d.aft_names.size(); ++j) {
    if (d.aft[j].size() != d.pps.size()) {
      throw std::runtime_error("liferegcpp: Mismatch between aft column size and pps size for column: " + d.aft_names[j]);
    }
  }
}


inline const std::vector<double>& get_numeric_column(const trial_data& d, const std::string& col_name) {
  if(col_name == "pps"){
    return d.pps;
  }

  for (size_t i = 0; i < d.aft_names.size(); ++i) {
    if (col_name == d.aft_names[i]) {
      return d.aft[i];
    }
  }

  throw std::runtime_error("Variable not found: " + col_name);
}

/*
//
//
//
*/

// score vector
std::vector<double> f_score_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param->z.size();
  int nvar = param->z[0].size();
  int person, i, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z[person][i];
    }
  }

  std::vector<double> sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  std::vector<double> score(p);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    //std::vector<double> z = (param->z[person] / sigma); //turned into the bottom 2 lines
    std::vector<double> z = param->z[person];      // copy the row
    for (auto &v : z) v /= sigma;     
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - std::exp(u));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -wt*(1 - std::exp(u));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*((1 - std::exp(u))*(-u) - 1);
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(u*u - 1);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c0 = 1 - 2*plogis_cpp(u,true);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(c0*u - 1);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*u;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(u*u - 1);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c0 = 1 - 2*plogis_cpp(u,true);
        double c1 = wt*c0;
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(c0*u - 1);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(std::exp(v - std::exp(v)) - std::exp(u - std::exp(u)))/
          (std::exp(-std::exp(v)) - std::exp(-std::exp(u)));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(std::exp(v - std::exp(v)) - std::exp(u - std::exp(u)))/
          (std::exp(-std::exp(v)) - std::exp(-std::exp(u)));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(std::exp(v - std::exp(v))*v - std::exp(u - std::exp(u))*u)/
          (std::exp(-std::exp(v)) - std::exp(-std::exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dnorm_cpp(v, false), d2 = dnorm_cpp(u, false);
        double q1 = pnorm_cpp(v, false), q2 = pnorm_cpp(u, false);
        double c1 = wt*(d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dlogis_cpp(v, false), d2 = dlogis_cpp(u, false);
        double q1 = plogis_cpp(v, false), q2 = plogis_cpp(u, false);
        double c1 = wt*(d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dnorm_cpp(v, false), d2 = dnorm_cpp(u, false);
        double q1 = pnorm_cpp(v, true), q2 = pnorm_cpp(u, true);
        double c1 = wt*(d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v,true), q2 = plogis_cpp(u,false);
        double c1 = wt*(d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(q1 - q2);
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-std::exp(u - std::exp(u))/(1 - std::exp(-std::exp(u))));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-std::exp(u - std::exp(u))/(1 - std::exp(-std::exp(u))));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-dnorm_cpp(u, false)/pnorm_cpp(u, true));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*(-plogis_cpp(u, true));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-dnorm_cpp(u, false)/pnorm_cpp(u, true));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*(-plogis_cpp(u, true));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*std::exp(v);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*std::exp(v);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "lognormal") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*(dnorm_cpp(v, false)/pnorm_cpp(v, true));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*plogis_cpp(v, true);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*(dnorm_cpp(v, false)/pnorm_cpp(v, true));
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = wt*plogis_cpp(v, true);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += c1*v;
      }
    }
  }

  return score;
}

/* * f_llik_1_cpp
 * 
 * This function computes the log-likelihood for a given set of parameters
 * in an accelerated failure time (AFT) model.
 * 
 * Parameters:
 *   p - number of parameters
 *   par - vector of parameters
 *   ex - pointer to aftparams_pure struct containing data and model information
 * 
 * Returns:
 *   log-likelihood value as a double
 * 
 * Notes:
 * - This function was custom built and redesigned to avoid Rcpp API calls
 *  this comes at a cost of rare tail cases where survival times being 0 or extremely large causing issues
 *  this will be addressed in future versions
 * - Potential adding of enum class to compare easier and make code easier to follow and possibly improve performance
 *  having them compare against switch cases with int for the event ==1, 0, 3 2 
 *  and the distributions being added to an enum class
 */
// log likelihood
double f_llik_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param->z.size();
  int nvar = param->z[0].size();
  int person, i, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z[person][i];
    }
  }

  std::vector<double> sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = exp(par[k]);
    }
  }

  double loglik = 0;
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];

    //for this part introduce replacemenet for'
    // R::dnorm, R:dlogis, R::pnorm, R::plogis
    if (param->status[person] == 1) { // event
      double logsig = std::log(sigma);
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(u - std::exp(u) - logsig);
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(dnorm_cpp(u,true) - logsig);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*(dlogis_cpp(u,true) - logsig);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*(dnorm_cpp(u,true) - logsig);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*(dlogis_cpp(u,true) - logsig);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(std::exp(-std::exp(v)) - std::exp(-std::exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,true) - pnorm_cpp(u,true));
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,true) - plogis_cpp(u,true));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,false) - pnorm_cpp(u,false));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,true) - plogis_cpp(u,true));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*std::log(1.0 - std::exp(-std::exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(u,false));
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(u,false));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(u,false));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(u,false));
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*(-std::exp(v));
      } else if (param->dist == "lognormal") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,true));
      } else if (param->dist == "loglogistic") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,true));
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,true));
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,true));
      }
    }
  }

  return loglik;
}

// observed information matrix
std::vector<std::vector<double>> f_info_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param->z.size();
  int nvar = param->z[0].size();
  int person, i, j, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z[person][i];
    }
  }

  std::vector<double> sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = std::exp(par[k]);
    }
  }

  std::vector<std::vector<double>> imat(p,std::vector<double>(p,0.0));
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    std::vector<double> z = param->z[person];
    for (auto& val : z) val /= sigma;
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*std::exp(u);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*std::exp(u);
        double c2 = wt*(std::exp(u)*u - (1 - std::exp(u)));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c2 = wt*2*u;
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += wt*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = wt*2*dlogis_cpp(u, false);
        double c2 = wt*(2*dlogis_cpp(u, false)*u +
                        1 - 2*dlogis_cpp(u, true));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c2 = wt*2*u;
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += wt*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = wt*2*dlogis_cpp(u, false);
        double c2 = wt*(2*dlogis_cpp(u, false)*u +
                        1 - 2*dlogis_cpp(u, true));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double w1 = std::exp(v), w2 = std::exp(u);
        double q1 = std::exp(-w1), q2 = std::exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = wt*(std::pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(1-w1) - d2*(1-w2))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double w1 = std::exp(v), w2 = std::exp(u);
        double q1 = std::exp(-w1), q2 = std::exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = wt*(std::pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(1-w1) - d2*(1-w2))/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/std::pow(q1 - q2, 2) +
                        (d1*(1 + (1-w1)*v) - d2*(1 + (1-w2)*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += wt*(std::pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1 + (1-w1)*v)*v - d2*(1 + (1-w2)*u)*u)/(q1 - q2));
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dnorm_cpp(v, false), d2 = dnorm_cpp(u, false);
        double q1 = pnorm_cpp(v,true), q2 = pnorm_cpp(u, true);
        double c1 = wt*(std::pow((d1 - d2)/(q1 -  q2), 2) +
                        (-d1*v + d2*u)/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/std::pow(q1 - q2, 2) +
                        (d1*(1 - v*v) - d2*(1 - u*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += wt*(std::pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1 - v*v)*v - d2*(1 - u*u)*u)/(q1 - q2));
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dlogis_cpp(v, false), d2 = dlogis_cpp(u, false);
        double q1 = plogis_cpp(v, true), q2 = plogis_cpp(u, true);
        double c1 = wt*(std::pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/std::pow(q1 - q2, 2) +
                        (d1*(1+(2*q1-1)*v) - d2*(1+(2*q2-1)*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += wt*(std::pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1+(2*q1-1)*v)*v - d2*(1+(2*q2-1)*u)*u)/(q1 - q2));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), d2 = dnorm_cpp(u,false);
        double q1 = pnorm_cpp(v,false), q2 = pnorm_cpp(u,false);
        double c1 = wt*(std::pow((d1 - d2)/(q1 -  q2), 2) +
                        (-d1*v + d2*u)/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/std::pow(q1 - q2, 2) +
                        (d1*(1 - v*v) - d2*(1 - u*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += wt*(std::pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1 - v*v)*v - d2*(1 - u*u)*u)/(q1 - q2));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v, false), q2 = plogis_cpp(u, false);
        double c1 = wt*(std::pow((d1 - d2)/(q1 - q2), 2) +
                        (d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2));
        double c2 = wt*((d1 - d2)*(d1*v - d2*u)/std::pow(q1 - q2, 2) +
                        (d1*(1+(2*q1-1)*v) - d2*(1+(2*q2-1)*u))/(q1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += wt*(std::pow((d1*v - d2*u)/(q1 - q2), 2) +
          (d1*(1+(2*q1-1)*v)*v - d2*(1+(2*q2-1)*u)*u)/(q1 - q2));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double w2 = std::exp(u), q2 = std::exp(-w2), d2 = w2*q2;
        double c1 = wt*(std::pow(d2/(1 - q2), 2) - d2*(1-w2)/(1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double w2 = std::exp(u), q2 = std::exp(-w2), d2 = w2*q2;
        double c1 = wt*(std::pow(d2/(1 - q2), 2) - d2*(1-w2)/(1 - q2));
        double c2 = wt*(std::pow(d2/(1 - q2), 2)*u - d2*(1 + (1-w2)*u)/(1 - q2));
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double d2 = dnorm_cpp(u,false), p2 = pnorm_cpp(u,true);
        double c1 = wt*(std::pow(d2/p2, 2) + d2*u/p2);
        double c2 = wt*(std::pow(d2/p2, 2)*u - d2*(1 - u*u)/p2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double d2 = dlogis_cpp(u,false), q2 = plogis_cpp(u,false);
        double c1 = wt*d2;
        double c2 = wt*(d2*u - q2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double d2 = dnorm_cpp(u,false), p2 = pnorm_cpp(u,true);
        double c1 = wt*(std::pow(d2/p2, 2) + d2*u/p2);
        double c2 = wt*(std::pow(d2/p2, 2)*u - d2*(1 - u*u)/p2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double d2 = dlogis_cpp(u,false), q2 = plogis_cpp(u,false);
        double c1 = wt*d2;
        double c2 = wt*(d2*u - q2);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*std::exp(v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
      } else if (param->dist == "weibull") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = wt*std::exp(v);
        double c2 = wt*std::exp(v)*(1+v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*v;
      } else if (param->dist == "lognormal") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), q1 = pnorm_cpp(v,true);
        double c1 = wt*(std::pow(d1/q1, 2) - d1*v/q1);
        double c2 = wt*(std::pow(d1/q1, 2)*v + d1*(1 - v*v)/q1);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*v;
      } else if (param->dist == "loglogistic") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), q1 = plogis_cpp(v,false);
        double c1 = wt*d1;
        double c2 = wt*(1-q1 + d1*v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), q1 = pnorm_cpp(v,false);
        double c1 = wt*(std::pow(d1/q1, 2) - d1*v/q1);
        double c2 = wt*(std::pow(d1/q1, 2)*v + d1*(1 - v*v)/q1);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), q1 = plogis_cpp(v,false);
        double c1 = wt*d1;
        double c2 = wt*(1-q1 + d1*v);
        for (i=0; i<nvar; i++) {
          for (j=0; j<=i; j++) {
            imat[i][j] += c1*z[i]*z[j];
          }
        }
        for (j=0; j<nvar; j++) {
          imat[k][j] += c2*z[j];
        }
        imat[k][k] += c2*v;
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat[i][j] = imat[j][i];
    }
  }

  return imat;
}

// score residual matrix
std::vector<std::vector<double>> f_ressco_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param->z.size();
  int nvar = param->z[0].size();
  int person, i, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z[person][i];
    }
  }

  std::vector<double> sig(n, 1.0);
  if (param->dist != "exponential") {
    for (person = 0; person < n; person++) {
      k = param->strata[person] + nvar - 1;
      sig[person] = std::exp(par[k]);
    }
  }

  std::vector<std::vector<double>> resid(n, std::vector<double>(p));
  for (person = 0; person < n; person++) {
    double sigma = sig[person];
    std::vector<double> z = param->z[person];
    for(double &v : z) {
      v /= sigma;
    }
    k = param->strata[person] + nvar - 1;

    if (param->status[person] == 1) { // event
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -(1 - std::exp(u));
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -(1 - std::exp(u));
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = (1 - std::exp(u))*(-u) - 1;
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        for (i=0; i<nvar; i++) {
          resid[person][i] = u*z[i];
        }
        resid[person][k] = u*u - 1;
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = 1 - 2*plogis_cpp(u,false);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u - 1;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        for (i=0; i<nvar; i++) {
          resid[person][i] = u*z[i];
        }
        resid[person][k] = u*u - 1;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = 1 - 2*plogis_cpp(u,false);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u - 1;
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double w1 = std::exp(v), w2 = std::exp(u);
        double q1 = std::exp(-w1), q2 = std::exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double w1 = std::exp(v), w2 = std::exp(u);
        double q1 = std::exp(-w1), q2 = std::exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), d2 = dnorm_cpp(u,false);
        double q1 = pnorm_cpp(v,false), q2 = pnorm_cpp(u,false);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v,false), q2 = plogis_cpp(u,false);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), d2 = dnorm_cpp(u,false);
        double q1 = pnorm_cpp(v,false), q2 = pnorm_cpp(u,false);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = (d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v,false), q2 = plogis_cpp(u,false);
        double c1 = (d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = (d1*v - d2*u)/(q1 - q2);
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double w2 = std::exp(u), q2 = std::exp(-w2), d2 = w2*q2;
        double c1 = -d2/(1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double w2 = std::exp(u), q2 = std::exp(-w2), d2 = w2*q2;
        double c1 = -d2/(1 - q2);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u;
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double c1 = (-dnorm_cpp(u,false)/pnorm_cpp(u,false));
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u;
      } else if (param->dist == "loglogistic") {
        double u = (log(param->tstop[person]) - eta[person])/sigma;
        double c1 = -plogis_cpp(u,false);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u;
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = -dnorm_cpp(u,false)/pnorm_cpp(u,true);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u;
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double c1 = -plogis_cpp(u,false);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*u;
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = std::exp(v);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
      } else if (param->dist == "weibull") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = std::exp(v);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*v;
      } else if (param->dist == "lognormal") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = dnorm_cpp(v)/pnorm_cpp(v,false);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*v;
      } else if (param->dist == "loglogistic") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double c1 = plogis_cpp(v,true);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*v;
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = dnorm_cpp(v,false)/pnorm_cpp(v,false);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*v;
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        double c1 = plogis_cpp(v,true);
        for (i=0; i<nvar; i++) {
          resid[person][i] = c1*z[i];
        }
        resid[person][k] = c1*v;
      }
    }
  }

  return resid;
}

// substitute information matrix guaranteed to be positive definite
std::vector<std::vector<double>> f_jj_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param->z.size();
  int person, i, j;

  std::vector<std::vector<double>> resid = f_ressco_1_cpp(p, par, param);
  std::vector<std::vector<double>> jj(p, std::vector<double>(p));
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        jj[i][j] += wt*resid[person][i]*resid[person][j];
      }
    }
  }

  return jj;
}
/*Initial call
List liferegloop(int p, NumericVector par, void *ex,
                 int maxiter, double eps,
                 IntegerVector colfit, int ncolfit) {

*/
liferegloopresult liferegloop_cpp(int p, std::vector<double> par, void *ex,
                 int maxiter, double eps,
                 std::vector<int> colfit, int ncolfit) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int i, j, iter, halving = 0;
  bool fail;

  int nstrata = param -> nstrata;
  int nsub = static_cast<int>(param -> z.size()); // reminder that z is a vector so grabbing the .size() gives size_t not int, not a massive difference 
  int nvar = static_cast<int>(param -> z[0].size()); // important
  std::vector<std::vector<double>> z1 = param -> z;

  // standardize the design matrix
  std::vector<double> mu(nvar, 0.0), sigma(nvar, 1.0);
  std::vector<std::vector<double>> z2(nsub, std::vector<double>(nvar, 0.0));

  for (i=0; i<nvar; i++) {
    std::vector<double> u(nsub);
    for (j=0; j < nsub; j++) {
      u[j] = z1[j][i];
    }

    bool binary = std::all_of(u.begin(), u.end(), [](double val) { return val == 0 || val == 1; });

    if (binary) { // no standardization
      mu[i] = 0.0;
      sigma[i] = 1.0;
    } else {
      mu[i] = mean_cpp(u);
      sigma[i] = sd_cpp(u);
    }
    if (sigma[i] <= 0.0) sigma[i] = 1e-12;  // avoid division by zero
  }

  for(int i = 0; i < nsub; i++) {
    for (j=0; j < nvar; j++) {
      z2[i][j] = (z1[i][j] - mu[j]) / sigma[j];
    }
  }

  // corresponding initial beta
  std::vector<double> beta(p), newbeta(p);
  beta[0] = par[0];
  for (i=1; i<nvar; i++) {
    beta[i] = par[i]*sigma[i];
    beta[0] += par[i]*mu[i];
  }
  if (param->dist != "exponential") {
    for (i=nvar; i<p; i++) {
      beta[i] = par[i];
    }
  }
  
  aftparams_pure para = {param->dist, param->strata, param->tstart, param->tstop,
                    param->status, param->weight, param->offset, z2, nstrata};

  double toler = 1e-12;
  double newlk = -1.0;
  double loglik = f_llik_1_cpp(p, beta, &para);
  // assign loglik and newlk values right now to avoid garbage values 
  std::vector<double> u(p);
  std::vector<std::vector<double>> imat(p, std::vector<double>(p));
  std::vector<std::vector<double>> jj(p, std::vector<double>(p));
  std::vector<double> u1(ncolfit);
  std::vector<std::vector<double>> imat1(ncolfit, std::vector<double>(ncolfit));
  std::vector<std::vector<double>> jj1(ncolfit, std::vector<double>(ncolfit));

  //loglik = f_llik_1_cpp(p, beta, &para);
  u = f_score_1_cpp(p, beta, &para);
  for (i=0; i<ncolfit; i++) {
    u1[i] = u[colfit[i]];
  }

  imat = f_info_1_cpp(p, beta, &para);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1[i][j] = imat[colfit[i]][colfit[j]];
    }
  }

  i = cholesky2_cpp(imat1, ncolfit, toler);
  if (i < 0) {
    jj = f_jj_1_cpp(p, beta, &para);
    for (i=0; i<ncolfit; i++) {
      for (j=0; j<ncolfit; j++) {
        jj1[i][j] = jj[colfit[i]][colfit[j]];
      }
    }

    i = cholesky2_cpp(jj1, ncolfit, toler);
    chsolve2_cpp(jj1, ncolfit, u1);
  } else {
    chsolve2_cpp(imat1, ncolfit, u1);
  }

  std::fill(u.begin(), u.end(), 0.0);
  for (i=0; i<ncolfit; i++) {
    u[colfit[i]] = u1[i];
  }

  // new beta
  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + u[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    // new log likelihood
    newlk = f_llik_1_cpp(p, newbeta, &para);

    // check convergence
    // can replace fail with a single check instead of
    // fail = !std::isfinite(newlk); 
    // worth considering 
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) {
      break;
    }

    if (fail || newlk < loglik) { // adjust step size if likelihood decreases
      halving++;
      for (i=0; i<p; i++) {
        newbeta[i] = (beta[i] + newbeta[i]) /2;
      }

      // special handling of sigmas
      if (halving == 1 && param->dist != "exponential") {
        for (i=0; i<nstrata; i++) {
          if (beta[nvar+i] - newbeta[nvar+i] > 1.1) {
            newbeta[nvar+i] = beta[nvar+i] - 1.1;
          }
        }
      }
    } else { // update beta normally
      halving = 0;

      for (i=0; i<p; i++) {
        beta[i] = newbeta[i];
      }
      loglik = newlk;

      u = f_score_1_cpp(p, beta, &para);
      for (i=0; i<ncolfit; i++) {
        u1[i] = u[colfit[i]];
      }

      imat = f_info_1_cpp(p, beta, &para);
      for (i=0; i<ncolfit; i++) {
        for (j=0; j<ncolfit; j++) {
          imat1[i][j] = imat[colfit[i]][colfit[j]];
        }
      }

      i = cholesky2_cpp(imat1, ncolfit, toler);
      if (i < 0) {
        jj = f_jj_1_cpp(p, beta, &para);
        for (i=0; i<ncolfit; i++) {
          for (j=0; j<ncolfit; j++) {
            jj1[i][j] = jj[colfit[i]][colfit[j]];
          }
        }

        i = cholesky2_cpp(jj1, ncolfit, toler);
        chsolve2_cpp(jj1, ncolfit, u1);
      } else {
        chsolve2_cpp(imat1, ncolfit, u1);
      }

      std::fill(u.begin(), u.end(), 0.0);
      for (i=0; i<ncolfit; i++) {
        u[colfit[i]] = u1[i];
      }

      for (i=0; i<p; i++) {
        newbeta[i] = beta[i] + u[i];
      }
    }
  }

  if (iter == maxiter) fail = 1;

  // parameter estimates on the original scale of the design matrix
  for (i=1; i<nvar; i++) {
    newbeta[i] = newbeta[i]/sigma[i];
    newbeta[0] = newbeta[0] - newbeta[i]*mu[i];
  }

  imat = f_info_1_cpp(p, newbeta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1[i][j] = imat[colfit[i]][colfit[j]];
    }
  }

  std::vector<std::vector<double>> var1 = invsympd_cpp(imat1, ncolfit, toler);
  std::vector<std::vector<double>> var(p,std::vector<double>(p));
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      var[colfit[i]][colfit[j]] = var1[i][j];
    }
  }

  liferegloopresult res;
  res.coef = newbeta;
  res.iter = iter;
  res.var = var;
  res.loglik = newlk;
  res.fail = fail;

  return res;
  /*return List::create(
    Named("coef") = newbeta,
    Named("iter") = iter,
    Named("var") = var,
    Named("loglik") = newlk,
    Named("fail") = fail);*/
}
/*call from the original code
 List fit1 = liferegcpp(
                      trial_data, "", "", "pps", "", "event", 
                      covariates_aft, "", "", "", dist, 0, 0, alpha, 
                      50, 1.0e-9);

*/
// first and second derivatives of log likelihood with respect to eta
f_der_eta_1_result f_der_eta_1_cpp(std::vector<double> eta, std::vector<double> sig, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param->z.size();
  int person;

  std::vector<double> dg(n), ddg(n);
  for (person = 0; person < n; person++) {
    double sigma = sig[person];
    if (param->status[person] == 1) { // event
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        dg[person] = -(1 - std::acosf(u))/sigma;
        ddg[person] = -std::exp(u)/(sigma*sigma);
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        dg[person] = u/sigma;
        ddg[person] = -1/(sigma*sigma);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        dg[person] = (1 - 2*plogis_cpp(u,false))/sigma;
        ddg[person] = -2*dlogis_cpp(u,false)/(sigma*sigma);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        dg[person] = u/sigma;
        ddg[person] = -1/(sigma*sigma);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        dg[person] = (1 - 2*plogis_cpp(u,false))/sigma;
        ddg[person] = -2*dlogis_cpp(u,false)/(sigma*sigma);
      }
    } else if (param->status[person] == 3) { // interval censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double w1 = std::exp(v), w2 = std::exp(u);
        double q1 = std::exp(-w1), q2 = std::exp(-w2);
        double d1 = w1*q1, d2 = w2*q2;
        dg[person] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[person] = -(d1*(1-w1) - d2*(1-w2))/(q1 - q2)/(sigma*sigma) - 
          std::pow(dg[person], 2);
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), d2 = dnorm_cpp(u,false);
        double q1 = pnorm_cpp(v,false), q2 = pnorm_cpp(u,false);
        dg[person] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[person] = (d1*v - d2*u)/(q1 - q2)/(sigma*sigma) - 
          std::pow(dg[person], 2);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v,false), q2 = plogis_cpp(u,false);
        dg[person] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[person] = -(d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2)/(sigma*sigma) - 
          std::pow(dg[person], 2);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dnorm_cpp(v,false), d2 = dnorm_cpp(u,false);
        double q1 = pnorm_cpp(v,false), q2 = pnorm_cpp(u,false);
        dg[person] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[person] = (d1*v - d2*u)/(q1 - q2)/(sigma*sigma) - 
          std::pow(dg[person], 2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v,false), q2 = plogis_cpp(u,false);
        dg[person] = (d1 - d2)/(q1 - q2)/sigma;
        ddg[person] = -(d1*(2*q1-1) - d2*(2*q2-1))/(q1 - q2)/(sigma*sigma) - 
          std::pow(dg[person], 2);
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        dg[person] = -std::exp(u - std::exp(u))/(1 - std::exp(-std::exp(u)))/sigma;
        ddg[person] = (1 - std::exp(u) - std::exp(-std::exp(u)))*std::exp(u - std::exp(u))/
          std::pow((1 - std::exp(-std::exp(u)))*sigma, 2);
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        dg[person] = -dnorm_cpp(u,false)/pnorm_cpp(u,false)/sigma;
        ddg[person] = -u*dnorm_cpp(u,false)/
          (pnorm_cpp(u,true)*sigma*sigma) - std::pow(dg[person], 2);
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        dg[person] = -plogis_cpp(u,false)/sigma;
        ddg[person] = -dlogis_cpp(u,false)/(sigma*sigma);
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        dg[person] = -dnorm_cpp(u,false)/pnorm_cpp(u,true)/sigma;
        ddg[person] = -u*dnorm_cpp(u,false)/
          (pnorm_cpp(u,true)*sigma*sigma) - std::pow(dg[person], 2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        dg[person] = -plogis_cpp(u,false)/sigma;
        ddg[person] = -dlogis_cpp(u,false)/(sigma*sigma);
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        dg[person] = std::exp(v)/sigma;
        ddg[person] = -std::exp(v)/(sigma*sigma);
      } else if (param->dist == "lognormal") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        dg[person] = dnorm_cpp(v,false)/pnorm_cpp(v,false)/sigma;
        ddg[person] = v*dnorm_cpp(v,false)/
          (pnorm_cpp(v,false)*sigma*sigma) - std::pow(dg[person], 2);
      } else if (param->dist == "loglogistic") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        dg[person] = plogis_cpp(v,true)/sigma;
        ddg[person] = -dlogis_cpp(v,false)/(sigma*sigma);
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        dg[person] = dnorm_cpp(v,false)/pnorm_cpp(v,false)/sigma;
        ddg[person] = v*dnorm_cpp(v,false)/
          (pnorm_cpp(v,false)*sigma*sigma) - std::pow(dg[person], 2);
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        dg[person] = plogis_cpp(v,true)/sigma;
        ddg[person] = -dlogis_cpp(v,false)/(sigma*sigma);
      }
    }
  }
  f_der_eta_1_result res;
  res.dg = dg;
  res.ddg = ddg;
  /*List result = List::create(
    Named("dg") = dg,
    Named("ddg") = ddg);
  */
  return res;
}

// confidence limit of profile likelihood method
double liferegplloop_cpp(int p, std::vector<double> par, void *ex,
                     int maxiter, double eps,
                     int k, int which, double l0) {
  aftparams_pure *param = (aftparams_pure *) ex;

  int i, j, iter;
  bool fail = 0;

  double toler = 1e-12;
  std::vector<double> beta(p), newbeta(p);
  double loglik, newlk;
  std::vector<double> u(p);
  std::vector<double> delta(p);
  std::vector<std::vector<double>> imat(p, std::vector<double>(p));
  std::vector<std::vector<double>> jj(p, std::vector<double>(p));
  std::vector<std::vector<double>> v(p, std::vector<double>(p));

  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  loglik = f_llik_1_cpp(p, beta, param);
  u = f_score_1_cpp(p, beta, param);

  imat = f_info_1_cpp(p, beta, param);
  jj = imat;

  i = cholesky2_cpp(jj, p, toler);
  if (i < 0) {
    jj = f_jj_1_cpp(p, beta, param);
    v = invsympd_cpp(jj, p, toler);
  } else {
    v = invsympd_cpp(imat, p, toler);
  }
  //v = -1.0*v;
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      v[i][j] = -1.0 * v[i][j];
    }
  }

  // Lagrange multiplier method as used in SAS PROC LOGISTIC
  double w = 0;
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      w += u[i]*v[i][j]*u[j];
    }
  }

  double underroot = 2*(l0 - loglik + 0.5*w)/v[k][k];
  double lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
  u[k] += lambda;

  //delta.fill(0.0);
  for(int i=0; i<p; i++) {
    delta[i] = 0.0;
  }
  for (i=0; i<p; i++) {
    for (j=0; j<p; j++) {
      delta[i] -= v[i][j]*u[j];
    }
  }

  // update beta
  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + delta[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    newlk = f_llik_1_cpp(p, newbeta, param);

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && fabs(newlk - l0) < eps && w < eps) {
      break;
    }

    for (i=0; i<p; i++) {
      beta[i] = newbeta[i];
    }
    loglik = newlk;

    u = f_score_1_cpp(p, beta, param);

    imat = f_info_1_cpp(p, beta, param);
    jj = imat;

    i = cholesky2_cpp(jj, p, toler);
    if (i < 0) {
      jj = f_jj_1_cpp(p, beta, param);
      v = invsympd_cpp(jj, p, toler);
    } else {
      v = invsympd_cpp(imat, p, toler);
    }
    //v = -1.0*v;
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        v[i][j] = -1.0 * v[i][j];
      }
    }

    // Lagrange multiplier method as used in SAS PROC LOGISTIC
    w = 0;
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        w += u[i]*v[i][j]*u[j];
      }
    }

    underroot = 2*(l0 - loglik + 0.5*w)/v[k][k];
    lambda = underroot < 0.0 ? 0.0 : which*sqrt(underroot);
    u[k] += lambda;

    //delta.fill(0.0);
    for (i=0; i<p; i++) {
      delta[i] = 0.0;
    }
    for (i=0; i<p; i++) {
      for (j=0; j<p; j++) {
        delta[i] -= v[i][j]*u[j];
      }
    }
    // update beta
    for (i=0; i<p; i++) {
      newbeta[i] = beta[i] + delta[i];
    }
  }

  if (iter == maxiter) fail = 1;

  if (fail) {
    //stop("The algorithm in liferegplloop did not converge");
    std::string errmsg = "The algorithm in liferegplloop did not converge";
    throw std::runtime_error(errmsg);
  }

  return newbeta[k];
}


// ##################[[Rcpp::export
List lifereg_purecpp(
    //const DataFrame data expects trial_data with 
    const trial_data data,
    const std::vector<std::string>& rep = {},
    const std::vector<std::string>& stratum = {},
    const std::string time = "time",
    const std::string time2 = "",
    const std::string event = "event",
    const std::vector<std::string>& covariates = {},
    const std::string weight = "",
    const std::string offset = "",
    const std::string id = "",
    const std::string dist = "weibull",
    const std::vector<double>& init = {},
    const bool robust = 0,
    const bool plci = 0,
    const double alpha = 0.05,
    const int maxiter = 50,
    const double eps = 1.0e-9) {

  /*Notes and things to consider:
    - The function is a direct translation of the R function liferegcpp
    - The function expects a trial_data object which contains the necessary data
    - The function handles various distributions and covariates
    - Error handling is done using exceptions
    - The function returns a List similar to the R version
    - The function does NOT need to handle hasVariable checks for time2, rep, and stratum
      because the higher level wrapper will handle that logic
  
  */

  Rcpp::Rcout << "[lifereg_purecpp] line 1662" << std::endl;

  int h, i, j, k, n = data.pps.size(); //.nrows();
  //covarities is covarities_aft from the capture list!
  int nvar = static_cast<int>(covariates.size()) + 1;
  if (nvar == 2 && (covariates[0] == "" || covariates[0] == "none")) {
    nvar = 1;
  }

  std::string dist1 = dist;
  std::for_each(dist1.begin(), dist1.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  if ((dist1 == "log-logistic") || (dist1 == "llogistic")) {
    dist1 = "loglogistic";
  } else if  ((dist1 == "log-normal") || (dist1 == "lnormal")) {
    dist1 = "lognormal";
  } else if (dist1 == "gaussian") {
    dist1 = "normal";
  }

  if (!((dist1 == "exponential") || (dist1 == "weibull") ||
      (dist1 == "lognormal") || (dist1 == "loglogistic") ||
      (dist1 == "normal") || (dist1 == "logistic"))) {
    std::string str1 = "dist must be exponential, weibull, lognormal,";
    std::string str2 = "loglogistic, normal, or logistic";
    std::string errmsg = str1 + " " + str2;
    throw std::runtime_error(errmsg);
  }

  //rework this part with original logic
  /*
  //group_index out;
  bool has_rep;
  //IntegerVector repn(n);
  std::vector<int> repn(n);
  std::vector<std::vector<int>> u_rep;
  int p_rep = static_cast<int>(rep.size());
  if (p_rep == 1 && (rep[0] == "" || rep[0] == "none")) {
    has_rep = 0;
    std::fill(repn.begin(), repn.end(), 1);
  } else {
    group_index out = group_by(rep);
    has_rep = 1;
    //swaps the pointer to the index and lookup vectors
    //std::move is used to avoid copying and is faster
    repn = std::move(out.index);
    u_rep = std::move(out.lookup);
  }


  std::vector<int> stratumn(n);
  int p_stratum = static_cast<int>(stratum.size());
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    std::fill(stratumn.begin(), stratumn.end(), 1);
  } else {
    group_index out = group_by(stratum);
    stratumn = std::move(out.index);
  }
  */
  //initially we check if empty, but because of the higher level wrapper they are always called empty so just follow empty logic

  //bool has_rep = 0;
  // TEMPORARILY INITIALIZED LIKE THIS BECAUSE TSESIMP NEVER CALLS THE ELSE FUNCTIONS!
  std::vector<int> repn(n, 1);

  std::vector<int> stratumn(n, 1);


  // originally meant to be unique values only but we can do it with a simple for loop
  std::vector<int> stratumn1;
  std::unordered_set<int> seen;
  stratumn1.reserve(stratumn.size());
  for (int x : stratumn) {
    if (seen.insert(x).second) {
      // x wasn’t in the set yet, so it’s a new unique value
      stratumn1.push_back(x);
    }
  }
  int nstrata = static_cast<int>(stratumn1.size());
  int p = dist1 == "exponential" ? nvar : (nvar+nstrata);

  if (dist1 == "exponential" && nstrata > 1) {
    throw std::runtime_error("Stratification is not valid with the exponential distribution");
  }

  const std::vector<double>& timenz = data.pps;
  std::vector<double> timen = timenz; // clone(timenz); // no need to clone, we can use the vector directly
  for (i=0; i<n; i++) {
    if (!std::isnan(timen[i]) && ((dist1 == "exponential") ||
        (dist1 == "weibull") || (dist1 == "lognormal") ||
        (dist1 == "loglogistic")) && (timen[i] <= 0)) {
      std::string str1 = "time must be positive for each subject for the";
      std::string str2 = "distribution";
      std::string errmsg = str1 + " " + dist1 + " " + str2;
      throw std::runtime_error(errmsg);
    }
  }

  //high level wrapper currently level passes anything for time2
  //for now this will be skipped or left empty but functionality will be added after testing entire pipeline
  //bool has_time2 = false; //hasVariable(data, time2);
  //NumericVector time2n(n);
  std::vector<double> time2n(n, 0.0);
  if (!data.time2.empty()) {
    std::vector<double> time2nz = data.time2;
    time2n = time2nz; // clone(time2nz); // no need to clone, we can use the vector directly
    for(i=0; i<n; i++) {
      if (!std::isnan(time2n[i]) && (
          (dist1 == "exponential") ||
          (dist1 == "weibull") || (dist1 == "lognormal") ||
          (dist1 == "loglogistic")
            )  && (time2n[i] <= 0.0)) {
        std::string str1 = "time2 must be positive for each subject for the";
        std::string str2 = "distribution";
        std::string errmsg = str1 + " " + dist1 + " " + str2;
        throw std::runtime_error(errmsg);
      }
    }
  }

  
  if(data.event.empty() && data.time2.empty()) {
    throw std::runtime_error("data must contain the event variable for right censored data");
  }

  std::vector<int> eventn(n);
  if (!data.event.empty()) {
    const std::vector<int>& eventnz = data.event;
    eventn = eventnz; // clone(eventnz);
    int sum = 0;
    for(std::size_t i = 0; i < eventn.size(); ++i) {
      sum += eventn[i];
      if (eventn[i] != 1 && eventn[i] != 0) {
        throw std::runtime_error("event must be 1 or 0 for each subject");
      }
    }
    if(sum == 0){
        throw std::runtime_error("at least 1 event is needed to fit the parametric model");
    }
  }


  std::vector<std::vector<double>> zn(n, std::vector<double>(nvar, 0.0));
  for (int i=0; i<n; i++) {
    zn[i][0] = 1; // intercept
  }

  for (j=0; j < nvar - 1; j++) {
    std::string zj = covariates[j];
    if (!has_col_aft(data, zj)) {
      throw std::runtime_error("data must contain the variables in covariates");
    }
    const std::vector<double>& u = get_numeric_column(data,zj);
    for (i=0; i<n; i++) {
      zn[i][j+1] = u[i];
    }
  }
  Rcpp::Rcout << "[lifereg_purecpp] line 1843, after zn matrix construction" << std::endl;

  std::vector<double> weightn(n, 1.0);
  //higher level wrapper does not pass anything for weight so we skip this temporarily
  if (!data.weight.empty()) {
    weightn = data.weight;
    if (std::any_of(weightn.begin(), weightn.end(), [](double w) { return w <= 0; })) {
      throw std::runtime_error("weight must be greater than 0");
    }
  }

  std::vector<double> offsetn(n);
  if (!data.offset.empty()) { // if it does have offset data
    offsetn = data.offset;
    //this is a deep copy so its same as rcpp .clone(...)
  }

  std::vector<int> idn(n);
  if (data.id_raw.empty() || id == "") {
    idn = seq_cpp(1,n);
  } else {
    //NOTE: the data.id_raw is redundant check but we keep it for consistency 
    //will remove later
    //error checking for length n for id column
    if (data.id_raw.size() != static_cast<std::size_t>(n)) {
      throw std::runtime_error("ID column must have length n");
    }
    //error checking for all empty values in id column
    if (!data.id_raw.empty()) {
      bool all_blank = std::all_of(
        data.id_raw.begin(), data.id_raw.end(),
        [](auto &s){ return s.empty(); }
      );
      if (all_blank){
        throw std::runtime_error("ID column must contain non-empty values");
      }
    }
    group_index g = group_by(data.id_raw);
    idn = std::move(g.index);
  }

  //sort the data by rep
  std::vector<int> order = seq_cpp(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  // reorder the vectors and matrix by the order
  reorder(repn,order);
  reorder(stratumn,order);
  reorder(timen,order);
  reorder(time2n,order);
  reorder(eventn,order);
  reorder(weightn,order);
  reorder(offsetn,order);
  reorder(idn,order);
  zn = subset_matrix_by_row_cpp(zn, order);

  // exclude observations with missing valuesbeca
  std::vector<bool> sub(n,1);
  for (i=0; i<n; i++) {
    if (//(repn[i] == NA_INTEGER) || (stratumn[i] == NA_INTEGER) ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        (std::isnan(weightn[i])) || (std::isnan(offsetn[i])) )//||
        //(idn[i] == NA_INTEGER)) 
        {
      sub[i] = 0;
    }
    for (j=0; j<nvar-1; j++) {
      if (std::isnan(zn[i][j+1])) sub[i] = 0;
    }
  }

  //order = which_cpp
  order = which_true(sub);
  reorder(repn,order);
  reorder(stratumn,order);
  reorder(timen,order);
  reorder(time2n,order);
  reorder(eventn,order);
  reorder(weightn,order);
  reorder(offsetn,order);
  reorder(idn,order);
  zn = subset_matrix_by_row_cpp(zn, order);
  n = sum_bool_cpp(sub);

  // identify the locations of the unique values of rep
  std::vector<int> idx(1,0);
  for (i=1; i<n; i++) {
    if (repn[i] != repn[i-1]) {
      idx.push_back(i);
    }
  }

  int nreps = static_cast<int>(idx.size());
  idx.push_back(n);

  // variables in the output data sets
  std::vector<int> rep01 = seq_cpp(1,nreps);
  std::vector<int> nobs(nreps), nevents(nreps);
  std::vector<double> loglik0(nreps,NAN), loglik1(nreps,NAN);
  std::vector<int> niter(nreps);
  std::vector<bool> fails(nreps, false);

  std::vector<int> rep0(nreps*p,-1); //-1 not NAN for integer vectors
  std::vector<std::string> par0(nreps*p, "");
  std::vector<double> beta0(nreps*p,NAN), sebeta0(nreps*p,NAN), rsebeta0(nreps*p,NAN);
  std::vector<std::vector<double>> vbeta0(nreps*p, std::vector<double>(p,NAN)), rvbeta0(nreps*p, std::vector<double>(p,NAN));
  std::vector<double> lb0(nreps*p,NAN), ub0(nreps*p,NAN), prob0(nreps*p,NAN);
  std::vector<std::string> clparm0(nreps*p,"");
  //print check statement
  Rcpp::Rcout << "[lifereg_purecpp] line 1964, reached pragma parallelization" << std::endl;   // continues

  //if the crash dissapears when its a single thread that means we can assume its a OMP Rcpp interaction issue
  // for more notes about this specifically check below liferegpure_cpp function signature
  //#if defined(_OPENMP)
  //  omp_set_num_threads(1); // TEMP: make it single-threaded to confirm
  //#endif
  Rcpp::Rcout << "[lifereg_purecpp] line 1952 BIG CHECK n=" << n
            << " nreps=" << nreps
            << " nvar=" << nvar
            << " nstrata=" << nstrata
            << " p=" << p
            << " beta0.len=" << (nreps*p) << std::endl;
  //#pragma omp parallel for schedule(static) private(i,j,k)
  for (h=0; h<nreps; h++) {
    //these are declared bc each thread needs its own copy
    //and we cannot use the outer scope variables directly
    bool fail = false;
    liferegloopresult out;

    std::vector<int> q1 = seq_cpp(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    std::vector<int> stratum1 = subset_by_idx(stratumn, q1);
    std::vector<double> time1 = subset_by_idx(timen, q1);
    std::vector<double> time21 = subset_by_idx(time2n, q1);
    std::vector<int> event1 = subset_by_idx(eventn, q1);
    std::vector<double> weight1 = subset_by_idx(weightn, q1);
    std::vector<double> offset1 = subset_by_idx(offsetn, q1);
    std::vector<int> id1 = subset_by_idx(idn, q1);
    std::vector<std::vector<double>> z1 = subset_matrix_by_row_cpp(zn, q1);

    // unify right censored data with interval censored data
    std::vector<double> tstart(n1), tstop(n1);
    if (data.time2.empty()) {
      tstart = time1;
      for (i=0; i<n1; i++) {
        tstop[i] = event1[i] == 1 ? tstart[i] : NA_REAL;
      }
    } else {
      tstart = time1;
      tstop = time21;
    }

    std::vector<int> status(n1);
    for (i=0; i<n1; i++) {
      if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) &&
          (tstart[i] == tstop[i])) {
        status[i] = 1; // event
      } else if (!std::isnan(tstart[i]) && !std::isnan(tstop[i]) && (tstart[i] < tstop[i])) {
        status[i] = 3; // interval censoring
      } else if (std::isnan(tstart[i]) && !std::isnan(tstop[i])) {
        status[i] = 2; // left censoring
      } else if (!std::isnan(tstart[i]) && std::isnan(tstop[i])) {
        status[i] = 0; // right censoring
      } else {
        status[i] = -1; // exclude the observation
      }
    }

    nobs[h] = n1;
    nevents[h] = static_cast<int>(
    std::count_if(status.begin(), status.end(),
                [](int s){ return s == 1; })
    );

    // exclude records with invalid status
    std::vector<int> q2;
    for (size_t i = 0; i < status.size(); i++) {
      if (status[i] != -1) {
        q2.push_back(i);
      }
    }

    int n2 = static_cast<int>(q2.size());

    if (n2 < n1) {
      stratum1 = subset_by_idx(stratum1, q2);
      tstart = subset_by_idx(tstart, q2);
      tstop = subset_by_idx(tstop, q2);
      status = subset_by_idx(status, q2);
      weight1 = subset_by_idx(weight1, q2);
      offset1 = subset_by_idx(offset1, q2);
      id1 = subset_by_idx(id1, q2);
      z1 = subset_matrix_by_row_cpp(z1,q2);
    }

    // intercept only model
    std::vector<double> time0(n2);
    for (i=0; i<n2; i++) {
      if (status[i] == 0 || status[i] == 1) { // right censoring or event
        time0[i] = tstart[i];
      } else if (status[i] == 2) { // left censoring
        time0[i] = tstop[i];
      } else if (status[i] == 3) { // interval censoring
        time0[i] = (tstart[i] + tstop[i])/2;
      }
    }

    std::vector<double> y0 = time0;
    if ((dist1 == "exponential") || (dist1 == "weibull") ||
        (dist1 == "lognormal") || (dist1 == "loglogistic")) {
      for (auto& val : y0) {
        val = std::log(val);
      }
    }

    double int0 = mean_cpp(y0); // use the sample mean as the initial value of intercept
    double logsig0 = std::log(sd_cpp(y0)); // use the sample sd as the initial value of log(scale)

    std::vector<double> bint0(p);
    int ncolfit0 = dist1 == "exponential" ? 1 : nstrata + 1;
    std::vector<int> colfit0(ncolfit0);
    if (dist1 == "exponential") {
      bint0[0] = int0;
      ncolfit0 = 1;
      colfit0[0] = 0;
    } else {
      bint0[0] = int0;
      for (i=0; i<nstrata; i++) {
        bint0[nvar+i] = logsig0;
      }

      colfit0[0] = 0;
      for (i=0; i<nstrata; i++) {
        colfit0[i+1] = nvar+i;
      }
    }

    // last left off here
    // parameter estimates and standard errors for the null model
    aftparams_pure param = {dist1, stratum1, tstart, tstop, status, weight1,
                          offset1, z1, nstrata};

    liferegloopresult outint = liferegloop_cpp(p, bint0, &param, maxiter, eps,
                              colfit0, ncolfit0);
    std::vector<double> bint = outint.coef;
    std::vector<std::vector<double>> vbint = outint.var;

    std::vector<double> b(p);
    std::vector<std::vector<double>> vb(p,std::vector<double>(p));

    if (nvar > 1) {
      std::vector<int> colfit = seq_cpp(0,p-1);

      auto all_finite = [](const std::vector<double>& v){
        return std::all_of(v.begin(),v.end(), [](double x){ return std::isfinite(x); });
      };
      bool have_init = (init.size() == static_cast<size_t> (p) && all_finite(init));

      if(have_init){
        out = liferegloop_cpp(p,init, &param, maxiter, eps,
                              colfit, p);
      }
      else{
        out = liferegloop_cpp(p,bint, &param, maxiter, eps,
                              colfit, p);
      }

      fail = out.fail; 
      if(fail){ 
        std::vector<double> y1; // = y0 - offset;
        for(size_t i = 0; i < y0.size(); i++){
          y1.push_back(y0[i] - offset1[i]);
        }
        std::vector<std::vector<double>> v1(nvar, std::vector<double>(nvar));
        std::vector<double> u1(nvar);

        for(int i = 0; i < n2; i++){
          for(int j = 0; j < nvar; j++){
            for(int k = 0; k < nvar; k++){
              v1[j][k] += weight1[i] * (z1[i][j] * z1[i][k]);
            }
            u1[j] += weight1[i] * z1[i][j] * y1[i];
          }
        }

      
        // update the parameter vector
        double toler = 1e-12;
        i = cholesky2_cpp(v1, nvar, toler);
        chsolve2_cpp(v1, nvar, u1);
      
        // transform back to the unstandardized design matrix scale
        std::vector<double> binit(p);
        for(j = 0; j < nvar; j++){
          binit[j] = u1[j];
        }

        if(dist1 != "exponential"){
          double s = 0.0;
          for(i = 0; i < n2; i++){
            double pred = std::inner_product(z1[i].begin(), z1[i].end(), u1.begin(), 0.0);
            double r = y1[i] - pred;
            s += weight1[i] * r * r;
          }
          s = 0.5* std::log(s/sum_cpp(weight1)*n2/(n2-nvar));

          for(j = nvar; j < p; j++){
            binit[j] = s;
          }
        }
    
        out = liferegloop_cpp(p, binit, &param, maxiter, eps, colfit, p);
        fail = out.fail;
      }

      //moved this labelling part here to run it before the fail check continues;
      for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;

      if (i==0) {
        par0[h*p+i] = "(Intercept)";
      } else if (i < nvar) {
        par0[h*p+i] = covariates[i-1];
      } else {
        if (nstrata == 1) {
          par0[h*p+i] = "Log(scale)";
        } else {
          std::string str1 = "Log(scale ";
          std::string str2 = ")";
          par0[h*p+i] = str1 + std::to_string(i-nvar+1) + str2;
        }
      }
    }

      if(fail){ //throw error normally but we cant do that in a parallel loop
        niter[h] = out.iter;
        fails[h] = fail;
        loglik0[h] = outint.loglik;
        loglik1[h] = NAN;
        continue; 
      }
      b = out.coef;
      vb = out.var;
    } else{
      b = bint;
      vb = vbint;
      out  = outint;
      //niter[h] = out.iter;
    }

    niter[h] = out.iter;
    fails[h] = out.fail;
    //b = std::move(out.coef);
    //vb = std::move(out.var);
    //OPTIONAL: fit outputs for this h index with NaNs or Zeros
    //if(out.fail) continue;
    //b = out.coef;
    //vb = out.var;
    //std::vector<double> seb(p);
    //for (j=0; j<p; j++) {
      //seb[j] = std::sqrt(vb[j][j]);
    //}

    /*
    for (i=0; i<p; i++) {
      rep0[h*p+i] = h+1;

      if (i==0) {
        par0[h*p+i] = "(Intercept)";
      } else if (i < nvar) {
        par0[h*p+i] = covariates[i-1];
      } else {
        if (nstrata == 1) {
          par0[h*p+i] = "Log(scale)";
        } else {
          std::string str1 = "Log(scale ";
          std::string str2 = ")";
          par0[h*p+i] = str1 + std::to_string(i-nvar+1) + str2;
        }
      }
      //beta0[h*p+i] = b[i];
      //sebeta0[h*p+i] = seb[i];
      //for (j=0; j<p; j++) {
      //  vbeta0[h*p+i][j] = vb[i][j];
      //}
    }*/


    if(out.fail){
      loglik0[h] = outint.loglik;
      loglik1[h] = NAN;
    }
    else{
      std::vector<double> seb(p);
      for (j=0; j<p; j++) {
        seb[j] = std::sqrt(vb[j][j]);
      }

      for (i = 0; i< p; i++){
        beta0[h*p+i] = b[i];
        sebeta0[h*p+i] = seb[i];
        for (j=0; j<p; j++) {
          vbeta0[h*p+i][j] = vb[i][j];
        }
      }

      // robust variance estimates
      std::vector<double> rseb(p);  // robust standard error for betahat
      if (robust) {
        std::vector<std::vector<double>> ressco = f_ressco_1_cpp(p, b, &param);

        int nr; // number of rows in the score residual matrix
        if (data.id_raw.empty() || id == "") {
          for (i=0; i<n2; i++) {
            for (j=0; j<p; j++) {
              ressco[i][j] = weight1[i]*ressco[i][j];
            }
          }
          nr = n2;
        } else { // need to sum up score residuals by id
          //IntegerVector order = seq(0, n2-1);
          std::vector<int> order = seq_cpp(0, n2-1);
          std::sort(order.begin(), order.end(), [&](int i, int j) {
            return id1[i] < id1[j];
          });

          //IntegerVector id2 = id1[order];
          std::vector<int> id2 = subset_by_idx(id1, order);
          //IntegerVector idx(1,0);
          std::vector<int> idx(1,0);
          for (i=1; i<n2; i++) {
            if (id2[i] != id2[i-1]) {
              idx.push_back(i);
            }
          }

          int nids = static_cast<int>(idx.size());
          idx.push_back(n2);

          //NumericVector weight2 = weight1[order];
          std::vector<double> weight2 = subset_by_idx(weight1, order);

          std::vector<std::vector<double>> ressco2(nids,std::vector<double>(p));
          for (i=0; i<nids; i++) {
            for (j=0; j<p; j++) {
              for (k=idx[i]; k<idx[i+1]; k++) {
                ressco2[i][j] += weight2[k]*ressco[order[k]][j];
              }
            }
          }

          ressco = ressco2;  // update the score residuals
          nr = nids;
        }

        std::vector<std::vector<double>> D(nr,std::vector<double>(p)); // DFBETA
        for (i=0; i<nr; i++) {
          for (j=0; j<p; j++) {
            for (k=0; k<p; k++) {
              D[i][j] += ressco[i][k]*vb[k][j];
            }
          }
        }

        std::vector<std::vector<double>> rvb(p,std::vector<double>(p)); // robust variance matrix for betahat
        for (j=0; j<p; j++) {
          for (k=0; k<p; k++) {
            for (i=0; i<nr; i++) {
              rvb[j][k] += D[i][j]*D[i][k];
            }
          }
        }

        for (i=0; i<p; i++) {
          rseb[i] = std::sqrt(rvb[i][i]);
        }

        for (i=0; i<p; i++) {
          rsebeta0[h*p+i] = rseb[i];
          for (j=0; j<p; j++) {
            rvbeta0[h*p+i][j] = rvb[i][j];
          }
        }
      }

      // profile likelihood confidence interval for regression coefficients
      std::vector<double> lb(p), ub(p), prob(p);
      std::vector<std::string> clparm(p);

      //double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
      double zcrit = qnorm_cpp(1-alpha/2, true,false);
      if (plci) {
        double lmax = f_llik_1_cpp(p, b, &param);
        //double l0 = lmax - 0.5*R::qchisq(1-alpha, 1, 1, 0);
        double l0 = lmax - 0.5*qchisq_cpp(1-alpha, 1, true,false);

        for (k=0; k<p; k++) {
          lb[k] = liferegplloop_cpp(p, b, &param, maxiter, eps, k, -1, l0);
          ub[k] = liferegplloop_cpp(p, b, &param, maxiter, eps, k, 1, l0);

          std::vector<int> colfit1(p-1);
          for (i=0; i<k; i++) {
            colfit1[i] = i;
          }
          for (i=k+1; i<p; i++) {
            colfit1[i-1] = i;
          }

          std::vector<double> b0(p);
          liferegloopresult out0 = liferegloop_cpp(p, b0, &param, maxiter, eps, colfit1, p-1);
          double lmax0 = out0.loglik; 
          //prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
          prob[k] = pchisq_cpp(-2*(lmax0 - lmax), 1.0, false, false);
          clparm[k] = "PL";
        }
      } else {
        for (k=0; k<p; k++) {
          if (!robust) {
            lb[k] = b[k] - zcrit*seb[k];
            ub[k] = b[k] + zcrit*seb[k];
            //prob[k] = R::pchisq(pow(b[k]/seb[k], 2), 1, 0, 0);
            prob[k] = pchisq_cpp(std::pow(b[k]/seb[k], 2), 1.0, false, false);
          } else {
            lb[k] = b[k] - zcrit*rseb[k];
            ub[k] = b[k] + zcrit*rseb[k];
            //prob[k] = R::pchisq(pow(b[k]/rseb[k], 2), 1, 0, 0);
            prob[k] = pchisq_cpp(std::pow(b[k]/rseb[k], 2), 1.0, false, false);
          }
          clparm[k] = "Wald";
        }
      }

      for (i=0; i<p; i++) {
        lb0[h*p+i] = lb[i];
        ub0[h*p+i] = ub[i];
        prob0[h*p+i] = prob[i];
        clparm0[h*p+i] = clparm[i];
      }

      // log-likelihoods
      //loglik0[h] = outint["loglik"];
      loglik0[h] = outint.loglik;
      //loglik1[h] = out["loglik"];
      loglik1[h] = out.loglik;
    }
  }
  Rcpp::Rcout << "[lifereg_purecpp] line 2425, after parallel for loop" << std::endl;

  //end of parllel for loop so AFTER this line rcpp  
  // convert the vectors to Rcpp types
  //NumericVector expbeta0 = exp(beta0);
  Rcpp::Rcout << "[liferegloop_cpp] beta0.size()=" << beta0.size() << std::endl;
  NumericVector expbeta0(beta0.size());
  for (size_t i=0; i<beta0.size(); i++) {
    if (i < 5) Rcpp::Rcout << "beta0[" << i << "]=" << beta0[i] << std::endl;
    expbeta0[i] = std::exp(beta0[i]);
  }

  Rcpp::Rcout << "[liferegloop_cpp] line 2434, expbeta0 done" << std::endl;

  Rcpp::Rcout << "[liferegloop_cpp] line 2436, converting beta0, sebeta0, rsebeta0, z0 to NumericVector" << std::endl;
  NumericVector beta0_r = NumericVector(beta0.begin(), beta0.end());
  NumericVector sebeta0_r = NumericVector(sebeta0.begin(), sebeta0.end());
  NumericVector rsebeta0_r = NumericVector(rsebeta0.begin(), rsebeta0.end());
  NumericVector z0(nreps*p);
  Rcpp::Rcout << "[liferegloop_cpp] line 2441, beta0_r, sebeta0_r, rsebeta0_r done" << std::endl;
  if (!robust) z0 = beta0_r/sebeta0_r;
  else z0 = beta0_r/rsebeta0_r;

  //these two bottom dataframes are the key outputs
  // anything inside these dataframes must be in rcpp types so we convert our stl containers into rcpp types again
  // we can completely eradicate this part whenever we bootstrap phregcpp to make it return stl types instead of dataframes
  
  //for loop checking if data wrangling error
  for (size_t r = 0; r < vbeta0.size(); ++r) {
    if (vbeta0[r].size() != (size_t)p) {
      Rcpp::stop("vbeta0 row %zu has length %zu, expected %d", r, vbeta0[r].size(), p);
    }
  }
  Rcpp::Rcout << "[liferegloop_cpp] line 2451, wrapping start" << std::endl;
  //for sumstat
  IntegerVector nobs_r = Rcpp::wrap(nobs);
  IntegerVector nevents_r = Rcpp::wrap(nevents);
  NumericVector loglik0_r = Rcpp::wrap(loglik0);
  NumericVector loglik1_r = Rcpp::wrap(loglik1);
  IntegerVector niter_r = Rcpp::wrap(niter);

  Rcpp::Rcout << "[liferegloop_cpp] line 2459, sumstat wrap done parest wrap start" << std::endl;
  //for parest
  NumericMatrix vbeta0_r = Rcpp::wrap(vbeta0);
  Rcpp::Rcout << "[liferegloop_cpp] line 2462, vbeta0_r wrap done" << std::endl;
  StringVector par0_r = Rcpp::wrap(par0);
  NumericVector lb0_r = Rcpp::wrap(lb0);
  NumericVector ub0_r = Rcpp::wrap(ub0);
  NumericVector prob0_r = Rcpp::wrap(prob0);
  StringVector clparm0_r = Rcpp::wrap(clparm0);
  NumericMatrix rvbeta0_r = Rcpp::wrap(rvbeta0);
  LogicalVector fails_r = Rcpp::wrap(fails);

  Rcpp::Rcout << "[liferegloop_cpp] line 2468, wrapping end" << std::endl;
  //end of conversion
  DataFrame sumstat = List::create(
    _["n"] = nobs_r,
    _["nevents"] = nevents_r,
    _["loglik0"] = loglik0_r,
    _["loglik1"] = loglik1_r,
    _["niter"] = niter_r,
    _["dist"] = dist1,
    _["p"] = p,
    _["nvar"] = nvar-1,
    _["robust"] = robust,
    _["fail"] = fails_r);
  DataFrame parest;
  if (!robust) {
    parest = DataFrame::create(
      _["param"] = par0_r,
      _["beta"] = beta0_r,
      _["sebeta"] = sebeta0_r,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = vbeta0_r,
      _["lower"] = lb0_r,
      _["upper"] = ub0_r,
      _["p"] = prob0_r,
      _["method"] = clparm0_r);
  } else {
    parest = DataFrame::create(
      _["param"] = par0_r,
      _["beta"] = beta0_r,
      _["sebeta"] = rsebeta0_r,
      _["z"] = z0,
      _["expbeta"] = expbeta0,
      _["vbeta"] = rvbeta0_r,
      _["lower"] = lb0_r,
      _["upper"] = ub0_r,
      _["p"] = prob0_r,
      _["method"] = clparm0_r,
      _["sebeta_naive"] = sebeta0_r,
      _["vbeta_naive"] = vbeta0_r);
  }

  //we can ignore this bit because has_rep will ALWAYS be false from higher level wrapper logic
  /*
  if (has_rep) {
    for (i=0; i<p_rep; i++) {
      String s = rep[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector repwi = u_rep[s];
        sumstat.push_back(repwi[rep01-1], s);
        parest.push_back(repwi[rep0-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector repwn = u_rep[s];
        sumstat.push_back(repwn[rep01-1], s);
        parest.push_back(repwn[rep0-1], s);
      } else if (TYPEOF(data[rep]) == STRSXP) {
        StringVector repwc = u_rep[s];
        sumstat.push_back(repwc[rep01-1], s);
        parest.push_back(repwc[rep0-1], s);
      }
    }
  }
  */
  //add functionality for this in another function that requires it but here it is not needed.

  List result = List::create(
    _["sumstat"] = sumstat,
    _["parest"] = parest);

  return result;
}