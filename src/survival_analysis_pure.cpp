#include <Rcpp.h>
#include <R_ext/Applic.h>
//#include "utilities.h"
//#include "survival_analysis.h"
using namespace Rcpp;

//custom includes (basically just utilities and survival analysis in pure cpp)
#include "survival_analysis_pure.h"
#include "utilities_pure.h"
//using namespace Rcpp;


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
/*
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
*/

int find_aft_index(const trial_data& d, const std::string& col_name) {
  for (size_t i = 0; i < d.aft_names.size(); ++i) {
    if (d.aft_names[i] == col_name) {
      return static_cast<int>(i);
    }
  }
  return -1; // not found
}

NumericMatrix to_matrix(const std::vector<std::vector<double>>& matrix) {
    int n_rows = matrix.size();
    int n_cols = matrix.front().size();

    //check for rectangle shape first to avoid issues
    for(int row = 0; row < n_rows; row++){
      int cols_here = matrix[row].size();
      if(cols_here != n_cols){
        stop("error in to_matrix");
      }
    }

    NumericMatrix out(n_rows, n_cols);

    for(int r = 0 ; r < n_rows; r++){
      for(int c = 0; c < n_cols; c++){
        out(r,c) = matrix[r][c];
      }
    }
  return out;
} 

/*
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
*/

 

/*
//
//
//
*/

// score vector
std::vector<double> f_score_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param-> nrows;
  int nvar = param-> ncols;
  int person, i, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
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
  std::vector<double> z(nvar);
  for (person = 0; person < n; person++) {
    double wt = param->weight[person];
    double sigma = sig[person];
    //std::vector<double> z = (param->z[person] / sigma); //turned into the bottom 2 lines
    //std::vector<double> z = param->z[person];      // copy the row
    for (int jj = 0; jj < nvar; ++jj) {
      z[jj] = param->z(person, jj);   // raw
    }
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
        double c0 = 1 - 2*plogis_cpp(u,true,false);
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
        double c0 = 1 - 2*plogis_cpp(u,true,false);
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
        double q1 = pnorm_cpp(v, false), q2 = pnorm_cpp(u, false);
        double c1 = wt*(d1 - d2)/(q1 - q2);
        for (i=0; i<nvar; i++) {
          score[i] += c1*z[i];
        }
        score[k] += wt*(d1*v - d2*u)/(q1 - q2);
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        double d1 = dlogis_cpp(v,false), d2 = dlogis_cpp(u,false);
        double q1 = plogis_cpp(v,false), q2 = plogis_cpp(u,false);
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
        double c1 = wt*(-plogis_cpp(u, false));
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
        double c1 = wt*(-plogis_cpp(u, false));
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
        double c1 = wt*(dnorm_cpp(v, false)/pnorm_cpp(v, false));
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
        double c1 = wt*(dnorm_cpp(v, false)/pnorm_cpp(v, false));
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
  int n = param-> nrows; 
  int nvar = param-> ncols; 
  int person, i, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
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
        loglik += wt*std::log(pnorm_cpp(v,false) - pnorm_cpp(u,false));
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,false) - plogis_cpp(u,false));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,false) - pnorm_cpp(u,false));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,false) - plogis_cpp(u,false));
      }
    } else if (param->status[person] == 2) { // upper used as left censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*std::log(1.0 - std::exp(-std::exp(u)));
      } else if (param->dist == "lognormal") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(u,true));
      } else if (param->dist == "loglogistic") {
        double u = (std::log(param->tstop[person]) - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(u,true));
      } else if (param->dist == "normal") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(u,true));
      } else if (param->dist == "logistic") {
        double u = (param->tstop[person] - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(u,true));
      }
    } else if (param->status[person] == 0) { // lower used as right censoring
      if (param->dist == "exponential" || param->dist == "weibull") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*(-std::exp(v));
      } else if (param->dist == "lognormal") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,false));
      } else if (param->dist == "loglogistic") {
        double v = (std::log(param->tstart[person]) - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,false));
      } else if (param->dist == "normal") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(pnorm_cpp(v,false));
      } else if (param->dist == "logistic") {
        double v = (param->tstart[person] - eta[person])/sigma;
        loglik += wt*std::log(plogis_cpp(v,false));
      }
    }
  }

  return loglik;
}

// observed information matrix
std::vector<std::vector<double>> f_info_1_cpp(int p, std::vector<double> par, void *ex) {
  aftparams_pure *param = (aftparams_pure *) ex;
  int n = param -> nrows;
  int nvar = param -> ncols;
  int person, i, j, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
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
  std::vector<double> z(nvar);
  for (person = 0; person < n; person++) {
    for (int jj = 0; jj < nvar; ++jj) {
      z[jj] = param->z(person, jj);
    }
    double wt = param->weight[person];
    double sigma = sig[person];
    //std::vector<double> z = param->z[person];
    
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
                        1 - 2*plogis_cpp(u, true, false));
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
                        1 - 2*plogis_cpp(u, true, false));
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
        double q1 = pnorm_cpp(v,false), q2 = pnorm_cpp(u, false);
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
  int n = param-> nrows;
  int nvar = param-> ncols;
  int person, i, k;

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = param->offset[person];
    for (i=0; i<nvar; i++) {
      eta[person] += par[i]*param->z(person,i);
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

  std::vector<double> z(nvar); // re use buffer :D
  for (person = 0; person < n; person++) {
    
    /*double sigma = sig[person];
    std::vector<double> z = param->z[person];
    for(double &v : z) {
      v /= sigma;
    }
    k = param->strata[person] + nvar - 1; */
    double sigma = sig[person];
    if (!std::isfinite(sigma) || sigma == 0.0) sigma = 1e-12;

    // fill reusable scaled row buffer once
    for (int ii = 0; ii < nvar; ++ii) {
      z[ii] = param->z(person, ii) / sigma;
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
        double c1 = (-dnorm_cpp(u,false)/pnorm_cpp(u,true));
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
        double c1 = dnorm_cpp(v,false)/pnorm_cpp(v,false);
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
  int n = param-> nrows;
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
  int nsub = static_cast<int>(param -> nrows); // reminder that z is a vector so grabbing the .size() gives size_t not int, not a massive difference 
  int nvar = static_cast<int>(param -> ncols); // important
  //std::vector<std::vector<double>> z1 = param -> z;
  //anywhere z1 is called use params -> z(index1, index2)

  // standardize the design matrix
  std::vector<double> mu(nvar, 0.0), sigma(nvar, 1.0);
  //std::vector<std::vector<double>> z2(nsub, std::vector<double>(nvar, 0.0));
  MatrixRM z2 (nsub,nvar); //{ nsub, nvar, std::vector<double>((size_t)nsub * nvar, 0.0) };

  for (i=0; i<nvar; i++) {
    std::vector<double> u(nsub);
    for (j=0; j < nsub; j++) {
      u[j] = param -> z(j,i);
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
      z2(i,j) = (param -> z(i,j) - mu[j]) / sigma[j];
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
  
  //aftparams_pure para = {param->dist, param->strata, param->tstart, param->tstop,
  //                  param->status, param->weight, param->offset, nstrata, z2};
  // Identity map
  std::vector<int> rows2(nsub);
  for (int i = 0; i < nsub; ++i) rows2[i] = i;

  aftparams_pure para;
  para.dist   = param->dist;
  para.strata = param->strata;
  para.tstart = param->tstart;
  para.tstop  = param->tstop;
  para.status = param->status;
  para.weight = param->weight;
  para.offset = param->offset;
  para.nstrata = nstrata;

  para.zbase = &z2;
  para.rows  = &rows2;
  para.nrows = nsub;
  para.ncols = nvar;

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

/*Notes and things to consider:
  - The function is a direct translation of the R function liferegcpp
  - The function expects a trial_data object which contains the necessary data
  - The function handles various distributions and covariates
  - Error handling is returns a "fail" liferegOut object
  - The function returns a struct now, to be purely C++
  - The function does NOT need to handle hasVariable checks for time2, rep, and stratum
    because the higher level wrapper will handle that logic
*/

// [[Rcpp::export]]
liferegOut lifereg_purecpp(
    //const DataFrame data expects trial_data with 
    const trial_data& data,
    const std::vector<std::string>& covariates = {},
    const std::string dist = "weibull",
    const std::vector<double>& init = {},
    const bool robust = false,
    const bool plci = false,
    const double alpha = 0.05,
    const int maxiter = 50,
    const double eps = 1.0e-9) {
  int h, i, j, k, n = data.pps.size(); //.nrows();

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

  // initialized like this because repn and stratumn will never pass arguments!
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
  std::vector<double> timen = timenz; 
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

  MatrixRM zn (n,nvar); //{n, nvar, std::vector<double>((size_t)n*nvar,0.0)};
  const MatrixRM& A = *data.aft;
  const std::vector<int>& rp = *data.order_pp;

  for (int i=0; i<n; i++) {
    zn(i,0) = 1; // intercept
  }
  if (nvar > 1) {
    for (int i = 0; i < n; ++i) {
      zn(i, 1) = static_cast<double>(data.swtrt[i]);
    }
  }
  if(static_cast<int>(rp.size())!= n){
    throw std::runtime_error("order_pp size mistmatched");
    //error handle
  }

  for (int col = 0; col < A.ncols(); ++col) {
    for (int i = 0; i < n; ++i) {
        // add 2 not 1 because of swtrt
        zn(i, 2 + col) = A(rp[i], col);
      }
  }

  std::vector<int> row_idx(n);
  std::iota(row_idx.begin(),row_idx.end(),0);

  std::vector<double> weightn(n, 1.0);
  //higher level wrapper does not pass anything for weight so we skip this temporarily
  if (!data.weight.empty()) {
    weightn = data.weight;
    if (std::any_of(weightn.begin(), weightn.end(), [](double w) { return w <= 0; })) {
      throw std::runtime_error("weight must be greater than 0");
    }
  }

  std::vector<double> offsetn(n);
  if (!data.offset.empty()) { 
    offsetn = data.offset;
    //this is a deep copy so its same as rcpp .clone(...)
  }

  std::vector<int> idn(n);
  if (data.id_raw.empty()) {
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
  //zn = subset_matrix_by_row_cpp(zn, order);
  //replaced with
  reorder(row_idx,order);

  // exclude observations with missing valuesbeca
  std::vector<bool> sub(n,1);
  for (i=0; i<n; i++) {
    if (//(repn[i] == NA_INTEGER) || (stratumn[i] == NA_INTEGER) ||
        (std::isnan(timen[i]) && std::isnan(time2n[i])) ||
        (std::isnan(weightn[i])) || (std::isnan(offsetn[i])) )//||
        //(idn[i] == NA_INTEGER)) 
        {
      sub[i] = false;
      continue;
    }
    // covariate NA checks: use original row in zn via row_idx
    int r = row_idx[i];
    for (j=0; j<nvar-1; j++) {
      if (std::isnan(zn(r,j+1))){ 
        sub[i] = false;
        break;
      }
    }
  }

  std::vector<int> keep = which_true(sub);

  repn     = subset_by_idx(repn, keep);
  stratumn = subset_by_idx(stratumn, keep);
  timen    = subset_by_idx(timen, keep);
  time2n   = subset_by_idx(time2n, keep);
  eventn   = subset_by_idx(eventn, keep);
  weightn  = subset_by_idx(weightn, keep);
  offsetn  = subset_by_idx(offsetn, keep);
  idn      = subset_by_idx(idn, keep);

  // NEW: subset row mapping
  row_idx  = subset_by_idx(row_idx, keep);

  n = (int)keep.size();
  if (n == 0) {
    // fail early (no data left)
    
    Rcpp::stop("All observations filtered out due to missing values.");
  }

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
  
  // Rcpp::Rcout << "[lifereg_purecpp] line 1905, before h loop" << std::endl;
  for (h=0; h<nreps; h++) {
    bool fail = false;
    liferegloopresult out;
    aftparams_pure param;

    std::vector<int> q1 = seq_cpp(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    std::vector<int> stratum1 = subset_by_idx(stratumn, q1);
    std::vector<double> time1 = subset_by_idx(timen, q1);
    std::vector<double> time21 = subset_by_idx(time2n, q1);
    std::vector<int> event1 = subset_by_idx(eventn, q1);
    std::vector<double> weight1 = subset_by_idx(weightn, q1);
    std::vector<double> offset1 = subset_by_idx(offsetn, q1);
    std::vector<int> id1 = subset_by_idx(idn, q1);
    //std::vector<std::vector<double>> z1 = subset_matrix_by_row_cpp(zn, q1);
    int start = idx[h];
    int len   = idx[h+1] - idx[h];
    std::vector<int> rows1;
    rows1.reserve(len);
    for (int t = start; t < start + len; ++t) {
      rows1.push_back(row_idx[t]);  // row_idx already accounts for earlier reorder/subset
    }
    //RowRangeViewT<RowIndexView> z1{ &zn_view, start, len };

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
      //z1 = subset_matrix_by_row_cpp(z1,q2);
      rows1 = subset_by_idx(rows1, q2);
      param.nrows = rows1.size();
      //std::vector<int> rep_rows;              // original zn rows for the rep block
      //rep_rows.reserve(len);
      //for (int t = start; t < start + len; ++t) rep_rows.push_back(row_idx[t]);

      // now filter within the rep block:
      //std::vector<int> rep_rows2;
      //rep_rows2.reserve(q2.size());
      //for (int k : q2) rep_rows2.push_back(rep_rows[k]);

      //RowIndexView z1{ &zn, &rep_rows2 };
    }
    //RowIndexView z1 {&A, &rows1};
    RowIndexView z1 {&zn, &rows1}; // z_design(i,j) -> zn(rows1[i], j)

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

    // parameter estimates and standard errors for the null model
    //aftparams_pure param = {dist1, stratum1, tstart, tstop, status, weight1,
    //                      offset1, nstrata, z1};
    //int start = idx[h];
    //int len   = idx[h+1] - idx[h];

    // Build mapping from local rep row i -> original zn row

    
    param.dist   = dist1;
    param.strata = std::move(stratum1);
    param.tstart = std::move(tstart);
    param.tstop  = std::move(tstop);
    param.status = std::move(status);
    param.weight = std::move(weight1);
    param.offset = std::move(offset1);
    param.nstrata = nstrata;

    param.zbase = &zn;
    param.rows  = &rows1;
    param.nrows = static_cast<int>(rows1.size());      // or n2 after q2 filtering
    param.ncols = nvar;     // p? careful: nvar for covariates, p for full parameter length; z cols should be nvar

    liferegloopresult outint = liferegloop_cpp(p, bint0, &param, maxiter, eps,
                              colfit0, ncolfit0);
    std::vector<double> bint = outint.coef;
    std::vector<std::vector<double>> vbint = outint.var;

    std::vector<double> b(p);
    std::vector<std::vector<double>> vb(p,std::vector<double>(p));

    if (nvar > 1) {
      std::vector<int> colfit = seq_cpp(0,p-1);

      auto all_finite = [](const std::vector<double>& v){
        return std::all_of(v.begin(),v.end(), [](double x){ return !std::isnan(x); });
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
        //bc we do std::move we can just instead reference param.reference incase we create "use after move issues"
        std::vector<double> y1; // = y0 - offset;
        for(size_t i = 0; i < y0.size(); i++){
          y1.push_back(y0[i] - param.offset[i]);
        }
        std::vector<std::vector<double>> v1(nvar, std::vector<double>(nvar));
        std::vector<double> u1(nvar);

        for(int i = 0; i < n2; i++){
          for(int j = 0; j < nvar; j++){
            for(int k = 0; k < nvar; k++){
              v1[j][k] += param.weight[i] * (z1(i,j) * z1(i,k));
            }
            u1[j] += param.weight[i] * z1(i,j) * y1[i];
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
            //double pred = std::inner_product(z1[i].begin(), z1[i].end(), u1.begin(), 0.0);
            //change bc of the view replacement
            double pred = 0.0;
            for (int j = 0; j < nvar; ++j) {
              pred += z1(i, j) * u1[j];
            }
            double r = y1[i] - pred;
            s += param.weight[i] * r * r;
          }
          s = 0.5* std::log(s/sumdouble_cpp(param.weight)*n2/(n2-nvar));

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
        if (data.id_raw.empty()) {
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
      loglik0[h] = outint.loglik;
      loglik1[h] = out.loglik;
    }
  }
  
  std::vector<double> expbeta0(beta0.size());
  for (size_t i=0; i<beta0.size(); i++) {
    expbeta0[i] = std::exp(beta0[i]);
  }

  std::vector<double> z0(nreps*p);
  for(int i = 0; i < nreps*p; i++){
    z0[i] = beta0[i]/sebeta0[i];
  }

  liferegOut result; 
  result.fail = false; // any error handling beforehand should have creates result and passed back fail;
  result.beta = beta0;
  result.sebeta = sebeta0;
  // add more stuff to return if needed/wanted

  return result;
}

// bottom are phregcpp functions and helpers

// score vector
std::vector<double> f_score_2_cpp(int p, const std::vector<double>& par, void *ex) {
  coxparams_purecpp *param = (coxparams_purecpp *) ex;
  // local aliases (readability + less pointer noise)
  const auto& n = param -> n;
  const auto& offset = *param->offset;
  const auto& weight = *param->weight;
  const auto& tstart = *param->tstart;
  const auto& tstop  = *param->tstop;
  const auto& event  = *param->event;
  const auto& stratum= *param->stratum;
  const auto& order  = *param->order;
  // we are leaving the call param ->z(i,j) because its a listview but this might change in the future ? 

  int i, k, person;

  std::vector<double> u(p,0.0);       // score vector
  double dtime;                   // distinct time
  int ndead = 0;                  // number of deaths at this time point
  double deadwt = 0;              // sum of weights for the deaths
  double meanwt;                  // average weight for the deaths
  double risk;                    // weighted risk, w*exp(zbeta)
  double denom = 0;               // s0(beta,k,t)
  double denom2 = 0;              // sum of weighted risks for the deaths
  std::vector<double> a(p,0.0);       // s1(beta,k,t)
  std::vector<double> a2(p,0.0);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;                    // zbar(beta,k,t)

  std::vector<double> eta(n,0.0);
  for (person = 0; person < n; person++) {
    eta[person] = offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int istrata = stratum[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < n; ) {
    if (stratum[person] != istrata) { // hit a new stratum
      istrata = stratum[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      std::fill(a.begin(), a.end(), 0.0);
      std::fill(a2.begin(), a2.end(), 0.0); // IMPORTANT
    }

    dtime = tstop[person];
    while (person < n && tstop[person] == dtime) {
      // walk through this set of tied times
      risk = weight[person]*std::exp(eta[person]);
      if (event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) {
          a[i] += risk*param->z(person,i);
        }
      } else {
        ndead++;
        deadwt += weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          u[i] += weight[person]*param->z(person,i);
        }
      }

      person++;

      if (person < n && stratum[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<n; i1++) {
      p1 = order[i1];
      if (tstart[p1] < dtime || stratum[p1] != istrata) break;
      risk = weight[p1]*std::exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
      }
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        if (!(denom > 0.0) || !std::isfinite(denom)) return std::vector<double>(p, std::numeric_limits<double>::quiet_NaN());
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          u[i] -= deadwt*xbar;
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          if (!(denom > 0.0) || !std::isfinite(denom)) return std::vector<double>(p, std::numeric_limits<double>::quiet_NaN());
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            u[i] -= meanwt*xbar;
          }
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
      }
    }
  }

  return u;
}


// observed information matrix
MatrixRM f_info_2_cpp(int p, const std::vector<double>& par, void *ex) {
  coxparams_purecpp *param = (coxparams_purecpp *) ex;
  // local aliases (readability + less pointer noise)
  const auto& n = param -> n;
  const auto& offset = *param->offset;
  const auto& weight = *param->weight;
  const auto& tstart = *param->tstart;
  const auto& tstop  = *param->tstop;
  const auto& event  = *param->event;
  const auto& stratum= *param->stratum;
  const auto& order  = *param->order;

  int i, j, k, person;

  MatrixRM imat (p,p);//{p,p,std::vector<double>(p*p)};  // information matrix
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths
  std::vector<double> a(p);       // s1(beta,k,t)
  std::vector<double> a2(p);      // sum of w*exp(zbeta)*z for the deaths
  double xbar;              // zbar(beta,k,t)
  MatrixRM cmat (p,p); //{p,p,std::vector<double>(p*p)};  // s2(beta,k,t)
  MatrixRM cmat2 (p,p); //{p,p,std::vector<double>(p*p)}; // sum of w*exp(zbeta)*z*z' for the deaths

  std::vector<double> eta(n);
  for (person = 0; person < n; person++) {
    eta[person] = offset[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int istrata = stratum[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < n; ) {
    if (stratum[person] != istrata) {  // hit a new stratum
      istrata = stratum[person]; // reset temporary variables
      i1 = person;
      denom = 0;
      for (i=0; i<p; i++) {
        a[i] = 0;
        for (j=0; j<=i; j++) {
          cmat(i,j) = 0;
        }
      }
    }

    dtime = tstop[person];
    while (person < n && tstop[person] == dtime) {
      // walk through this set of tied times
      risk = weight[person]*std::exp(eta[person]);
      if (event[person] == 0) {
        denom += risk;
        for (i=0; i<p; i++) {
          a[i] += risk*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat(i,j) += risk*param->z(person,i)*param->z(person,j);
          }
        }
      } else {
        ndead++;
        deadwt += weight[person];
        denom2 += risk;
        for (i=0; i<p; i++) {
          a2[i] += risk*param->z(person,i);
          for (j=0; j<=i; j++) {
            cmat2(i,j) += risk*param->z(person,i)*param->z(person,j);
          }
        }
      }

      person++;

      if (person < n && stratum[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1<n; i1++) {
      p1 = order[i1];
      if (tstart[p1] < dtime || stratum[p1] != istrata) break;
      risk = weight[p1]*std::exp(eta[p1]);
      denom -= risk;
      for (i=0; i<p; i++) {
        a[i] -= risk*param->z(p1,i);
        for (j=0; j<=i; j++) {
          cmat(i,j) -= risk*param->z(p1,i)*param->z(p1,j);
        }
      }
    }

    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        for (i=0; i<p; i++) {
          a[i] += a2[i];
          xbar = a[i]/denom;
          for (j=0; j<=i; j++) {
            cmat(i,j) += cmat2(i,j);
            imat(i,j) += deadwt*(cmat(i,j) - xbar*a[j])/denom;
          }
        }
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          for (i=0; i<p; i++) {
            a[i] += a2[i]/ndead;
            xbar = a[i]/denom;
            for (j=0; j<=i; j++) {
              cmat(i,j) += cmat2(i,j)/ndead;
              imat(i,j) += meanwt*(cmat(i,j) - xbar*a[j])/denom;
            }
          }
        }
      }

      // reset for next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
      for (i=0; i<p; i++) {
        a2[i] = 0;
        for (j=0; j<=i; j++) {
          cmat2(i,j) = 0;
        }
      }
    }
  }

  for (i=0; i<p-1; i++) {
    for (j=i+1; j<p; j++) {
      imat(i,j) = imat(j,i);
    }
  }

  return imat;
}

// define functions in likelihood inference for the Cox model
// algorithms adapted from the coxph function in the survival package

// we can remove allocation of eta later to reuse a scratch buffer to avoid copies/repeated allocations
// log likelihood
double f_llik_2_cpp(int p, const std::vector<double>& par, const void *ex) {
  coxparams_purecpp *param = (coxparams_purecpp *) ex;
  int i, k, person;

  double loglik = 0;        // log-likelihood value
  double dtime;             // distinct time
  int ndead = 0;            // number of deaths at this time point
  double deadwt = 0;        // sum of weights for the deaths
  double meanwt;            // average weight for the deaths
  double risk;              // weighted risk, w*exp(zbeta)
  double denom = 0;         // s0(beta,k,t)
  double denom2 = 0;        // sum of weighted risks for the deaths

  std::vector<double> eta(param->n);
  for (person = 0; person < param->n; person++) {
    eta[person] = (*param -> offset)[person];
    for (i=0; i<p; i++) {
      eta[person] += par[i]*param->z(person,i);
    }
  }

  int istrata = (*param->stratum)[0]; // current stratum
  int i1 = 0;                     // index of descending start time
  int p1;                         // corresponding location in input data
  for (person = 0; person < param->n; ) {
    if ( (*param->stratum)[person] != istrata) { // hit a new stratum
      istrata = (*param->stratum)[person]; // reset temporary variables
      i1 = person;
      denom = 0;
    }

    dtime = (*param->tstop)[person];
    while (person < param->n && (*param->tstop)[person] == dtime) {
      // walk through this set of tied times
      risk = (*param->weight)[person]*std::exp(eta[person]);
      if ( (*param->event)[person] == 0) {
        denom += risk;
      } else {
        ndead++;
        deadwt += (*param->weight)[person];
        denom2 += risk;
        loglik += (*param->weight)[person]*eta[person];
      }

      person++;

      if (person < param->n && (*param->stratum)[person] != istrata) break;
    }

    // remove subjects no longer at risk
    for (; i1 < param-> n; i1++) {
      p1 = (*param->order)[i1];
      if ( (*param->tstart)[p1] < dtime || (*param->stratum)[p1] != istrata) break;
      risk = (*param->weight)[p1]*std::exp(eta[p1]);
      denom -= risk;
    }

    
    // add to the main terms
    if (ndead > 0) {
      if (param->method == 0 || ndead == 1) {
        denom += denom2;
        if (denom <= 0.0 || !std::isfinite(denom)) return std::numeric_limits<double>::quiet_NaN();
        loglik -= deadwt*std::log(denom);
      } else {
        meanwt = deadwt/ndead;
        for (k=0; k<ndead; k++) {
          denom += denom2/ndead;
          if (denom <= 0.0 || !std::isfinite(denom)) return std::numeric_limits<double>::quiet_NaN();
          loglik -= meanwt*std::log(denom);
        }
      }

      // reset for the next death time
      ndead = 0;
      deadwt = 0;
      denom2 = 0;
    }
  }

  return loglik;
}

// underlying optimization algorithm for phreg
phregLoopOut phregloop_cpp(int p, const std::vector<double>& par, void *ex,
               int maxiter, double eps, bool firth,
               const std::vector<int>& colfit, int ncolfit) {
  coxparams_purecpp *param = (coxparams_purecpp *) ex;
  int i, j, iter, halving = 0;
  bool fail;

  const RowIndexViewView& z = param->z;

  double toler = 1e-12;
  std::vector<double> beta(p), newbeta(p);
  double loglik, newlk = 0;
  std::vector<double> u(p);
  MatrixRM imat (p,p); //{p,p, std::vector<double>(p*p)};
  std::vector<double> u1(ncolfit);
  MatrixRM imat1 (ncolfit,ncolfit); //{ncolfit, ncolfit, std::vector<double>(ncolfit*ncolfit)};

  // initial beta and log likelihood
  for (i=0; i<p; i++) {
    beta[i] = par[i];
  }

  if (!firth) {
    loglik = f_llik_2_cpp(p, beta, param);
    u = f_score_2_cpp(p, beta, param);
  } else{
    #if TRT_PHREG_ENABLE_FIRTH
    loglik = f_pen_llik_2(p, beta, param);
    u = f_pen_score_2(p, beta, param);
    #endif
  }
  for (i=0; i<ncolfit; i++) {
    u1[i] = u[colfit[i]];
  }

  imat = f_info_2_cpp(p, beta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  i = cholesky2_cpp(imat1, ncolfit, toler);

  chsolve2_cpp(imat1, ncolfit, u1);

  std::fill(u.begin(),u.end(),0.0);
  for (i=0; i<ncolfit; i++) {
    u[colfit[i]] = u1[i];
  }

  for (i=0; i<p; i++) {
    newbeta[i] = beta[i] + u[i];
  }

  for (iter=0; iter<maxiter; iter++) {
    // new log likelihood
    if (!firth) newlk = f_llik_2_cpp(p, newbeta, param);
    else{ 
      #if TRT_PHREG_ENABLE_FIRTH
      newlk = f_pen_llik_2(p, newbeta, param);
      #endif
    }

    // check convergence
    fail = std::isnan(newlk) || std::isinf(newlk) == 1;

    if (!fail && halving == 0 && fabs(1 - (loglik/newlk)) < eps) {
      break;
    }

    if (fail || newlk < loglik) { // adjust step size if likelihood decreases
      halving++;
      for (i=0; i<p; i++) {
        newbeta[i] = (beta[i] + newbeta[i])/2;
      }
    } else { // update beta normally
      halving = 0;

      for (i=0; i<p; i++) {
        beta[i] = newbeta[i];
      }
      loglik = newlk;

      if (!firth) u = f_score_2_cpp(p, beta, param);
      else{ 
        #if TRT_PHREG_ENABLE_FIRTH
        u = f_pen_score_2(p, beta, param); 
        #endif
      }
      for (i=0; i<ncolfit; i++) {
        u1[i] = u[colfit[i]];
      }

      imat = f_info_2_cpp(p, beta, param);
      for (i=0; i<ncolfit; i++) {
        for (j=0; j<ncolfit; j++) {
          imat1(i,j) = imat(colfit[i], colfit[j]);
        }
      }

      i = cholesky2_cpp(imat1, ncolfit, toler);

      chsolve2_cpp(imat1, ncolfit, u1);
      std::fill(u.begin(),u.end(),0.0); //u.fill(0);
      for (i=0; i<ncolfit; i++) {
        u[colfit[i]] = u1[i];
      }

      for (i=0; i<p; i++) {
        newbeta[i] = beta[i] + u[i];
      }
    }
  }

  if (iter == maxiter) fail = 1;

  imat = f_info_2_cpp(p, newbeta, param);
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      imat1(i,j) = imat(colfit[i], colfit[j]);
    }
  }

  MatrixRM var1 = invsympd_cpp(imat1, ncolfit, toler);
  MatrixRM var (p,p); //{p,p, std::vector<double>(p*p)};
  for (i=0; i<ncolfit; i++) {
    for (j=0; j<ncolfit; j++) {
      var(colfit[i], colfit[j]) = var1(i,j);
    }
  }

  phregLoopOut out;
  out.coef   = std::move(newbeta);
  out.iter   = iter;
  out.var    = std::move(var);
  out.loglik = newlk;
  out.fail   = fail;
  return out;
}


/*Notes and things to consider:
  - The function is a direct translation of the R function phregcpp
  - The function expects a coxfitout object which contains the necessary data
  - Important notice, from tsesimp it will never pass true to any of the bools below, so a compiler flag was used to disable certain branches that will never be reached.
      These branches will be re-enabled later, but for the purposes of creating a proof of concept quickly this is a shortcut I took. 
*/

// [[Rcpp::export]]
coxfitout phreg_purecpp(
        const coxdata data,
        const std::vector<std::string> covariates,
        const std::string ties = "efron",
        const std::vector<double>& init = {},
        const bool robust = 0,
        const bool est_basehaz = 1,
        const bool est_resid = 1,
        const bool firth = 0,
        const bool plci = 0,
        const double alpha = 0.05,
        const int maxiter = 50,
        const double eps = 1.0e-9) {
          
  int h, i, j, k, n = data.n;
  int p = static_cast<int>(data.p);
  if (p == 0) {
    //error handling needed
  }

  bool has_rep = false;
  std::vector<int> repn(n,1);

  bool has_stratum;
  std::vector<int> stratumn(n,1); // alternativly we could do stratumn(n,1) but filling it later keeps it clearer which is nice
  int n_stratum = static_cast<int>(data.stratum.size());
  // if no stratum
  if (n_stratum == 0) {
    has_stratum = false;
  } else if (n_stratum == n){
    has_stratum = true;
    stratumn = data.stratum;
  } else {
    // n_stratum != n 
    has_stratum = true;
    // pass stratum so this is now redundant 
    // return fail here    
  }

  /*
  // ERROR CHECK before calling phregcpp!
  bool has_time = hasVariable(data, time);
  if (!has_time) stop("data must contain the time variable");
  */
  std::vector<double> timenz = data.time;
  std::vector<double> timen = data.time;
  const int n_time = static_cast<int>(data.time.size());
  if(n_time == 0){
    // fail (or stop :P )
  }
  for(int i = 0; i < n_time; i++){
    const double t = data.time[i];
    if(std::isnan(t) || t < 0.0){
      //return fail
    }
  }
  bool has_time = true;

  std::vector<double> time2n(n,1.0);
  bool has_time2 = !data.time2.empty();
  /*
  // THIS IS TEMPORARILY UNSUPPORTED, to test numerical parity first
  if (has_time2) {
    NumericVector time2nz = data[time2];
    time2n = clone(time2nz);
    if (is_true(any(time2n <= timen))) {
      stop("time2 must be greater than time for each observation");
    }
  }
  bool has_time2 = hasVariable(data, time2);
  */

  bool has_event = !data.event.empty();
  if(!has_event){
    // return fail "data must contain the event variable"
  } 


  std::vector<int> eventnz = data.event;
  std::vector<int> eventn = eventnz;

  int n_event = static_cast<int>(data.event.size());
  int sumEvent = 0;

  for(int i = 0; i < n_event; i++){
    const int e = data.event[i];
    if(e != 1 && e != 0){
      // return fail "event must be 1 or 0 for each observation"
    }
    sumEvent += e;
  }
  if(sumEvent == 0){
    // return fail "at least 1 event is needed to fit the Cox model"
  }

  MatrixRMView zbase{n, p, data.x.data()};
  std::vector<int> row_idx(n);
  std::iota(row_idx.begin(), row_idx.end(), 0);

  RowIndexViewView zn{&zbase, &row_idx};

  bool has_weight = !data.weights.empty();
  std::vector<double> weightn(n,1.0);
  
  /*
  // THIS IS TEMPORARILY UNSUPPORTED, to test numerical parity 
  NumericVector weightn(n, 1.0);
  if (has_weight) {
    NumericVector weightnz = data[weight];
    weightn = clone(weightnz);
    if (is_true(any(weightn <= 0))) {
      stop("weight must be greater than 0");
    }
  }
  */

  bool has_offset = !data.offset.empty();

  std::vector<double> offsetn(n);
  
  /*
  // THIS IS TEMPORARILY UNSUPPORTED, to test numerical parity
  if (has_offset) {
    NumericVector offsetnz = data[offset];
    offsetn = clone(offsetnz);
  }
  */


  bool has_id = !data.id.empty();

  std::vector<int> idn(n,1);
  idn = seq_cpp(1,n);
  /*
  // THIS IS TEMPORARILY SUPPORTED, to test numerical parity
  //
  if (!has_id) {
    idn = seq(1,n);
  } else {
    if (TYPEOF(data[id]) == INTSXP) {
      IntegerVector idv = data[id];
      IntegerVector idwi = unique(idv);
      idwi.sort();
      idn = match(idv, idwi);
    } else if (TYPEOF(data[id]) == REALSXP) {
      NumericVector idv = data[id];
      NumericVector idwn = unique(idv);
      idwn.sort();
      idn = match(idv, idwn);
    } else if (TYPEOF(data[id]) == STRSXP) {
      StringVector idv = data[id];
      StringVector idwc = unique(idv);
      idwc.sort();
      idn = match(idv, idwc);
    } else {
      stop("incorrect type for the id variable in the input data");
    }
  }
    */

  if (robust && has_time2 && !has_id) {
    // fail ("id is needed for counting process data with robust variance");
  }


  std::string meth = ties;
  std::for_each(meth.begin(), meth.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });

  int method = meth == "efron" ? 1 : 0;


  // sort the data by rep
  std::vector<int> order = seq_cpp(0, n-1);
  std::stable_sort(order.begin(), order.end(), [&](int i, int j) {
    return repn[i] < repn[j];
  });

  reorder(repn,order);
  reorder(stratumn,order);
  reorder(timen,order); 
  reorder(time2n,order); 
  reorder(eventn, order); 
  reorder(weightn, order);
  reorder(offsetn,order);
  reorder(idn,order);
  if (p > 0) reorder(row_idx,order); //zn = subset_matrix_by_row(zn, order);

  // exclude observations with missing values
  std::vector<bool> sub(n,1);
  for (i=0; i<n; i++) {
    if ((std::isnan(timen[i])) || (eventn[i] == NA_INTEGER) ||
        (std::isnan(weightn[i])) || (std::isnan(offsetn[i])) 
        ) {
      sub[i] = 0;
    }
    for (j=0; j<p; j++) {
      if (std::isnan(zn(i,j))) sub[i] = 0;
    }
  }

  order = which_true(sub);
  reorder(repn,order);
  reorder(stratumn,order);
  reorder(timen,order);
  reorder(time2n,order);
  reorder(eventn,order);
  reorder(weightn,order);
  reorder(offsetn,order);
  reorder(idn,order);
  if (p > 0) reorder(row_idx, order); //zn = subset_matrix_by_row(zn, order);
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
  std::vector<std::vector<double>> loglik(nreps, std::vector<double>(2));
  std::vector<std::vector<double>> regloglik(nreps, std::vector<double>(2));
  std::vector<double> scoretest(nreps);
  std::vector<int> niter(nreps);
  std::vector<bool> fails(nreps);

  std::vector<int> rep0(nreps*p);
  std::vector<std::string> par0(nreps*p);
  std::vector<double> beta0(nreps*p), sebeta0(nreps*p), rsebeta0(nreps*p,0.0);
  std::vector<std::vector<double>> vbeta0(nreps*p, std::vector<double>(p)), rvbeta0(nreps*p, std::vector<double>(p) );
  std::vector<double> lb0(nreps*p), ub0(nreps*p), prob0(nreps*p);
  std::vector<std::string> clparm0(nreps*p);

  // baseline hazards
  int N = 2*n; // account for additional time 0 rows
  std::vector<int> drep(N, NA_INTEGER), dstratum(N);
  std::vector<double> dtime(N), dnrisk(N), dnevent(N), dncensor(N);
  std::vector<double> dhaz(N), dvarhaz(N);
  std::vector<std::vector<double>> dgradhaz(N, std::vector<double>(p));
  int n0 = 0;
  int bign0 = 0;
  double toler = 1e-12;

  DataFrame basehaz2;
  std::vector<double> resmart(n);

  for (h=0; h<nreps; h++) {
    std::vector<int> q1 = seq_cpp(idx[h], idx[h+1]-1);
    int n1 = static_cast<int>(q1.size());

    std::vector<int> stratum1 = subset_by_idx(stratumn, q1); 
    std::vector<double> time1 = subset_by_idx(timen,q1); 
    std::vector<double> time21 = subset_by_idx(time2n,q1); 
    std::vector<int> event1 = subset_by_idx(eventn,q1); 
    std::vector<double> weight1 = subset_by_idx(weightn,q1); 
    std::vector<double> offset1 = subset_by_idx(offsetn,q1); 
    std::vector<int> id1 = subset_by_idx(idn,q1); 

    MatrixRM z1 (n1,p);

    if (p > 0) {
      for (int ii = 0; ii < n1; ++ii) {
        int src_row = idx[h] + ii; // because q1 = Range(idx[h], idx[h+1]-1)
        for (int jj = 0; jj < p; ++jj) {
          z1(ii, jj) = zn(src_row, jj);
        }
      }
    }

    nobs[h] = n1;
    nevents[h] = sumint_cpp(event1);

    // unify right censored data with counting process data
    std::vector<double> tstart(n1), tstop(n1);
    if (!has_time2) {
      tstop = time1;
    } else {
      tstart = time1;
      tstop = time21;
    }

    // sort by stratum
    std::vector<int> order0 = seq_cpp(0, n1-1);
    std::sort(order0.begin(), order0.end(), [&](int i, int j) {
      return stratum1[i] < stratum1[j];
    });

    std::vector<int> stratum1z = subset_by_idx(stratum1,order0); 
    std::vector<double> tstartz = subset_by_idx(tstart,order0); 
    std::vector<double> tstopz =  subset_by_idx(tstop,order0);
    std::vector<int> event1z = subset_by_idx(event1,order0); 

    // locate the first observation within each stratum
    std::vector<int> istratum(1,0);
    for (i=1; i<n1; i++) {
      if (stratum1z[i] != stratum1z[i-1]) {
        istratum.push_back(i);
      }
    }

    int nstrata = static_cast<int>(istratum.size());
    istratum.push_back(n1);
    // ignore subjects not at risk for any event time
    
    std::vector<int> ignore1z(n1);
    for (int s = 0; s < nstrata; s++) {
      int L = istratum[s];
      int R = istratum[s+1];

      std::vector<double> etime;
      etime.reserve(R - L);

      for (int j = L; j < R; ++j) {
        if (event1z[j] == 1) etime.push_back(tstopz[j]);
      }

      std::sort(etime.begin(), etime.end());
      etime.erase(std::unique(etime.begin(), etime.end()), etime.end());

      if (etime.empty()) {
        for (int j = L; j < R; ++j) ignore1z[j] = 1;
        continue;
      }

      for (int j = L; j < R; ++j) {
        int idx1 = (int)(std::upper_bound(etime.begin(), etime.end(), tstartz[j]) - etime.begin());
        int idx2 = (int)(std::upper_bound(etime.begin(), etime.end(), tstopz[j])  - etime.begin());
        ignore1z[j] = (idx1 == idx2) ? 1 : 0;
      }
    }

    std::vector<int> ignore1(n1);
    for (i=0; i<n1; i++) {
      ignore1[order0[i]] = ignore1z[i];
    }

    int nused = n1 - sumint_cpp(ignore1);

    // sort by stopping time in descending order within each stratum
    std::vector<int> order2 = seq_cpp(0, n1-1);
    std::sort(order2.begin(), order2.end(), [&](int i, int j) {
      return (ignore1[i] < ignore1[j]) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] < stratum1[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] > tstop[j])) ||
        ((ignore1[i] == ignore1[j]) && (stratum1[i] == stratum1[j]) &&
        (tstop[i] == tstop[j]) && (event1[i] < event1[j]));
    });

    std::vector<int> stratum1a = subset_by_idx(stratum1,order2); 
    std::vector<double> tstarta = subset_by_idx(tstart,order2); 
    std::vector<double> tstopa = subset_by_idx(tstop,order2); 
    std::vector<int> event1a = subset_by_idx(event1,order2); 
    std::vector<double> weight1a = subset_by_idx(weight1,order2); 
    std::vector<double> offset1a = subset_by_idx(offset1,order2); 
    std::vector<int> id1a = subset_by_idx(id1,order2);  
    std::vector<int> ignore1a = subset_by_idx(ignore1,order2);  

    //part of subset_matrix replacement
    int start = idx[h];
    std::vector<int> idx_z1a = compose_block_rows(row_idx, start, order2); 
    RowIndexViewView z1a;            // default
    if (p > 0) {
      z1a = RowIndexViewView{ &zbase, &idx_z1a };
    }

    // sort by starting time in descending order within each stratum
    std::vector<int> order1 = seq_cpp(0, n1-1);
    std::sort(order1.begin(), order1.end(), [&](int i, int j) {
      return (ignore1a[i] < ignore1a[j]) ||
        ((ignore1a[i] == ignore1a[j]) && (stratum1a[i] < stratum1a[j])) ||
        ((ignore1a[i] == ignore1a[j]) && (stratum1a[i] == stratum1a[j]) &&
        (tstarta[i] > tstarta[j]));
    });

    coxparams_purecpp param = {nused, &stratum1a, &tstarta, &tstopa, &event1a,
                       &weight1a, &offset1a, z1a, &order1, method};

    std::vector<double> bint(p);

    #if TRT_PHREG_ENABLE_BASEHAZ
    // we disable everything here also because paramx is passed to est_basehaz check only so not needed for tsesimp!
    // prepare the data for estimating the baseline hazards at all time points
    List basehaz1;

    // sort by stopping time in descending order within each stratum
    std::vector<int> order2x = seq_cpp(0, n1-1);
    std::sort(order2x.begin(), order2x.end(), [&](int i, int j) {
      return (stratum1[i] > stratum1[j]) ||
        ((stratum1[i] == stratum1[j]) && (tstop[i] > tstop[j]));
    });

    std::vector<int> stratum1x = subset_by_idx(stratum1,order2x); //stratum1[order2x];
    std::vector<double> tstartx = subset_by_idx(tstart,order2x); //tstart[order2x];
    std::vector<double> tstopx = subset_by_idx(tstop,order2x); //tstop[order2x];
    std::vector<int> event1x = subset_by_idx(event1,order2x); //event1[order2x];
    std::vector<double> weight1x = subset_by_idx(weight1,order2x); //weight1[order2x];
    std::vector<double> offset1x = subset_by_idx(offset1,order2x); //offset1[order2x];
    std::vector<int> id1x = subset_by_idx(id1,order2x); //id1[order2x];
    std::vector<int> ignore1x = subset_by_idx(ignore1,order2x); //ignore1[order2x];
    //NumericMatrix z1x(n1,p);
    //if (p > 0) z1x = subset_matrix_by_row(z1, order2x);
    // need to figure out how to do this part 

    // sort by starting time in descending order within each stratum
    std::vector<int>IntegerVector order1x = seq(0, n1-1);
    std::sort(order1x.begin(), order1x.end(), [&](int i, int j) {
      return (stratum1x[i] > stratum1x[j]) ||
        ((stratum1x[i] == stratum1x[j]) && (tstartx[i] > tstartx[j]));
    });

    coxparams paramx = {n1, stratum1x, tstartx, tstopx, event1x,
                        weight1x, offset1x, z1x, order1x, method};
    #endif 

    std::vector<double> b(p);
    MatrixRM vb (p,p); 
    if (p > 0) {
      phregLoopOut out;
      std::vector<int> colfit = seq_cpp(0,p-1);
      auto init_ok = (init.size() == (size_t)p);
      if (init_ok) {
        for (double v : init) {
          if (!std::isnan(v)) { init_ok = false; break; }
        }
      }
      if (init_ok) {
        out = phregloop_cpp(p, init, &param, maxiter, eps, firth, colfit, p);
      } else {
        out = phregloop_cpp(p, bint, &param, maxiter, eps, firth, colfit, p);
      }
      
      bool fail = out.fail; //["fail"];
      if (fail){
        // throw fail error warning("The algorithm in phregr did not converge");
      } 

      b = out.coef; //["coef"];
      vb = out.var; //as<NumericMatrix>(out["var"]);

      std::vector<double> seb(p);
      for (j=0; j<p; j++) {
        seb[j] = std::sqrt(vb(j,j));
      }

      for (i=0; i<p; i++) {
        rep0[h*p+i] = h+1;
        par0[h*p+i] = covariates[i];
        beta0[h*p+i] = b[i];
        sebeta0[h*p+i] = seb[i];
        for (j=0; j<p; j++) {
          vbeta0[h*p+i][j] = vb(i,j);
        }
      }

      // score statistic
      std::vector<double> scoreint = f_score_2_cpp(p, bint, &param);
      MatrixRM infobint = f_info_2_cpp(p, bint, &param);
      MatrixRM vbint = invsympd_cpp(infobint, p, toler);
      for (i=0; i<p; i++) {
        for (j=0; j<p; j++) {
          scoretest[h] += scoreint[i]*vbint(i,j)*scoreint[j];
        }
      }

      niter[h] = out.iter; //["iter"];
      fails[h] = out.fail; //["fail"];
      
      // robust variance estimates
      //NumericVector rseb(p);  // robust standard error for betahat
      std::vector<double> rseb(p,0.0);
      #if TRT_PHREG_ENABLE_ROBUST
      if (robust) {
        NumericMatrix ressco = f_ressco_2(p, b, &param);

        int nr; // number of rows in the score residual matrix
        if (!has_id) {
          for (i=0; i<n1; i++) {
            for (j=0; j<p; j++) {
              ressco(i,j) = weight1a[i]*ressco(i,j);
            }
          }
          nr = n1;
        } else { // need to sum up score residuals by id
          IntegerVector order = seq(0, n1-1);
          std::sort(order.begin(), order.end(), [&](int i, int j) {
            return id1a[i] < id1a[j];
          });

          IntegerVector id2 = id1a[order];
          IntegerVector idx(1,0);
          for (i=1; i<n1; i++) {
            if (id2[i] != id2[i-1]) {
              idx.push_back(i);
            }
          }

          int nids = static_cast<int>(idx.size());
          idx.push_back(n1);

          NumericVector weight2 = weight1a[order];

          NumericMatrix ressco2(nids,p);
          for (i=0; i<nids; i++) {
            for (j=0; j<p; j++) {
              for (k=idx[i]; k<idx[i+1]; k++) {
                ressco2(i,j) += weight2[k]*ressco(order[k],j);
              }
            }
          }

          ressco = ressco2;  // update the score residuals
          nr = nids;
        }


        NumericMatrix D(nr,p); // DFBETA
        for (i=0; i<nr; i++) {
          for (j=0; j<p; j++) {
            for (k=0; k<p; k++) {
              D(i,j) += ressco(i,k)*vb(k,j);
            }
          }
        }

        NumericMatrix rvb(p,p); // robust variance matrix for betahat
        for (j=0; j<p; j++) {
          for (k=0; k<p; k++) {
            for (i=0; i<nr; i++) {
              rvb(j,k) += D(i,j)*D(i,k);
            }
          }
        }

        for (i=0; i<p; i++) {
          rseb[i] = sqrt(rvb(i,i));
        }

        for (i=0; i<p; i++) {
          rsebeta0[h*p+i] = rseb[i];
          for (j=0; j<p; j++) {
            rvbeta0(h*p+i,j) = rvb(i,j);
          }
        }
      }
      #endif

      // profile likelihood confidence interval for regression coefficients
      std::vector<double> lb(p), ub(p), prob(p);
      std::vector<String> clparm(p);

      //double zcrit = R::qnorm(1-alpha/2,0,1,1,0);
      double zcrit = qnorm_cpp(1.0 - alpha/2.0,true,false);
      if (plci) {
        #if TRT_PHREG_ENABLED_PLCI
        double lmax;
        if (firth) {
          // ALWAYS false bc of the pipeline call
          //lmax = f_pen_llik_2(p, b, &param);
        } else {
          lmax = f_llik_2(p, b, &param);
        }
        double l0 = lmax - 0.5*R::qchisq(1-alpha, 1, 1, 0);

        for (k=0; k<p; k++) {
          lb[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, -1, l0);
          ub[k] = phregplloop(p, b, &param, maxiter, eps, firth, k, 1, l0);

          IntegerVector colfit1(p-1);
          for (i=0; i<k; i++) {
            colfit1[i] = i;
          }
          for (i=k+1; i<p; i++) {
            colfit1[i-1] = i;
          }

          NumericVector b0(p);
          List out0 = phregloop(p, b0, &param, maxiter, eps, firth,
                                colfit1, p-1);
          double lmax0 = out0["loglik"];
          prob[k] = R::pchisq(-2*(lmax0 - lmax), 1, 0, 0);
          clparm[k] = "PL";
        }
        #endif
      } if(!plci) {
        //formally an else statement
        for (k=0; k<p; k++) {
          if (!robust) {
            lb[k] = b[k] - zcrit*seb[k];
            ub[k] = b[k] + zcrit*seb[k];
            double z = b[k]/seb[k];
            prob[k] = pchisq_cpp(z*z, 1.0, false, false); 
          } else {
            lb[k] = b[k] - zcrit*rseb[k];
            ub[k] = b[k] + zcrit*rseb[k];
            double z = b[k]/rseb[k];
            prob[k] = pchisq_cpp(z*z, 1.0, false, false); 
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
    }

    // log-likelihoods
    #if TRT_PHREG_ENABLE_FIRTH
    if (firth) {
      loglik(h,0) = f_pen_llik_2(p, bint, &param);
      loglik(h,1) = f_pen_llik_2(p, b, &param);

      regloglik(h,0) = f_llik_2(p, bint, &param);
      regloglik(h,1) = f_llik_2(p, b, &param);
    } else {
      loglik(h,0) = f_llik_2(p, bint, &param);
      loglik(h,1) = f_llik_2(p, b, &param);
    }
    #endif


    // estimate baseline hazard
    #if TRT_PHREG_ENABLE_BASEHAZ
    if (est_basehaz) {
      basehaz1 = f_basehaz(p, b, &paramx);

      IntegerVector dstratum1 = basehaz1["stratum"];
      NumericVector dtime1 = basehaz1["time"];
      NumericVector dnrisk1 = basehaz1["nrisk"];
      NumericVector dnevent1 = basehaz1["nevent"];
      NumericVector dncensor1 = basehaz1["ncensor"];
      NumericVector dhaz1 = basehaz1["haz"];
      NumericVector dvarhaz1 = basehaz1["varhaz"];
      int J = static_cast<int>(dstratum1.size());

      // add to output data frame
      for (j=0; j<J; j++) {
        drep[n0+j] = h+1;
        dstratum[n0+j] = dstratum1[j];
        dtime[n0+j] = dtime1[j];
        dnrisk[n0+j] = dnrisk1[j];
        dnevent[n0+j] = dnevent1[j];
        dncensor[n0+j] = dncensor1[j];
        dhaz[n0+j] = dhaz1[j];
        dvarhaz[n0+j] = dvarhaz1[j];

        if (p > 0) {
          NumericMatrix dgradhaz1 = basehaz1["gradhaz"];
          for (i=0; i<p; i++) {
            dgradhaz(n0+j,i) = dgradhaz1(j,i);
          }
        }
      }

      n0 += J;
    }
    #endif

    // martingale residuals
    #if TRT_PHREG_ENABLE_RESID
    if (est_resid) {
      NumericVector resid = f_resmart(p, b, &param);

      for (i=0; i<n1; i++) {
        resmart[bign0 + order2[i]] = resid[i];
      }

      bign0 += n1;
    }
    #endif
  }

  #if TRT_PHREG_ENABLE_BASEHAZ
  if (est_basehaz) {
    IntegerVector sub = which(!is_na(drep));
    drep = drep[sub];
    dstratum = dstratum[sub];
    dtime = dtime[sub];
    dnrisk = dnrisk[sub];
    dnevent = dnevent[sub];
    dncensor = dncensor[sub];
    dhaz = dhaz[sub];
    dvarhaz = dvarhaz[sub];
    if (p > 0) dgradhaz = subset_matrix_by_row(dgradhaz, sub);
  }
  #endif
  // returns results after this, previously had code to turn vectors/matrixes back into R objects but not needed anymore
  
  coxfitout output;
  output.fail = false; // should not have failed at this point, (when fail conditions are met it should instead instantiate coxfitout type and returning it w result.fail = true)
  output.beta = beta0;
  output.sebeta = sebeta0;
  output.p = prob0;
  return output;

  //add any additional parameters below that aren't hard required

}