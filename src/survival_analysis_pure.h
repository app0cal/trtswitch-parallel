#include <vector>
#include <string>
#include "utilities_pure.h"

//structs to tightly connect every process in this function, to extend this or enable other features simply edit these or 
//what we get sent from out test dummy of tsesimp_test.cpp
struct trial_data {
  std::vector<double>              pps; // time
  std::vector<int>                 event; // event :P
  std::vector<int>                 swtrt; // switch treatment
  std::vector<std::string>         aft_names; // covariates for the outcome model
  std::vector<std::vector<double>> aft;  // columns matching aft names

  //bottom below are optional for survival analysis parts
  std::vector<double>             time2; // time variable
  std::vector<double>             weight; // weights for the Cox model
  std::vector<double>             offset; // covariates for the Cox model
  std::vector<std::string>       id_raw; // numeric id variable
};

/*
Holds exactly what bygroup(...) returned in Rcpp:
//  - index[i] = group ID (1…G) for row i
//  - lookup[g] = vector of all row-indices i where index[i]==g

*/
struct group_index {
  std::vector<int> index;
  std::vector<std::vector<int>> lookup;
};

struct aftparams_pure {
  std::string dist;
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<int> status;
  std::vector<double> weight;
  std::vector<double> offset;
  int nstrata;
  //std::vector<std::vector<double>> z;
  // replaced w experimental view for better performance
  const MatrixRM* zbase = nullptr;
  const std::vector<int>* rows = nullptr; // mapping
  int nrows = 0, ncols = 0;

  double z(int i, int j) const {
    return (*zbase)((*rows)[i], j);
  }
  
  
};

struct liferegloopresult{
  std::vector<double> coef;
  int iter;
  std::vector<std::vector<double>> var;
  double loglik;
  bool fail;
};

struct f_der_eta_1_result {
  std::vector<double> dg;
  std::vector<double> ddg;
};


/*Function Signatures:

Notes:
A lot of these functions are easy to write in for loops but just to make code more readable and maintainable, these functions are introduced.

*/


//this function is kind of not used so i might remove when finalizing
/*
template<typename STL, typename RCPPTYPE>
STL to_std(
  const RCPPTYPE& rvec
);
*/

group_index group_by(
  const std::vector<std::string>& factor
);

bool has_col_aft(const trial_data& d, const std::string& col_name);

void safety_check_col_aft(const trial_data& d);

bool has_variable(
  const trial_data& d, 
  const std::string& var_name
);

std::vector<double> f_score_1_cpp(
  int p, 
  std::vector<double> par, 
  void *ex
);

double f_llik_1_cpp(
  int p, 
  std::vector<double> par, 
  void *ex
);

std::vector<std::vector<double>> f_info_1_cpp(
  int p, 
  std::vector<double> par, 
  void *ex
);

std::vector<std::vector<double>> f_ressco_1_cpp(
  int p, 
  std::vector<double> par, 
  void *ex
);

std::vector<std::vector<double>> f_jj_1_cpp(
  int p, 
  std::vector<double> par, 
  void *ex
);

liferegloopresult liferegloop_cpp(
  int p, std::vector<double> par, 
  void *ex,
  int maxiter, 
  double eps,
  std::vector<int> colfit, 
  int ncolfit
);

f_der_eta_1_result f_der_eta_1_cpp(
  std::vector<double> eta, 
  std::vector<double> sig, 
  void *ex
);

double liferegplloop_cpp(
  int p, std::vector<double> par, 
  void *ex,
  int maxiter, 
  double eps,
  int k, 
  int which, 
  double l0
);

List lifereg_purecpp(
  //const DataFrame data expects trial_data with ...
  const trial_data data,
  const std::vector<std::string>& rep,
  const std::vector<std::string>& stratum,
  const std::string time,
  const std::string time2,
  const std::string event,
  const std::vector<std::string>& covariates,
  const std::string weight,
  const std::string offset,
  const std::string id,
  const std::string dist,
  const std::vector<double>& init,
  const bool robust,
  const bool plci,
  const double alpha,
  const int maxiter,
  const double eps
);

template<typename STL, typename RCPPTYPE>
STL to_std(const RCPPTYPE& rvec) {
  return Rcpp::as<STL>(rvec);
}