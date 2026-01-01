#include <vector>
#include <string>

// tsesimp functions below: 

//structs to tightly connect every process in this function, to extend this or enable other features simply edit these or 
//what we get sent from out test dummy of tsesimp_test.cpp
struct trial_data {
  std::vector<double>              pps; // time
  std::vector<int>                 event; // event 
  std::vector<int>                 swtrt; // switch treatment
  std::vector<std::string>         aft_names; // covariates for the outcome model
  std::vector<std::vector<double>> aft;  // columns matching aft names

  //bottom below are optional for survival analysis parts
  std::vector<double>             time2; // time variable
  std::vector<double>             weight; // weights for the Cox model
  std::vector<double>             offset; // covariates for the Cox model
  std::vector<std::string>       id_raw; // numeric id variable
};


// survival_analysis structs below:

struct aftparams_pure {
  std::string dist;
  std::vector<int> strata;
  std::vector<double> tstart;
  std::vector<double> tstop;
  std::vector<int> status;
  std::vector<double> weight;
  std::vector<double> offset;
  std::vector<std::vector<double>> z;
  int nstrata;
};

struct coxdata {
  std::vector<double> time;     // t_star
  std::vector<int> event;       // d_star
  std::vector<int> stratum;     // length n,nullptr if none
  std::vector<double> x;        // row-major n*p 
  int n;
  int p;

  //OPTIONAL
  std::vector<double> time2;     // empty -> none
  std::vector<double> weights; // empty -> none
  std::vector<double> offset; // empty -> none
  std::vector<int> id; // empty -> none
};

struct coxfitout {
  bool fail;
  std::vector<double> beta;     // size p (here p=1 for treated, or more if you include more covars)
  std::vector<double> sebeta;   // size p
  std::vector<double> p;        // size p (optional)
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

// helper functions below:

/*
Holds exactly what bygroup(...) returned in Rcpp:
//  - index[i] = group ID (1…G) for row i
//  - lookup[g] = vector of all row-indices i where index[i]==g

*/
struct group_index {
  std::vector<int> index;             //n, 1 .. G 
  std::vector<std::string> levels;    //G+1, levels[g] label
  std::vector<std::vector<int>> rows; //optional
};

struct group_lookup_table {
  std::vector<int> index; // length n
  // one entry per group, each group row holds the stratum columns as strings
  std::vector<std::vector<std::string>> group_values; // [g][col]
  // or: group_values[g] is a vector<string> of length p_stratum
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

List phreg_purecpp(
const DataFrame data,
const StringVector& rep = "",
const StringVector& stratum = "",
const std::string time = "time",
const std::string time2 = "",
const std::string event = "event",
const StringVector& covariates = "",
const std::string weight = "",
const std::string offset = "",
const std::string id = "",
const std::string ties = "efron",
const NumericVector& init = NA_REAL,
const bool robust = 0,
const bool est_basehaz = 1,
const bool est_resid = 1,
const bool firth = 0,
const bool plci = 0,
const double alpha = 0.05,
const int maxiter = 50,
const double eps = 1.0e-9);

template<typename STL, typename RCPPTYPE>
STL to_std(const RCPPTYPE& rvec) {
  return Rcpp::as<STL>(rvec);
}