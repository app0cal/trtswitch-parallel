#include <vector>
#include <string>
#include "utilities_pure.h"

// tsesimp structs to be used

// results from each bootstrap function call
struct RepResult {
  double psi0hat = NA_REAL, psi0lower = NA_REAL, psi0upper = NA_REAL;
  double psi1hat = NA_REAL, psi1lower = NA_REAL, psi1upper = NA_REAL;
  double hrhat   = NA_REAL, hrlower   = NA_REAL, hrupper   = NA_REAL;
  double pvalue  = NA_REAL;
  int fail = 0;
};

// lifereg structs
 
// input for lifereg_purecpp
struct trial_data {
  std::vector<double>              pps; // time/nickname for time
  std::vector<int>                 event; // event :P
  std::vector<int>                 swtrt; // switch treatment
  std::vector<std::string>         aft_names; // covariates for the outcome model
  const MatrixRM* aft;  // columns matching aft names
  const std::vector<int>* order_pp; 

  //bottom below are optional for survival analysis parts
  std::vector<std::string>        rep; // never passed
  std::vector<double>             stratum; // never passed
  std::vector<double>             time2; // time variable
  std::vector<double>             weight; // weights for the Cox model
  std::vector<double>             offset; // covariates for the Cox model
  std::vector<std::string>       id_raw; // numeric id variable
};

// output for liferegloop()
struct liferegloopresult{
  std::vector<double> coef;
  int iter;
  std::vector<std::vector<double>> var;
  double loglik;
  bool fail;
};

//output struct for lifereg_purecpp
struct liferegOut{
  bool fail;

  // basically all we need for our purposes xd 

  // sum stat
  //add rest if needed

  // parest
  std::vector<double> beta; //beta0
  std::vector<double> sebeta; //rsebeta0
  // add rest if needed
};

// phreg structs

// input for phreg_purecpp
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

struct coxparams_purecpp {
  int n = 0;

  const std::vector<int>* stratum = nullptr;
  const std::vector<double>* tstart = nullptr;
  const std::vector<double>* tstop  = nullptr;
  const std::vector<int>* event     = nullptr;
  const std::vector<double>* weight = nullptr;
  const std::vector<double>* offset = nullptr;

  // This is your “X matrix” view (row-reordered)
  RowIndexViewView z;

  // The row traversal order used by the math routines (e.g., order1/order1x)
  const std::vector<int>* order = nullptr;

  int method = 0; // ties: 0=breslow, 1=efron 
};

struct PhregFit {
  std::vector<double> coef;   // length p
  int iter = 0;
  std::vector<double> var;    // row-major p*p
  double loglik = NAN;
  bool fail = false;
};

struct phregLoopOut {
  std::vector<double> coef;  // size p
  int iter = 0;
  MatrixRM var;              // p x p
  double loglik = std::numeric_limits<double>::quiet_NaN();
  bool fail = false;
};

struct coxfitout {
  bool fail;
  std::vector<double> beta;     // size p (here p=1 for treated, or more if you include more covars)
  std::vector<double> sebeta;   // size p
  std::vector<double> p;        // size p (optional)
};

// list of functions and other struct headings below

/*
// Note: Planning to replace this because this does way more work than needed when we can use pointers for a flat array and pass a view or just & 
Holds exactly what bygroup(...) returned in Rcpp:
//  - index[i] = group ID (1…G) for row i
//  - lookup[g] = vector of all row-indices i where index[i]==g
*/
struct group_index {
  std::vector<int> index;
  std::vector<std::vector<int>> lookup;
};


group_index group_by(
  const std::vector<std::string>& factor
);

double f_llik_2_cpp(
  int p, 
  const std::vector<double>& par, 
  const void *ex
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

double liferegplloop_cpp(
  int p, std::vector<double> par, 
  void *ex,
  int maxiter, 
  double eps,
  int k, 
  int which, 
  double l0
);

liferegOut lifereg_purecpp(
  //const DataFrame data expects trial_data with ...
  const trial_data& data,
  const std::vector<std::string>& covariates,
  const std::string dist,
  const std::vector<double>& init,
  const bool robust,
  const bool plci,
  const double alpha,
  const int maxiter,
  const double eps
);

coxfitout phreg_purecpp(
  const coxdata data,
  const std::vector<std::string> covariates,
  const std::string ties,
  const std::vector<double>& init,
  const bool robust,
  const bool est_basehaz,
  const bool est_resid,
  const bool firth,
  const bool plci,
  const double alpha,
  const int maxiter,
  const double eps 
);

template<typename STL, typename RCPPTYPE>
STL to_std(const RCPPTYPE& rvec) {
  return Rcpp::as<STL>(rvec);
}