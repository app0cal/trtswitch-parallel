// [[Rcpp::plugins(cpp14)]]
// [[asdRcpp::plugins(openmp)]]
//#error "COMPILE CHECK: tsesimp_mt.cpp is being compiled"

#include "utilities.h"
#include "survival_analysis.h"

using namespace Rcpp;

//add custom includes
#include "utilities_pure.h"
#include "survival_analysis_pure.h"
#include <omp.h>
#include <R_ext/Print.h>


struct sharedInputs{
  int n, q, p, p2; 
  std::vector<std::string> covariates;
  std::vector<std::string> covariates_aft;
  std::vector<std::string> cov_names_cox;
  std::string dist;
  bool recensor, swtrt_control_only; 
  double alpha, zcrit; 
  std::string ties; 
  double offset;
};

struct replicateInputs{ 
  std::vector<int> &idb, &stratumb, &eventb, &treatb;
  std::vector<int> &pdb, &swtrtb;
  std::vector<double> &timeb, &censor_timeb, &pd_timeb, &swtrt_timeb;
  const MatrixRM &zn, &zn_aft;
  std::vector<int> &row;
};

// pure cpp version of lambda f() only computes core outputs
// parameters for function are based off lambda capture list + parameter list below
RepResult run_replicate_cpp(  
  int k,  // current index 
  // const refs to base data + boot indices + config 
  const sharedInputs& shared,
  replicateInputs& in){

  // AFTER NUMERICAL PARITY CHECKS HERE REMOVE AUTO KEYWORDS, auto makes it much more work to debug so we can just swap back to regular types
  // shared alias
  auto& covariates_aft = shared.covariates_aft;
  auto& n = shared.n, & q = shared.q, &p = shared.p, &p2 = shared.p2;
  auto& swtrt_control_only = shared.swtrt_control_only;
  auto& recensor = shared.recensor;
  auto& offset = shared.offset;
  auto& alpha = shared.alpha, &zcrit = shared.zcrit;
  std::string ties = shared.ties, dist = shared.dist;
  std::vector<std::string> cov_names_cox = shared.cov_names_cox;

  // replicate alias
  auto& idb = in.idb;
  auto& stratumb = in.stratumb;
  auto& timeb = in.timeb;
  auto& eventb = in.eventb;
  auto& treatb = in.treatb;
  auto& censor_timeb = in.censor_timeb;
  auto& pdb = in.pdb;
  auto& pd_timeb = in.pd_timeb;
  auto& swtrtb = in.swtrtb;
  auto& swtrt_timeb = in.swtrt_timeb;
  auto& zb = in.zn;
  auto& zb_aft = in.zn_aft;
  auto& row = in.row;
  
  int h, i, j;
  bool fail = false; // whether any model fails to converge
  //NumericVector init(1, NA_REAL);
  std::vector<double> init_cpp(1,std::numeric_limits<double>::quiet_NaN()); 

  //List fit1 = List::create(Named("junk") = "dog");
  liferegOut fit1_cpp; // native cpp version
  //DataFrame data1;
  //List fit_outcome = List::create(Named("junk2") = "dog2");
  coxfitout fit_outcome_cpp; // native cpp version

  bool uselifepure = true;
  bool usephregpure = true;

  // order data by treat
  //IntegerVector order = seq(0, n-1);
  std::vector<int> ordercpp = seq_cpp(0,n-1);
  std::sort(ordercpp.begin(), ordercpp.end(), [&](int i, int j) {
    return (treatb[i] < treatb[j]);
  });
  //std::vector<int> ordercpp = seq_cpp(0,n-1);

  // checks if row and order passed are identity (so no need to reorder)
  bool is_identity = true;
  for (int i = 0; i < n; ++i) {
    if (row[i] != i || ordercpp[i] != i){ 
      is_identity = false; 
      break; 
    }
  }
  // when called from the bootstrap check should auto sort but we have this as a backup just in case!
  if (!is_identity) {
    //Rcpp::Rcout << "identity fail" << std::endl;
    // reorder vectors AND reorder matrices (see 3B)
    reorder(idb,ordercpp); //idb = idb[order];
    reorder(stratumb,ordercpp); //stratumb = stratumb[order];
    reorder(timeb,ordercpp); //timeb = timeb[order];
    reorder(eventb,ordercpp); //eventb = eventb[order];
    reorder(treatb,ordercpp); //= treatb[order];
    reorder(censor_timeb,ordercpp); //  censor_timeb = censor_timeb[order];
    reorder(pdb,ordercpp);  //  pdb = pdb[order];
    reorder(pd_timeb,ordercpp); // = pd_timeb[order];
    reorder(swtrtb,ordercpp); //swtrtb = swtrtb[order];
    reorder(swtrt_timeb,ordercpp); // = swtrt_timeb[order];
    reorder(row,ordercpp);
  }

  //zb = subset_matrix_by_row(zb, order);
  //zb_aft = subset_matrix_by_row(zb_aft, order);

  // time and event adjusted for treatment switching
  std::vector<double> t_star = timeb; //clone(timeb);
  std::vector<int> d_star = eventb; //clone(eventb);

  double psi0hat = 0, psi0lower = 0, psi0upper = 0;
  double psi1hat = 0, psi1lower = 0, psi1upper = 0;

  //DataFrame data_outcome;

  // treat arms that include patients who switched treatment
  //IntegerVector treats(1);
  //treats[0] = 0;
  std::vector<int> treats(1,0);
  if (!swtrt_control_only) {
    treats.push_back(1);
  }
  
  int K = static_cast<int>(treats.size());
  for (h=0; h<K; h++) {
    // post progression data
    std::vector<int> l;
    l.reserve(treatb.size());
    for(int i = 0; i < (int)treatb.size(); i++){
      // change treatb[i] == h to treatb[i] == treats[h] later if needed
      if(treatb[i] == h && pdb[i] == 1){
        l.push_back(i);
      }
    }
    int ls = (int)l.size();
    std::vector<int> rows_pp; rows_pp.reserve(ls);
    for (int ii = 0; ii < l.size(); ++ii) { 
      //1 index so we are subtracting -1 for now 
      rows_pp.push_back(row[l[ii]]);  // map logical index -> base row
    }
    std::vector<int> id2;       id2.reserve(ls);
    std::vector<double> time2;  time2.reserve(ls);
    std::vector<int> event2;    event2.reserve(ls);
    std::vector<int> swtrt2;    swtrt2.reserve(ls);
    for(int i : l){
      id2.push_back(idb[i]);
      time2.push_back(timeb[i] - pd_timeb[i] + offset);
      event2.push_back(eventb[i]);
      swtrt2.push_back(swtrtb[i]);
    }
 
    // rcpp containers -> stl container -> rcpp List
    trial_data td;
    td.pps = time2;
    td.event = event2;
    td.swtrt = swtrt2;
    td.aft = &zb_aft;
    td.order_pp = &rows_pp;

    //Rcpp::Rcout << "[tsesimp_mt] line 853 before lifereg call" << std::endl;
    fit1_cpp = lifereg_purecpp(
      td, covariates_aft, dist, 
      init_cpp, 0, 0, alpha, 50, 1.0e-9);
    //Rcpp::Rcout << "[tsesimp_mt] line 857 after lifereg call" << std::endl;
    
    bool fail1 = fit1_cpp.fail; //sumstat1["fail"];
    if (fail1 == true) fail = true;
    
    //DataFrame parest1 = DataFrame(fit1["parest"]);
    std::vector<double> beta1 = fit1_cpp.beta; //["beta"];
    std::vector<double> sebeta1 = fit1_cpp.sebeta; //["sebeta"];

    //int isw = k;
    double psihat = -beta1[1];
    double psilower = -(beta1[1] + zcrit*sebeta1[1]);
    double psiupper = -(beta1[1] - zcrit*sebeta1[1]);
    
    // calculate counter-factual survival times
    double a = std::exp(psihat);
    double c0 = std::min(1.0, a);
    for (i=0; i<n; i++) {
      if (treatb[i] == h) {
        double b2, u_star, c_star;
        if (swtrtb[i] == 1) {
          b2 = pdb[i] == 1 ? pd_timeb[i] : swtrt_timeb[i];
          b2 = b2 - offset;
          u_star = b2 + (timeb[i] - b2)*a;
        } else {
          u_star = timeb[i];
        }
        
        if (recensor) {
          c_star = censor_timeb[i]*c0;
          t_star[i] = std::min(u_star, c_star);
          d_star[i] = c_star < u_star ? 0 : eventb[i];
        } else {
          t_star[i] = u_star;
          d_star[i] = eventb[i];
        }
      }
    }
    
    // update treatment-specific causal parameter estimates
    if (h == 0) {
      psi0hat = psihat;
      psi0lower = psilower;
      psi0upper = psiupper;
    } else {
      psi1hat = psihat;
      psi1lower = psilower;
      psi1upper = psiupper;
    }
  }

  //
  coxdata data_outcome;
  //List fit_outcome;
  data_outcome.n = t_star.size();
  data_outcome.p = 1 + zb.ncols(); // + 1 bc it gives n of columns and we want to iterate over p 

  data_outcome.time = t_star;
  data_outcome.event = d_star;
  data_outcome.stratum = stratumb;
  // Pack x row-major: [treated | zb]
  data_outcome.x.resize((size_t)data_outcome.n * (size_t)data_outcome.p);

  for (int i = 0; i < data_outcome.n; ++i) {
    data_outcome.x[(size_t)i * data_outcome.p + 0] = (double)treatb[i];   // treated is column 0
    for (int j = 0; j < p; ++j) {
      data_outcome.x[(size_t)i * data_outcome.p + (size_t)(1 + j)] = zb(row[i], j);
    }
  }
  
  fit_outcome_cpp = phreg_purecpp(data_outcome,cov_names_cox,ties,init_cpp,
                              0,0,0,0,0,alpha,50,1.0e-9);

  

  //DataFrame sumstat_cox = DataFrame(fit_outcome["sumstat"]);
  bool fail_cox = fit_outcome_cpp.fail; //["fail"];
  if (fail_cox == 1) fail = 1;

  //DataFrame parest = DataFrame(fit_outcome["parest"]);
  std::vector<double> beta = fit_outcome_cpp.beta; //["beta"];
  std::vector<double> sebeta = fit_outcome_cpp.sebeta; //parest["sebeta"]; 
  std::vector<double> pval = fit_outcome_cpp.p; //parest["p"];
  double hrhat = exp(beta[0]);
  double hrlower = exp(beta[0] - zcrit*sebeta[0]);
  double hrupper = exp(beta[0] + zcrit*sebeta[0]);
  double pvalue = pval[0];

  // CORE RESULTS FOR run_replicate in bootstrap loop
  RepResult result;
  result.fail = fail;
  result.hrhat = hrhat;
  result.psi0hat = psi0hat;
  result.psi1hat = psi1hat;

  // any additional outputs, make changes to the RepResult struct and add below :)


  return result;
};



// [[Rcpp::export]]
List tsesimpcpp_mt(const DataFrame data,
                const std::string id = "id",
                const StringVector& stratum = "",
                const std::string time = "time",
                const std::string event = "event",
                const std::string treat = "treat",
                const std::string censor_time = "censor_time",
                const std::string pd = "pd",
                const std::string pd_time = "pd_time",
                const std::string swtrt = "swtrt",
                const std::string swtrt_time = "swtrt_time",
                const StringVector& base_cov = "",
                const StringVector& base2_cov = "",
                const std::string aft_dist = "weibull",
                const bool strata_main_effect_only = 1,
                const bool recensor = 1,
                const bool admin_recensor_only = 1,
                const bool swtrt_control_only = 1,
                const double alpha = 0.05,
                const std::string ties = "efron",
                const double offset = 1,
                const bool boot = 1,
                const int n_boot = 1000,
                const int seed = NA_INTEGER) {
  //Rcpp::Rcout << "[tsesimp] line 313 start!" << std::endl;
  Rcpp::RNGScope scope;
  Rcpp::Rcout << "omp max threads: " << omp_get_max_threads() << "\n";
  int i, j, k, n = data.nrow();
  
  int p = static_cast<int>(base_cov.size());
  if (p == 1 && (base_cov[0] == "" || base_cov[0] == "none")) p = 0;
  
  int p2 = static_cast<int>(base2_cov.size());
  if (p2 == 1 && (base2_cov[0] == "" || base2_cov[0] == "none")) p2 = 0;
  
  int p_stratum = static_cast<int>(stratum.size());
  
  bool has_stratum;
  IntegerVector stratumn(n);
  DataFrame u_stratum;
  IntegerVector d(p_stratum);
  IntegerMatrix stratan(n,p_stratum);
  if (p_stratum == 1 && (stratum[0] == "" || stratum[0] == "none")) {
    has_stratum = 0;
    stratumn.fill(1);
    d[0] = 1;
    stratan(_,0) = stratumn;
  } else {
    List out = bygroup(data, stratum);
    has_stratum = 1;
    stratumn = out["index"];
    u_stratum = DataFrame(out["lookup"]);
    d = out["nlevels"];
    stratan = as<IntegerMatrix>(out["indices"]);
  }
  
  IntegerVector stratumn_unique = unique(stratumn);
  int nstrata = static_cast<int>(stratumn_unique.size());
  
  bool has_id = hasVariable(data, id);
  bool has_time = hasVariable(data, time);
  bool has_event = hasVariable(data, event);
  bool has_treat = hasVariable(data, treat);
  bool has_censor_time = hasVariable(data, censor_time);
  bool has_pd = hasVariable(data, pd);
  bool has_pd_time = hasVariable(data, pd_time);
  bool has_swtrt = hasVariable(data, swtrt);
  bool has_swtrt_time = hasVariable(data, swtrt_time);
  
  if (!has_id) {
    stop("data must contain the id variable");
  }
  
  IntegerVector idn(n);
  IntegerVector idwi;
  NumericVector idwn;
  StringVector idwc;
  if (TYPEOF(data[id]) == INTSXP) {
    IntegerVector idv = data[id];
    idwi = unique(idv);
    idwi.sort();
    idn = match(idv, idwi);
  } else if (TYPEOF(data[id]) == REALSXP) {
    NumericVector idv = data[id];
    idwn = unique(idv);
    idwn.sort();
    idn = match(idv, idwn);
  } else if (TYPEOF(data[id]) == STRSXP) {
    StringVector idv = data[id];
    idwc = unique(idv);
    idwc.sort();
    idn = match(idv, idwc);
  } else {
    stop("incorrect type for the id variable in the input data");
  }
  
  if (!has_time) {
    stop("data must contain the time variable");
  }
  
  if (TYPEOF(data[time]) != INTSXP && TYPEOF(data[time]) != REALSXP) {
    stop("time must take numeric values");
  }
  
  NumericVector timenz = data[time];
  NumericVector timen = clone(timenz);
  if (is_true(any(timen <= 0.0))) {
    stop("time must be positive");
  }
  
  if (!has_event) {
    stop("data must contain the event variable");
  }
  
  if (TYPEOF(data[event]) != INTSXP && TYPEOF(data[event]) != LGLSXP) {
    stop("event must take integer or logical values");
  }
  
  IntegerVector eventnz = data[event];
  IntegerVector eventn = clone(eventnz);
  if (is_true(any((eventn != 1) & (eventn != 0)))) {
    stop("event must be 1 or 0");
  }
  
  if (is_true(all(eventn == 0))) {
    stop("at least 1 event is needed");
  }
  
  if (!has_treat) {
    stop("data must contain the treat variable");
  }
  
  // create the numeric treat variable
  if (!has_treat) stop("data must contain the treat variable");
  IntegerVector treatn(n);
  IntegerVector treatwi;
  NumericVector treatwn;
  StringVector treatwc;
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    IntegerVector treatv = data[treat];
    treatwi = unique(treatv);
    if (treatwi.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwi == 0) | (treatwi == 1)))) {
      treatwi = IntegerVector::create(1,0);
      treatn = 2 - treatv;
    } else {
      treatwi.sort();
      treatn = match(treatv, treatwi);
    }
  } else if (TYPEOF(data[treat]) == REALSXP) {
    NumericVector treatv = data[treat];
    treatwn = unique(treatv);
    if (treatwn.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    // special handling for 1/0 treatment coding
    if (is_true(all((treatwn == 0) | (treatwn == 1)))) {
      treatwn = NumericVector::create(1,0);
      treatn = 2 - as<IntegerVector>(treatv);
    } else {
      treatwn.sort();
      treatn = match(treatv, treatwn);
    }
  } else if (TYPEOF(data[treat]) == STRSXP) {
    StringVector treatv = data[treat];
    treatwc = unique(treatv);
    if (treatwc.size() != 2) {
      stop("treat must have two and only two distinct values");
    }
    treatwc.sort();
    treatn = match(treatv, treatwc);
  } else {
    stop("incorrect type for the treat variable in the input data");
  }
  
  treatn = 2 - treatn; // use the 1/0 treatment coding
  
  if (!has_censor_time) {
    stop("data must contain the censor_time variable");
  }
  
  if (TYPEOF(data[censor_time]) != INTSXP &&
      TYPEOF(data[censor_time]) != REALSXP) {
    stop("censor_time must take numeric values");
  }
  
  NumericVector censor_timenz = data[censor_time];
  NumericVector censor_timen = clone(censor_timenz);
  if (is_true(any(censor_timen < timen))) {
    stop("censor_time must be greater than or equal to time");
  }
  
  if (!admin_recensor_only) {
    for (i=0; i<n; i++) {
      if (eventn[i] == 0) { // use the actual censoring time for dropouts
        censor_timen[i] = timen[i];
      }
    }
  }
  
  if (!has_pd) {
    stop("data must contain the pd variable");
  }
  
  if (TYPEOF(data[pd]) != INTSXP && TYPEOF(data[pd]) != LGLSXP) {
    stop("pd must take integer or logical values");
  }
  
  IntegerVector pdnz = data[pd];
  IntegerVector pdn = clone(pdnz);
  if (is_true(any((pdn != 1) & (pdn != 0)))) {
    stop("pd must be 1 or 0");
  }
  
  if (!has_pd_time) {
    stop("data must contain the pd_time variable");
  }
  
  if (TYPEOF(data[pd_time]) != INTSXP && TYPEOF(data[pd_time]) != REALSXP) {
    stop("pd_time must take numeric values");
  }
  
  NumericVector pd_timenz = data[pd_time];
  NumericVector pd_timen = clone(pd_timenz);
  for (i=0; i<n; i++) {
    if (pdn[i] == 1 && pd_timen[i] < 0.0) {
      stop("pd_time must be nonnegative");
    }
    
    if (pdn[i] == 1 && pd_timen[i] > timen[i]) {
      stop("pd_time must be less than or equal to time");
    }
  }
  
  if (!has_swtrt) {
    stop("data must contain the swtrt variable");
  }
  
  if (TYPEOF(data[swtrt]) != INTSXP && TYPEOF(data[swtrt]) != LGLSXP) {
    stop("swtrt must take integer or logical values");
  }
  
  IntegerVector swtrtnz = data[swtrt];
  IntegerVector swtrtn = clone(swtrtnz);
  if (is_true(any((swtrtn != 1) & (swtrtn != 0)))) {
    stop("swtrt must be 1 or 0");
  }
  
  if (is_false(any((pdn == 1) & (swtrtn == 1) & (treatn == 0)))) {
    stop("at least 1 pd and swtrt is needed in the control group");
  }
  
  if (!swtrt_control_only) {
    if (is_false(any((pdn == 1) & (swtrtn == 1) & (treatn == 1)))) {
      stop("at least 1 pd and swtrt is needed in the treatment group");
    }
  }
  
  if (!has_swtrt_time) {
    stop("data must contain the swtrt_time variable");
  }
  
  if (TYPEOF(data[swtrt_time]) != INTSXP &&
      TYPEOF(data[swtrt_time]) != REALSXP) {
    stop("swtrt_time must take numeric values");
  }
  
  NumericVector swtrt_timenz = data[swtrt_time];
  NumericVector swtrt_timen = clone(swtrt_timenz);
  for (i=0; i<n; i++) {
    if (swtrtn[i] == 1 && swtrt_timen[i] < 0.0) {
      stop("swtrt_time must be nonnegative");
    }
    
    if (swtrtn[i] == 1 && swtrt_timen[i] > timen[i]) {
      stop("swtrt_time must be less than or equal to time");
    }
  }
  
  // if the patient switched before pd, set pd time equal to switch time
  for (i=0; i<n; i++) {
    if (pdn[i] == 1 && swtrtn[i] == 1 && swtrt_timen[i] < pd_timen[i]) {
      pd_timen[i] = swtrt_timen[i];
    }
  }
  
  // make sure offset is less than or equal to observed time variables
  for (i=0; i<n; i++) {
    if (pdn[i] == 1 && pd_timen[i] < offset) {
      stop("pd_time must be great than or equal to offset");
    }
    if (swtrtn[i] == 1 && swtrt_timen[i] < offset) {
      stop("swtrt_time must be great than or equal to offset");
    }
  }
  
  // covariates for the Cox model containing treat and base_cov
  std::vector<std::string> covariates(p+1);
  MatrixRM zn_cpp (n,p);//{n,p,std::vector<double>(n*p)};
  covariates[0] = "treated";
  for (j=0; j<p; j++) {
    String zj = base_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base_cov");
    }

    NumericVector u = data[zj];
    covariates[j+1] = std::string(zj.get_cstring());
    //zn(_,j) = u;
    for(int i = 0; i < n; i ++){
      zn_cpp(i,j) = u[i];
    }
  }

  // VERSION WRITTEN WITH R OBJECTS! covariates for the Cox model containing treat and base_cov
  StringVector covariates_r(p+1);
  NumericMatrix zn(n,p);
  covariates_r[0] = "treated";
  for (j=0; j<p; j++) {
    String zj = base_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base_cov");
    }
    NumericVector u = data[zj];
    covariates_r[j+1] = zj;
    zn(_,j) = u;
  }

  // ADDED AS PART OF THE PURE CPP VERSION
  std::vector<std::string> cov_names_cox;
  cov_names_cox.reserve((size_t)1 + (size_t)p);
  cov_names_cox.push_back("treated");
  for (int j = 0; j < p; ++j) {
    cov_names_cox.push_back(Rcpp::as<std::string>(base_cov[j]));
  }
  
  // covariates for the accelerated failure time model for control with pd
  // including stratum and base2_cov
  int q; // number of columns corresponding to the strata effects
  if (strata_main_effect_only) {
    q = sum(d - 1);
  } else {
    q = nstrata - 1;
  }
  
  std::vector<std::string> covariates_aft(q+p2+1);
  MatrixRM zncpp_aft (n,q+p2); //{n,q+p2,std::vector<double>(n * (q+p2))};
  covariates_aft[0] = "swtrt";
  if (strata_main_effect_only) {
    k = 0;
    for (i=0; i<p_stratum; i++) {
      for (j=0; j<d[i]-1; j++) {
        covariates_aft[k+j+1] = "stratum_" + std::to_string(i+1) +
          "_level_" + std::to_string(j+1);
        //zn_aft(_,k+j) = 1.0*(stratan(_,i) == j+1);
        for(int row = 0; row < n; row++){
          zncpp_aft(row,k+j) = (stratan(row,i) == j + 1) ? 1.0 : 0.0;
        }
      }
      k += d[i]-1;
    }
  } else {
    for (j=0; j<nstrata-1; j++) {
      covariates_aft[j+1] = "stratum_" + std::to_string(j+1);
      //zn_aft(_,j) = 1.0*(stratumn == j+1);
      for(int row = 0; row < n; row ++){
        zncpp_aft(row,j) = (stratumn[row] == (j+1)) ? 1.0 : 0.0;
      }
    }
  }

  // R VERSION for debugging
  StringVector covariates_aft_r(q+p2+1);
  NumericMatrix zn_aft(n,q+p2);
  covariates_aft[0] = "swtrt";
  if (strata_main_effect_only) {
    k = 0;
    for (i=0; i<p_stratum; i++) {
      for (j=0; j<d[i]-1; j++) {
        covariates_aft[k+j+1] = "stratum_" + std::to_string(i+1) +
          "_level_" + std::to_string(j+1);
        zn_aft(_,k+j) = 1.0*(stratan(_,i) == j+1);
      }
      k += d[i]-1;
    }
  } else {
    for (j=0; j<nstrata-1; j++) {
      covariates_aft[j+1] = "stratum_" + std::to_string(j+1);
      zn_aft(_,j) = 1.0*(stratumn == j+1);
    }
  }
  
  for (j=0; j<p2; j++) {
    String zj = base2_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base2_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base2_cov");
    }
    NumericVector u = data[zj];
    covariates_aft[q+j+1] = zj;
    zn_aft(_,q+j) = u;
  }
  
  for (j=0; j<p2; j++) {
    String zj = base2_cov[j];
    if (!hasVariable(data, zj)) {
      stop("data must contain the variables in base2_cov");
    }
    if (zj == treat) {
      stop("treat should be excluded from base2_cov");
    }
    NumericVector u = data[zj];
    covariates_aft[q+j+1] = std::string(zj.get_cstring());
    //zn_aft(_,q+j) = u;
    for(int row = 0; row < n; row++){
      zncpp_aft(row,q+j) = u[row];
    }
  }

  std::vector<int> rows(n);
  std::iota(rows.begin(),rows.end(),0);
  
  std::string dist = aft_dist;
  std::for_each(dist.begin(), dist.end(), [](char & c) {
    c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
  });
  
  // error handling for arguments passed
  if ((dist == "log-logistic") || (dist == "llogistic")) {
    dist = "loglogistic";
  } else if  ((dist == "log-normal") || (dist == "lnormal")) {
    dist = "lognormal";
  }
  
  if (!((dist == "exponential") || (dist == "weibull") ||
      (dist == "lognormal") || (dist == "loglogistic"))) {
    stop("aft_dist must be exponential, weibull, lognormal, or loglogistic");
  }
  
  if (alpha <= 0.0 || alpha >= 0.5) {
    stop("alpha must lie between 0 and 0.5");
  }
  
  if (ties != "efron" && ties != "breslow") {
    stop("ties must be efron or breslow");
  }
  
  if (offset < 0.0) {
    stop("offset must be nonnegative");
  }
  
  if (n_boot < 100) {
    stop("n_boot must be greater than or equal to 100");
  }
  
  DataFrame lr = lrtest(data, "", stratum, treat, time, event, 0, 0);
  double logRankPValue = lr["logRankPValue"];
  double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);
  
  // Current goals:
  // - lambda function can be replaced with seperate function call based of f()
  // - fix per thread RNG
  // Important ideas/changes:
  // - remove reliance on names to extract from dataframe because we will explicitly create it in the C++ struct
  k = -1; // indicate the observed data
  auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, 
            treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,
            n, q, p, p2, covariates, covariates_aft, dist, 
            recensor, swtrt_control_only, alpha, zcrit, ties, offset,cov_names_cox
          //r stuff below
          ,covariates_aft_r, covariates_r
          ](
                IntegerVector& idb, IntegerVector& stratumb, 
                NumericVector& timeb, IntegerVector& eventb, 
                IntegerVector& treatb, NumericVector& censor_timeb, 
                IntegerVector& pdb, NumericVector& pd_timeb, 
                IntegerVector& swtrtb, NumericVector& swtrt_timeb, 
                const MatrixRM& zb, const MatrixRM& zb_aft, std::vector<int>& row
              // r stuff below
              )->List {
                  int h, i, j;
                  bool fail = false; // whether any model fails to converge
                  NumericVector init(1, NA_REAL);
                  std::vector<double> init_cpp(1,std::numeric_limits<double>::quiet_NaN()); 
                  
                  List fit1 = List::create(Named("junk") = "dog");
                  liferegOut fit1_cpp; // native cpp version
                  DataFrame data1;
                  List fit_outcome = List::create(Named("junk2") = "dog2");
                  coxfitout fit_outcome_cpp; // native cpp version

                  bool uselifepure = true;
                  bool usephregpure = true;

                  // order data by treat
                  IntegerVector order = seq(0, n-1);
                  std::sort(order.begin(), order.end(), [&](int i, int j) {
                    return (treatb[i] < treatb[j]);
                  });
                  std::vector<int> ordercpp (order.begin(),order.end());
                  
                  // checks if row and order passed are identity (so no need to reorder)
                  bool is_identity = true;
                  for (int i = 0; i < n; ++i) {
                    if (row[i] != i && order[i] != i){ 
                      is_identity = false; 
                      break; 
                    }
                  }
                  if (!is_identity) {
                    Rcpp::Rcout << "identity fail" << std::endl;
                    // reorder vectors AND reorder matrices (see 3B)
                  
                    idb = idb[order];
                    stratumb = stratumb[order];
                    timeb = timeb[order];
                    eventb = eventb[order];
                    treatb = treatb[order];
                    censor_timeb = censor_timeb[order];
                    pdb = pdb[order];
                    pd_timeb = pd_timeb[order];
                    swtrtb = swtrtb[order];
                    swtrt_timeb = swtrt_timeb[order];
                    reorder(row,ordercpp);
                  }
                  if(is_identity){
                    Rcpp::Rcout << "identity check"<< std::endl;
                  }

                  //zb = subset_matrix_by_row(zb, order);
                  //zb_aft = subset_matrix_by_row(zb_aft, order);
                  
                  // time and event adjusted for treatment switching
                  NumericVector t_star = clone(timeb);
                  IntegerVector d_star = clone(eventb);
                  
                  double psi0hat = 0, psi0lower = 0, psi0upper = 0;
                  double psi1hat = 0, psi1lower = 0, psi1upper = 0;
                  
                  // initialize data_aft and fit_aft
                  List data_aft(2), fit_aft(2);
                  if (k == -1) {
                    for (h=0; h<2; h++) {
                      List data_x = List::create(
                        Named("data") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        data_x[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        data_x[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        data_x[treat] = treatwc[1-h];
                      }
                      
                      data_aft[h] = data_x;
                      
                      List fit_x = List::create(
                        Named("fit") = R_NilValue,
                        Named(treat) = R_NilValue
                      );
                      
                      if (TYPEOF(data[treat]) == LGLSXP ||
                          TYPEOF(data[treat]) == INTSXP) {
                        fit_x[treat] = treatwi[1-h];
                      } else if (TYPEOF(data[treat]) == REALSXP) {
                        fit_x[treat] = treatwn[1-h];
                      } else if (TYPEOF(data[treat]) == STRSXP) {
                        fit_x[treat] = treatwc[1-h];
                      }
                      
                      fit_aft[h] = fit_x;
                    }
                  }
                  
                  //DataFrame data_outcome;
                  
                  // treat arms that include patients who switched treatment
                  IntegerVector treats(1);
                  treats[0] = 0;
                  if (!swtrt_control_only) {
                    treats.push_back(1);
                  }
                  
                  int K = static_cast<int>(treats.size());
                  for (h=0; h<K; h++) {
                    // post progression data
                    IntegerVector l = which((treatb == h) & (pdb == 1));
                    std::vector<int> rows_pp(static_cast<size_t>(l.size()));
                    for (int ii = 0; ii < l.size(); ++ii) { 
                      //1 index so we are subtracting -1 for now 
                      rows_pp[ii] = row[l[ii]];  // map logical index -> base row
                    }
                    IntegerVector id2 = idb[l];
                    NumericVector time2 = timeb[l] - pd_timeb[l] + offset;
                    IntegerVector event2 = eventb[l];
                    IntegerVector swtrt2 = swtrtb[l];
                    
                    // prepares all the outputs
                  
                    // 
                    if(uselifepure){
                      // rcpp containers -> stl container -> rcpp List
                      trial_data td;
                      td.pps = to_std<std::vector<double>>(time2);
                      td.event = to_std<std::vector<int>>(event2);
                      td.swtrt = to_std<std::vector<int>>(swtrt2);
                      td.aft = &zb_aft;
                      td.order_pp = &rows_pp;

                      //safety_check_col_aft(td);
                      //Rcpp::Rcout << "[tsesimp_mt] line 853 before lifereg call" << std::endl;
                      fit1_cpp = lifereg_purecpp(
                        td, covariates_aft, dist, 
                        init_cpp, 0, 0, alpha, 50, 1.0e-9);
                      //Rcpp::Rcout << "[tsesimp_mt] line 857 after lifereg call" << std::endl;
                    } 
                    // creates dataframe data1 either as output for k == -1 or to call the regular liferegcpp
                    if(!uselifepure || k == -1){
                      data1 = DataFrame::create(
                        Named("pps") = time2,
                        Named("event") = event2,
                        Named("swtrt") = swtrt2);
                      for (j=0; j<q+p2; j++) {
                        String zj = covariates_aft[j+1];
                        NumericVector u; 
                        for(int row = 0; row < zb_aft.nrows(); row++){ 
                          u.push_back(zb_aft(row,j));
                        }
                        data1.push_back(u[l], zj);
                      }

                      if(!uselifepure){
                        /*
                        Rcpp::Rcout << "using og lifereg" << std::endl;
                      fit1 = liferegcpp(
                      data1, "", "", "pps", "", "event", 
                      covariates_aft_r, "", "", "", dist, init, 
                      0, 0, alpha, 50, 1.0e-9); */
                      }
                    }

                    // we can build this with fit1_cpp later if we really care about this
                    //DataFrame sumstat1 = DataFrame(fit1["sumstat"]);
                    bool fail1 = fit1_cpp.fail; //sumstat1["fail"];
                    if (fail1 == 1) fail = 1;
                    
                    //DataFrame parest1 = DataFrame(fit1["parest"]);
                    std::vector<double> beta1 = fit1_cpp.beta; //["beta"];
                    std::vector<double> sebeta1 = fit1_cpp.sebeta; //["sebeta"];

                    //int isw = k;
                    double psihat = -beta1[1];
                    double psilower = -(beta1[1] + zcrit*sebeta1[1]);
                    double psiupper = -(beta1[1] - zcrit*sebeta1[1]);
                    
                    // calculate counter-factual survival times
                    double a = exp(psihat);
                    double c0 = std::min(1.0, a);
                    for (i=0; i<n; i++) {
                      if (treatb[i] == h) {
                        double b2, u_star, c_star;
                        if (swtrtb[i] == 1) {
                          b2 = pdb[i] == 1 ? pd_timeb[i] : swtrt_timeb[i];
                          b2 = b2 - offset;
                          u_star = b2 + (timeb[i] - b2)*a;
                        } else {
                          u_star = timeb[i];
                        }
                        
                        if (recensor) {
                          c_star = censor_timeb[i]*c0;
                          t_star[i] = std::min(u_star, c_star);
                          d_star[i] = c_star < u_star ? 0 : eventb[i];
                        } else {
                          t_star[i] = u_star;
                          d_star[i] = eventb[i];
                        }
                      }
                    }
                    
                    // update treatment-specific causal parameter estimates
                    if (h == 0) {
                      psi0hat = psihat;
                      psi0lower = psilower;
                      psi0upper = psiupper;
                    } else {
                      psi1hat = psihat;
                      psi1lower = psilower;
                      psi1upper = psiupper;
                    }
                    
                    // update data_aft and fit_aft
                    if (k == -1) {
                      IntegerVector stratum2 = stratumb[l];
                      
                      if (has_stratum) {
                        for (i=0; i<p_stratum; i++) {
                          String s = stratum[i];
                          if (TYPEOF(data[s]) == INTSXP) {
                            IntegerVector stratumwi = u_stratum[s];
                            data1.push_back(stratumwi[stratum2-1], s);
                          } else if (TYPEOF(data[s]) == REALSXP) {
                            NumericVector stratumwn = u_stratum[s];
                            data1.push_back(stratumwn[stratum2-1], s);
                          } else if (TYPEOF(data[s]) == STRSXP) {
                            StringVector stratumwc = u_stratum[s];
                            data1.push_back(stratumwc[stratum2-1], s);
                          }
                        }
                      }
                      
                      if (TYPEOF(data[id]) == INTSXP) {
                        data1.push_front(idwi[id2-1], id);
                      } else if (TYPEOF(data[id]) == REALSXP) {
                        data1.push_front(idwn[id2-1], id);
                      } else if (TYPEOF(data[id]) == STRSXP) {
                        data1.push_front(idwc[id2-1], id);
                      }
                      
                      List data_x = data_aft[h];
                      data_x["data"] = data1;
                      data_aft[h] = data_x;
                      
                      List fit_x = fit_aft[h];
                      fit_x["fit"] = fit1;
                      fit_aft[h] = fit_x;
                    }
                  }
                  
                  DataFrame data_outcome_r;
                  // Cox model for hypothetical treatment effect estimate
                  data_outcome_r = DataFrame::create(
                    Named("uid") = idb,
                    Named("t_star") = t_star,
                    Named("d_star") = d_star,
                    Named("treated") = treatb);
                  
                  data_outcome_r.push_back(stratumb, "ustratum");
                  
                  
                  // creates the data_outcome in r form for debugging/outputting during k == -1
                  for (j=0; j<p; j++) {
                    String zj = covariates[j+1];
                    NumericVector u; //= zb(_,j);
                    for(int row = 0; row < zb.rows; row++){
                      u.push_back(zb(row,j));
                    }
                    data_outcome_r.push_back(u, zj);
                  }
                  if(!usephregpure){/*
                    fit_outcome = phregcpp(
                    data_outcome_r, "", "ustratum", "t_star", "", "d_star",
                    covariates_r, "", "", "", ties, init, 
                    0, 0, 0, 0, 0, alpha, 50, 1.0e-9); */
                  }
                 
                  //
                  if(usephregpure){
                    coxdata data_outcome;
                    //List fit_outcome;
                    data_outcome.n = t_star.size();
                    data_outcome.p = 1 + zb.ncols(); // + 1 bc it gives n of columns and we want to iterate over p 

                    data_outcome.time = to_std<std::vector<double>>(t_star);
                    data_outcome.event = to_std<std::vector<int>>(d_star);
                    data_outcome.stratum = to_std<std::vector<int>>(stratumb);
                    // Pack x row-major: [treated | zb]
                    data_outcome.x.resize((size_t)data_outcome.n * (size_t)data_outcome.p);

                    for (int i = 0; i < data_outcome.n; ++i) {
                      data_outcome.x[(size_t)i * data_outcome.p + 0] = (double)treatb[i];   // treated is column 0
                      for (int j = 0; j < p; ++j) {
                        data_outcome.x[(size_t)i * data_outcome.p + (size_t)(1 + j)] = zb(row[i], j);
                      }
                    }
                    
                    fit_outcome_cpp = phreg_purecpp(data_outcome,cov_names_cox,ties,init_cpp,
                                                0,0,0,0,0,alpha,50,1.0e-9);
                  
                  }

                  //DataFrame sumstat_cox = DataFrame(fit_outcome["sumstat"]);
                  bool fail_cox = fit_outcome_cpp.fail; //["fail"];
                  if (fail_cox == 1) fail = 1;
                  
                  //DataFrame parest = DataFrame(fit_outcome["parest"]);
                  std::vector<double> beta = fit_outcome_cpp.beta; //["beta"];
                  std::vector<double> sebeta = fit_outcome_cpp.sebeta; //parest["sebeta"]; 
                  std::vector<double> pval = fit_outcome_cpp.p; //parest["p"];
                  double hrhat = exp(beta[0]);
                  double hrlower = exp(beta[0] - zcrit*sebeta[0]);
                  double hrupper = exp(beta[0] + zcrit*sebeta[0]);
                  double pvalue = pval[0];
                  
                  List out;
                  if (k == -1) {
                    out = List::create(
                      Named("data_aft") = data_aft,
                      Named("fit_aft") = fit_aft,
                      Named("data_outcome") = data_outcome_r,
                      Named("fit_outcome") = fit_outcome,
                      Named("psihat") = psi0hat,
                      Named("psilower") = psi0lower,
                      Named("psiupper") = psi0upper,
                      Named("psi1hat") = psi1hat,
                      Named("psi1lower") = psi1lower,
                      Named("psi1upper") = psi1upper,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail);
                  } else {
                    out = List::create(
                      Named("psihat") = psi0hat,
                      Named("psilower") = psi0lower,
                      Named("psiupper") = psi0upper,
                      Named("psi1hat") = psi1hat,
                      Named("psi1lower") = psi1lower,
                      Named("psi1upper") = psi1upper,
                      Named("hrhat") = hrhat,
                      Named("hrlower") = hrlower,
                      Named("hrupper") = hrupper,
                      Named("pvalue") = pvalue,
                      Named("fail") = fail);
                  }
                  
                  return out;
                };
  
  List out = f(idn, stratumn, timen, eventn, treatn, censor_timen,
               pdn, pd_timen, swtrtn, swtrt_timen, zn_cpp, zncpp_aft,rows);
  Rcpp::Rcout << "[tsesimp_mt] line 1127, after initial f() call" << std::endl;

  List data_aft = out["data_aft"];
  List fit_aft = out["fit_aft"];
  DataFrame data_outcome = DataFrame(out["data_outcome"]);
  List fit_outcome = out["fit_outcome"];

  IntegerVector uid = data_outcome["uid"];
  if (TYPEOF(data[id]) == INTSXP) {
    data_outcome.push_front(idwi[uid-1], id);
  } else if (TYPEOF(data[id]) == REALSXP) {
    data_outcome.push_front(idwn[uid-1], id);
  } else if (TYPEOF(data[id]) == STRSXP) {
    data_outcome.push_front(idwc[uid-1], id);
  }
  
  IntegerVector treated = data_outcome["treated"];
  if (TYPEOF(data[treat]) == LGLSXP || TYPEOF(data[treat]) == INTSXP) {
    data_outcome.push_back(treatwi[1-treated], treat);
  } else if (TYPEOF(data[treat]) == REALSXP) {
    data_outcome.push_back(treatwn[1-treated], treat);
  } else if (TYPEOF(data[treat]) == STRSXP) {
    data_outcome.push_back(treatwc[1-treated], treat);
  }
  
  if (has_stratum) {
    for (i=0; i<p_stratum; i++) {
      String s = stratum[i];
      if (TYPEOF(data[s]) == INTSXP) {
        IntegerVector stratumwi = u_stratum[s];
        data_outcome.push_back(stratumwi[stratumn-1], s);
      } else if (TYPEOF(data[s]) == REALSXP) {
        NumericVector stratumwn = u_stratum[s];
        data_outcome.push_back(stratumwn[stratumn-1], s);
      } else if (TYPEOF(data[s]) == STRSXP) {
        StringVector stratumwc = u_stratum[s];
        data_outcome.push_back(stratumwc[stratumn-1], s);
      }
    }
  }
  
  double psihat = out["psihat"];
  double psilower = out["psilower"];
  double psiupper = out["psiupper"];
  double psi1hat = out["psi1hat"];
  double psi1lower = out["psi1lower"];
  double psi1upper = out["psi1upper"];
  double hrhat = out["hrhat"];
  double hrlower = out["hrlower"];
  double hrupper = out["hrupper"];
  double pvalue = out["pvalue"];
  bool fail = out["fail"];
  String psi_CI_type = "AFT model";
  
  // construct the confidence interval for HR
  String hr_CI_type;
  //NumericVector hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  //LogicalVector fails(n_boot);
  std::vector<double> hrhats(n_boot), psihats(n_boot), psi1hats(n_boot);
  std::vector<bool> fails(n_boot);

  //debug variables 
  int dbg_n_ok = NA_INTEGER;
  double dbg_sdloghr = NA_REAL;
  double dbg_tcrit = NA_REAL;
  double dbg_M2 = NA_REAL;
  double dbg_min_hr = NA_REAL;
  double dbg_max_hr = NA_REAL;
  double dbg_alpha = alpha;   // capture what C++ thinks alpha is
  if (!boot) { // use Cox model to construct CI for HR if no boot
    hr_CI_type = "Cox model";
  } else { // bootstrap the entire process to construct CI for HR
    //Rcpp::Rcout << "[tsesimp_mt] line 1130, before bootstrapping" << std::endl;
    if (seed != NA_INTEGER) set_seed(seed);
    
    
    //std::vector<int> idb(n), stratumb(n), eventb(n), treatb(n);
    //std::vector<int> pdb(n), swtrtb(n);
    //std::vector<double> timeb(n), censor_timeb(n), pd_timeb(n), swtrt_timeb(n);
    //std::vector<std::vector<double>> zb(n,std::vector<double>(p)), zb_aft(n,std::vector<double>(q+p2));

    //IntegerVector idb(n), stratumb(n), eventb(n), treatb(n);
    //IntegerVector pdb(n), swtrtb(n);
    //NumericVector timeb(n), censor_timeb(n), pd_timeb(n), swtrt_timeb(n);
    //MatrixRM zb (zn_cpp.nrows(),zn_cpp.ncols()), zb_aft (zncpp_aft.nrows(),zncpp_aft.ncols());
    std::vector<int> boot_row(n);
    std::iota(boot_row.begin(),boot_row.end(),0);

    // sort data by treatment group
    IntegerVector idx0 = which(treatn == 0);
    IntegerVector idx1 = which(treatn == 1);
    int n0 = static_cast<int>(idx0.size());
    int n1 = static_cast<int>(idx1.size());
    IntegerVector order(n);
    for (i=0; i<n0; i++) {
      order[i] = idx0[i];
    }
    for (i=0; i<n1; i++){
      order[n0+i] = idx1[i];
    }
    std::vector<int> ordercpp (order.begin(),order.end());
    
    idn = idn[order];
    stratumn = stratumn[order];
    timen = timen[order];
    eventn = eventn[order];
    treatn = treatn[order];
    censor_timen = censor_timen[order];
    pdn = pdn[order];
    pd_timen = pd_timen[order];
    swtrtn = swtrtn[order];
    swtrt_timen = swtrt_timen[order];
    MatrixRM zn_sorted = permute_rows(zn_cpp,ordercpp);
    MatrixRM zn_aft_sorted = permute_rows(zncpp_aft,ordercpp);
    //reorder(boot_row, ordercpp);
    //zn = subset_matrix_by_row(zn, order);
    //zn_aft = subset_matrix_by_row(zn_aft, order);
    std::vector<int> idn_cpp(idn.begin(), idn.end());
    std::vector<int> stratumn_cpp(stratumn.begin(), stratumn.end());
    std::vector<double> timen_cpp(timen.begin(), timen.end());
    std::vector<int> eventn_cpp(eventn.begin(),eventn.end());
    std::vector<int> treatn_cpp(treatn.begin(),treatn.end());
    std::vector<double> censor_timen_cpp(censor_timen.begin(),censor_timen.end());
    std::vector<int> pdn_cpp(pdn.begin(),pdn.end());
    std::vector<double> pd_timen_cpp(pd_timen.begin(),pd_timen.end());
    std::vector<int> swtrtn_cpp(swtrtn.begin(),swtrtn.end());
    std::vector<double> swtrt_timen_cpp(swtrt_timen.begin(),swtrt_timen.end());

    sharedInputs input1 = {n,q,p,p2,
      covariates,covariates_aft,cov_names_cox,
      dist, recensor, swtrt_control_only,
      alpha, zcrit, ties, offset};

    // moved computation of the RNG before worker thread
    std::vector<std::vector<int>> boot_index(n_boot, std::vector<int>(n));
    for (int k = 0; k < n_boot; ++k) {
      for (int i = 0; i < n0; ++i) {
        double u = R::runif(0.0, 1.0);
        int j = static_cast<int>(std::floor(u * n0));
        if (j == n0) j = n0 - 1; // paranoia clamp
        boot_index[k][i] = j;
      }
      for (int i = n0; i < n; ++i) {
        double u = R::runif(0.0, 1.0);
        int j = n0 + static_cast<int>(std::floor(u * n1));
        if (j == n) j = n - 1;
        boot_index[k][i] = j;
      }
    }


    Rcpp::Rcout << "[tsesimp_mt] line 1266, before main pragma loop" << std::endl;

    #pragma omp parallel for schedule(static) default(none) \
      shared(n_boot, n, p, q, p2, boot_index, \
        idn_cpp, stratumn_cpp, timen_cpp, eventn_cpp, treatn_cpp, censor_timen_cpp, pdn_cpp, pd_timen_cpp, swtrtn_cpp, swtrt_timen_cpp, \
        zn_sorted, zn_aft_sorted, fails, hrhats, psihats, psi1hats, input1) 
    for (int k=0; k<n_boot; k++) {
      //thread local buffers
      std::vector<int> idb(n), stratumb(n), eventb(n), treatb(n);
      std::vector<int> pdb(n), swtrtb(n);
      std::vector<double> timeb(n), censor_timeb(n), pd_timeb(n), swtrt_timeb(n);
      MatrixRM zb (n,p), zb_aft (n,q + p2);
      // sample the data with replacement by treatment group
      for (int i=0; i<n; i++) {
        // pre calculated now so we can remove this
        //double u = R::runif(0,1);
        //if (i < n0) {
        //  j = static_cast<int>(std::floor(u*n0));
        //} else {
        //  j = n0 + static_cast<int>(std::floor(u*n1));
        //}
        int j = boot_index[k][i];

        // check scopoe of j here 
        idb[i] = idn_cpp[j];
        stratumb[i] = stratumn_cpp[j];
        timeb[i] = timen_cpp[j];
        eventb[i] = eventn_cpp[j];
        treatb[i] = treatn_cpp[j];
        censor_timeb[i] = censor_timen_cpp[j];
        pdb[i] = pdn_cpp[j];
        pd_timeb[i] = pd_timen_cpp[j];
        swtrtb[i] = swtrtn_cpp[j];
        swtrt_timeb[i] = swtrt_timen_cpp[j];

        const int base_j = j;
        for(int col = 0; col < p; col++){
          zb(i,col) = zn_sorted(base_j,col);
        }
        for(int col = 0; col < q + p2; col++){
          zb_aft(i,col) = zn_aft_sorted(base_j,col);
        }
        //zb(i,_) = zn(j,_);
        //zb_aft(i,_) = zn_aft(j,_);
      }

      std::vector<int> row_rep(n);
      std::iota(row_rep.begin(), row_rep.end(), 0);
      //Rcpp::Rcout << "current boot index:" << k<< std::endl;
      
      replicateInputs input2 = {idb,stratumb, eventb, treatb, 
        pdb, swtrtb, 
        timeb, censor_timeb, pd_timeb, swtrt_timeb,
        zb, zb_aft,row_rep};
      
      //std::vector<int> row_rep(n);
      //std::iota(row_rep.begin(), row_rep.end(), 0);
      //Rcpp::Rcout << "current boot index:" << k<< std::endl;
      
      RepResult output = run_replicate_cpp(k,input1,input2);

      fails[k] = output.fail;
      hrhats[k] = output.hrhat;
      psihats[k] = output.psi0hat;
      psi1hats[k] = output.psi1hat;
    }
    Rcpp::Rcout << "[tsesimp_mt] line 1345, after bootstrapping loop" << std::endl;
    
    // obtain bootstrap confidence interval for HR
    //const double loghr = std::log(hrhat);

    // Build canonical ok: non-fail AND valid HR

    // double checks its not improperly built
    /*for(int i = 0; i < n_boot; ++i){
      if(!fails[i]){
        double h = hrhats[i];
        if (!std::isfinite(h) || h <= 0.0){
          fails[i] = true;
        }
      }
    }
    LogicalVector ok(n_boot);
    for (int i = 0; i < n_boot; ++i) {
      ok[i] = !fails[i];
    }
    const int n_ok = Rcpp::sum(ok);

    NumericVector loghrs = log(hrhats[ok]);
    //int idx_lh = 0;
    //for (int i = 0; i < n_boot; ++i) {
      //  if (!ok[i]) continue;
      //  loghrs[idx_lh++] = std::log(hrhats[i]);
    //}

    // sample standard deviation with n_ok - 1 denominator
    double sdloghr = sd(loghrs);
    double tcrit   = R::qt(1.0 - alpha / 2.0, n_ok - 1, 1, 0);
    hrlower = std::exp(log(hrhat) - tcrit * sdloghr);
    hrupper = std::exp(log(hrhat) + tcrit * sdloghr);

    Rcpp::Rcout << "boot debug: n_boot" << n_boot << std::endl;
    Rcpp::Rcout << "n_ok =" << n_ok << std::endl; 
    Rcpp::Rcout << "sum_fails=" << std::accumulate(fails.begin(),fails.end(),0)
    << " sdloghr=" << sdloghr
    << " tcrit=" << tcrit
    << " hrhat=" << hrhat
    << " hr_CI_type=" << hr_CI_type.get_cstring()
    << std::endl;*/
    // obtain bootstrap confidence interval for HR
    const double loghr = std::log(hrhat);

    // Build canonical ok: non-fail AND valid HR
    LogicalVector ok(n_boot);
    for (int i = 0; i < n_boot; ++i) {
      const double h = hrhats[i];
      ok[i] = (!fails[i]) && std::isfinite(h) && (h > 0.0);
    }
    const int n_ok = Rcpp::sum(ok);

    double sdloghr = NA_REAL;
    double tcrit = NA_REAL;
    int count = 0;
    if (n_ok >= 2) {
    double mean = 0.0, M2 = 0.0;

    for (int i = 0; i < n_boot; ++i) {
      if (!ok[i]) continue;
      const double x = std::log(hrhats[i]);
      ++count;
      const double delta = x - mean;
      mean += delta / count;
      const double delta2 = x - mean;
      M2 += delta * delta2;
    }

    sdloghr = std::sqrt(M2 / (count - 1));
    tcrit = R::qt(1.0 - alpha / 2.0, count - 1, 1, 0);

    dbg_n_ok = count;
    dbg_sdloghr = sdloghr;
    dbg_tcrit = tcrit;
    dbg_M2 = M2;

    hrlower = std::exp(loghr - tcrit * sdloghr);
    hrupper = std::exp(loghr + tcrit * sdloghr);
    hr_CI_type = "bootstrap";

    const double tstat = (sdloghr > 0.0) ? (loghr / sdloghr) : 0.0;
    pvalue = (sdloghr > 0.0)
      ? 2.0 * (1.0 - R::pt(std::fabs(tstat), count - 1, 1, 0))
      : NA_REAL;

    } else {
      hrlower = hrhat;
      hrupper = hrhat;
      hr_CI_type = "bootstrap (n_ok < 2)";
    }

    double min_h = std::numeric_limits<double>::infinity();
    double max_h = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_boot; ++i) {
      if (!ok[i]) continue;
      min_h = std::min(min_h, hrhats[i]);
      max_h = std::max(max_h, hrhats[i]);
    }
    dbg_min_hr = std::isfinite(min_h) ? min_h : NA_REAL;
    dbg_max_hr = std::isfinite(max_h) ? max_h : NA_REAL;

    Rcpp::Rcout << "boot debug: n_boot" << n_boot << std::endl;
    Rcpp::Rcout << "n_ok =" << n_ok << std::endl; 
    Rcpp::Rcout << "sum_fails=" << std::accumulate(fails.begin(),fails.end(),0)
    << " sdloghr=" << sdloghr
    << " tcrit=" << tcrit
    << " hrhat=" << hrhat
    << " hr_CI_type=" << hr_CI_type.get_cstring()
    << std::endl;

    Rcpp::NumericVector psihats1(count);
    Rcpp::NumericVector psi1hats1(count);

    int pos = 0;
    for (int i = 0; i < n_boot; ++i) {
      if (!ok[i]) continue;
      psihats1[pos]  = psihats[i];    // std::vector<double>
      psi1hats1[pos] = psi1hats[i];   // std::vector<double>
      ++pos;
    }

    //NumericVector psihats1 = psihats[ok];
    double sdpsi = sd(psihats1);
    psilower = psihat - tcrit*sdpsi;
    psiupper = psihat + tcrit*sdpsi;
    psi_CI_type = "bootstrap";
    
    //NumericVector psi1hats1 = psi1hats[ok];
    double sdpsi1 = sd(psi1hats1);
    psi1lower = psi1hat - tcrit*sdpsi1;
    psi1upper = psi1hat + tcrit*sdpsi1;
  }
  
  List settings = List::create(
    Named("aft_dist") = aft_dist,
    Named("strata_main_effect_only") = strata_main_effect_only,
    Named("recensor") = recensor,
    Named("admin_recensor_only") = admin_recensor_only,
    Named("swtrt_control_only") = swtrt_control_only,
    Named("alpha") = alpha,
    Named("ties") = ties,
    Named("offset") = offset,
    Named("boot") = boot,
    Named("n_boot") = n_boot,
    Named("seed") = seed);
  
  List result = List::create(
    Named("psi") = psihat,
    Named("psi_CI") = NumericVector::create(psilower, psiupper),
    Named("psi_CI_type") = psi_CI_type,
    Named("logrank_pvalue") = 2*std::min(logRankPValue, 1-logRankPValue),
    Named("cox_pvalue") = pvalue,
    Named("hr") = hrhat,
    Named("hr_CI") = NumericVector::create(hrlower, hrupper),
    Named("hr_CI_type") = hr_CI_type,
    Named("data_aft") = data_aft,
    Named("fit_aft") = fit_aft,
    Named("data_outcome") = data_outcome,
    Named("fit_outcome") = fit_outcome,
    Named("fail") = fail,
    Named("settings") = settings,
  
    //debug result stuff not needed
    Named("debug_alpha") = dbg_alpha,
    Named("debug_n_ok") = dbg_n_ok,
    Named("debug_sdloghr") = dbg_sdloghr,
    Named("debug_tcrit") = dbg_tcrit,
    Named("debug_M2") = dbg_M2,
    Named("debug_min_hr_boot") = dbg_min_hr,
    Named("debug_max_hr_boot") = dbg_max_hr);
    result.push_back("DEBUG_TAG_2026_01_01", "__debug_tag");

  if (!swtrt_control_only) {
    result.push_back(psi1hat, "psi_trt");
    NumericVector psi1_CI = NumericVector::create(psi1lower, psi1upper);
    result.push_back(psi1_CI, "psi_trt_CI");
  }
  
  if (boot) {
    result.push_back(fails, "fail_boots");
    result.push_back(hrhats, "hr_boots");
    result.push_back(psihats, "psi_boots");
    if (!swtrt_control_only) {
      result.push_back(psi1hats, "psi_trt_boots");
    }
  }
  
  return result;
}

