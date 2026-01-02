[1mdiff --git a/src/tsesimp_mt.cpp b/src/tsesimp_mt.cpp[m
[1mindex a13ab94..05d5be6 100644[m
[1m--- a/src/tsesimp_mt.cpp[m
[1m+++ b/src/tsesimp_mt.cpp[m
[36m@@ -9,6 +9,7 @@[m [musing namespace Rcpp;[m
 //add custom includes[m
 #include "utilities_pure.h"[m
 #include "survival_analysis_pure.h"[m
[32m+[m[32m#include <R_ext/Print.h>[m
 [m
 // [[Rcpp::export]][m
 List tsesimpcpp_mt(const DataFrame data,[m
[36m@@ -37,12 +38,13 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                 const int seed = NA_INTEGER) {[m
   [m
   // checking it gets to the C++ code[m
[31m-  Rcpp::Rcout << "[tsesimpcpp_mt] time arg = '" << time << "'\n";[m
[31m-  Rcpp::CharacterVector nms = data.names();[m
[31m-  Rcpp::Rcout << "[tsesimpcpp_mt] names(data): ";[m
[31m-  for (auto s : nms) Rcpp::Rcout << Rcpp::as<std::string>(s) << " ";[m
[31m-  Rcpp::Rcout << "\n";[m
[32m+[m[32m  //Rcpp::Rcout << "[tsesimpcpp_mt] time arg = '" << time << "'\n";[m
[32m+[m[32m  //Rcpp::CharacterVector nms = data.names();[m
[32m+[m[32m  //Rcpp::Rcout << "[tsesimpcpp_mt] names(data): ";[m
[32m+[m[32m  //for (auto s : nms) Rcpp::Rcout << Rcpp::as<std::string>(s) << " ";[m
[32m+[m[32m  //Rcpp::Rcout << "\n";[m
 [m
[32m+[m[32m  Rcpp::RNGScope scope;[m
   int i, j, k, n = data.nrow();[m
   [m
   int p = static_cast<int>(base_cov.size());[m
[36m@@ -411,6 +413,11 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
   double zcrit = R::qnorm(1-alpha/2, 0, 1, 1, 0);[m
   [m
   k = -1; // indicate the observed data[m
[32m+[m[32m  // Current goals:[m
[32m+[m[32m  // - entire lambda function is to be rewritten in pure C++ to allow multi-threadding[m
[32m+[m[32m  // - fix per thread RNG[m
[32m+[m[32m  // Important ideas/changes:[m
[32m+[m[32m  // - remove reliance on names to extract from dataframe because we will explicitly create it in the C++ struct[m
   auto f = [&k, data, has_stratum, stratum, p_stratum, u_stratum, [m
             treat, treatwi, treatwn, treatwc, id, idwi, idwn, idwc,[m
             n, q, p, p2, covariates, covariates_aft, dist, [m
[36m@@ -489,7 +496,7 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                     }[m
                   }[m
                   [m
[31m-                  DataFrame data_outcome;[m
[32m+[m[32m                  //DataFrame data_outcome;[m
                   [m
                   // treat arms that include patients who switched treatment[m
                   IntegerVector treats(1);[m
[36m@@ -507,33 +514,13 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                     IntegerVector event2 = eventb[l];[m
                     IntegerVector swtrt2 = swtrtb[l];[m
                     [m
[31m-                    //DataFrame data1 = DataFrame::create([m
[31m-                      //Named("pps") = time2,[m
[31m-                      //Named("event") = event2,[m
[31m-                      //Named("swtrt") = swtrt2);[m
[31m-                    [m
[31m-                    /*DataFrame data1 = DataFrame::create([m
[31m-                      Named("pps") = time2,[m
[31m-                      Named("event") = event2,[m
[31m-                      Named("swtrt") = swtrt2);[m
[31m-                    [m
[31m-                    for (j=0; j<q+p2; j++) {[m
[31m-                      String zj = covariates_aft[j+1];[m
[31m-                      NumericVector u = zb_aft(_,j);[m
[31m-                      data1.push_back(u[l], zj);[m
[31m-                    }[m
[31m-                    */[m
[31m-                    //initialize this to get behavior ready for some of the code[m
[31m-                    //[m
                     List fit1;[m
                     DataFrame data1;[m
[31m-[m
[31m-                    if(k>=0){[m
[32m+[m[32m                    if(true){[m
[32m+[m[32m                    //if(k>=0){[m
                       //modify this from our regular standard function into[m
                       // rcpp containers -> stl container -> rcpp container[m
                       // a lot of this is a temporary band aid fix until finalization of the pipeline[m
[31m-                      Rcpp::Rcout << "[tsesimp_mt] enter\n";[m
[31m-                      Rcpp::Rcout << "[tsesimp_mt] building td...\n";[m
                       std::vector<std::string> covariates_aft_std = to_std<std::vector<std::string>>(covariates_aft);[m
                       trial_data td;[m
                       td.pps = to_std<std::vector<double>>(time2);[m
[36m@@ -556,31 +543,21 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                         std::vector<double> col; [m
                         col.reserve(l.size());[m
                         for(int ii = 0; ii < l.size(); ++ii) {[m
[31m-                          int i = l[ii] - 1; // -1 for 0-based indexing from R's 1-based[m
[32m+[m[32m                          int i = l[ii]; // -1 for 0-based indexing from R's 1-based but not needed[m
                           col.push_back(zb_aft(i,j));[m
                         }[m
                         td.aft_names.push_back(name);[m
                         td.aft.push_back(std::move(col));[m
                       }[m
[31m-                      std::vector<double> init_cpp; // leave empty on purpose[m
[32m+[m[32m                      std::vector<double> init_cpp;[m[41m [m
 [m
[31m-                      //ensures we have all the correct data and correct sizes[m
                       safety_check_col_aft(td);[m
[31m-                      // after filling td:[m
[31m-                      Rcpp::Rcout << "[tsesimp_mt] td sizes pps=" << td.pps.size()[m
[31m-                      << " event=" << td.event.size()[m
[31m-                      << " swtrt=" << td.swtrt.size()[m
[31m-                      << " aft_names=" << td.aft_names.size()[m
[31m-                      << " aft_cols=" << td.aft.size() << std::endl;[m
[31m-[m
[31m-                      // right before calling lifereg_purecpp:[m
[31m-                      Rcpp::Rcout << "[tsesimp_mt] calling lifereg_purecpp()\n";[m
                       fit1 = lifereg_purecpp([m
                         td, {}, {}, "pps", "", "event", [m
                         covariates_aft_std, "", "", "", dist, init_cpp, 0, 0, alpha, [m
                         50, 1.0e-9);[m
[31m-                      Rcpp::Rcout << "[tsesimp_mt] lifereg_purecpp() returned\n";[m
[31m-                    } else {[m
[32m+[m[32m                    }[m[41m [m
[32m+[m[32m                    if(k == -1){[m
                       data1 = DataFrame::create([m
                         Named("pps") = time2,[m
                         Named("event") = event2,[m
[36m@@ -590,12 +567,8 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                         NumericVector u = zb_aft(_,j);[m
                         data1.push_back(u[l], zj);[m
                       }[m
[31m-                      fit1 = liferegcpp([m
[31m-                        data1, "", "", "pps", "", "event", [m
[31m-                        covariates_aft, "", "", "", dist, init, [m
[31m-                        0, 0, alpha, 50, 1.0e-9);[m
                     }[m
[31m-                    [m
[32m+[m
                     DataFrame sumstat1 = DataFrame(fit1["sumstat"]);[m
                     bool fail1 = sumstat1["fail"];[m
                     if (fail1 == 1) fail = 1;[m
[36m@@ -603,24 +576,11 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                     DataFrame parest1 = DataFrame(fit1["parest"]);[m
                     NumericVector beta1 = parest1["beta"];[m
                     NumericVector sebeta1 = parest1["sebeta"];[m
[31m-                    [m
 [m
[31m-                    //BUG FIX LOGIC:[m
[31m-                    CharacterVector names1 = parest1.names();[m
[31m-                    int isw = -1;[m
[31m-                    for(int ii = 0; ii < names1.size(); ++ii) {[m
[31m-                      if(names1[ii] == "swtrt") {[m
[31m-                        isw = ii;[m
[31m-                        break;[m
[31m-                      }[m
[31m-                    }[m
[31m-[m
[31m-                    if(isw < 0) {[m
[31m-                      stop("swtrt not found in parest1 names");[m
[31m-                    }[m
[31m-                    double psihat = -beta1[1];[m
[31m-                    double psilower = -(beta1[1] + zcrit*sebeta1[1]);[m
[31m-                    double psiupper = -(beta1[1] - zcrit*sebeta1[1]);[m
[32m+[m[32m                    int isw = k;[m
[32m+[m[32m                    double psihat = -beta1[isw];[m
[32m+[m[32m                    double psilower = -(beta1[isw] + zcrit*sebeta1[isw]);[m
[32m+[m[32m                    double psiupper = -(beta1[isw] - zcrit*sebeta1[isw]);[m
                     [m
                     // calculate counter-factual survival times[m
                     double a = exp(psihat);[m
[36m@@ -696,6 +656,7 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                     }[m
                   }[m
                   [m
[32m+[m[32m                  /*[m
                   // Cox model for hypothetical treatment effect estimate[m
                   data_outcome = DataFrame::create([m
                     Named("uid") = idb,[m
[36m@@ -715,7 +676,23 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
                     data_outcome, "", "ustratum", "t_star", "", "d_star",[m
                     covariates, "", "", "", ties, init, [m
                     0, 0, 0, 0, 0, alpha, 50, 1.0e-9);[m
[31m-                  [m
[32m+[m[32m                  */[m
[32m+[m[41m                 [m
[32m+[m[32m                  //[m
[32m+[m[32m                  coxdata data_outcome;[m
[32m+[m[32m                  List fit_outcome;[m
[32m+[m[32m                  data_outcome.n = t_star.size();[m
[32m+[m[32m                  data_outcome.p = 1 + zb.ncol(); // + 1 bc it gives n of columns and we want to iterate over p[m[41m [m
[32m+[m[32m                  data_outcome.time = to_std<std::vector<double>>(t_star);[m
[32m+[m[32m                  data_outcome.event = to_std<std::vector<int>>(d_star);[m
[32m+[m[32m                  data_outcome.stratum = to_std<std::vector<int>>(stratumb);[m
[32m+[m
[32m+[m[32m                  if(true){[m
[32m+[m
[32m+[m[32m                    fit_outcome = phreg_purecpp(data_outcome);[m
[32m+[m[32m                  }[m
[32m+[m
[32m+[m
                   DataFrame sumstat_cox = DataFrame(fit_outcome["sumstat"]);[m
                   bool fail_cox = sumstat_cox["fail"];[m
                   if (fail_cox == 1) fail = 1;[m
[36m@@ -895,7 +872,114 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
     }[m
     [m
     // obtain bootstrap confidence interval for HR[m
[31m-    double loghr = log(hrhat);[m
[32m+[m[32m    double loghr = std::log(hrhat);[m
[32m+[m[32m    double sdloghr = -1, tcrit = -1;[m
[32m+[m[32m    LogicalVector ok = 1 - fails;[m
[32m+[m[41m    [m
[32m+[m[32m    /* Beta modified 1[m
[32m+[m[32m    int n_ok = 0;[m[41m [m
[32m+[m[32m    double sum_loghr = 0.0, sumsq_log_hr = 0.0;[m
[32m+[m
[32m+[m[32m    for(int z = 0; z < n_boot; z++){[m
[32m+[m[32m      if(!fails[z]){[m
[32m+[m[32m        double z1 = hrhats[z];[m
[32m+[m[32m        if(std::isfinite(z1) && z1 > 0.0){[m
[32m+[m[32m          double lh = std::log(z1);[m
[32m+[m[32m          sum_loghr += lh;[m
[32m+[m[32m          sumsq_log_hr += lh * lh;[m
[32m+[m[32m          n_ok++;[m
[32m+[m[32m        }[m
[32m+[m[32m      }[m
[32m+[m[32m    }[m
[32m+[m
[32m+[m[32m    double sdloghr = -1, tcrit = -1;[m
[32m+[m
[32m+[m[32m    if(n_ok >= 2){[m
[32m+[m[32m      double mean = sum_loghr / n_ok;[m
[32m+[m[32m      double var = (sumsq_log_hr - n_ok * mean * mean) / (n_ok -1);[m
[32m+[m[32m      sdloghr = std::sqrt(std::max(0.0,var));[m
[32m+[m[32m      tcrit = R::qt( 1- alpha /2, n_ok -1, 1, 0);[m
[32m+[m
[32m+[m[32m      hrlower = std::exp(loghr - tcrit * sdloghr);[m
[32m+[m[32m      hrupper = std::exp(loghr + tcrit * sdloghr);[m
[32m+[m[32m      hr_CI_type = "bootstrap";[m
[32m+[m[32m    } else{[m
[32m+[m[32m      hrlower = hrhat;[m
[32m+[m[32m      hrupper = hrhat;[m
[32m+[m[32m      hr_CI_type = "bootstrap (n_ok) < 2";[m
[32m+[m[32m    }[m
[32m+[m[32m    */[m
[32m+[m
[32m+[m[32m    int n_ok = 0;[m
[32m+[m[32m    double mean = 0.0, M2 = 0.0;[m
[32m+[m[32m    for(int i = 0; i < n_boot; i++){[m
[32m+[m[32m      if(!fails[i]){[m
[32m+[m[32m        double h = hrhats[i];[m
[32m+[m[32m        /*TEST PART 3 PART 1[m
[32m+[m[32m        n_ok++;[m
[32m+[m[32m        double delta = h - mean;[m
[32m+[m[32m        mean += delta / n_ok;[m[41m [m
[32m+[m[32m        M2 += delta / (h - mean);[m
[32m+[m[32m        */[m
[32m+[m[41m        [m
[32m+[m[32m        if(std::isfinite(h) && h > 0.0){[m
[32m+[m[32m          double x = std::log(h);[m
[32m+[m[32m          n_ok++;[m
[32m+[m[32m          double delta = x - mean;[m
[32m+[m[32m          mean += delta/n_ok;[m
[32m+[m[32m          double delta2 = x - mean;[m
[32m+[m[32m          M2 += delta * delta2;[m
[32m+[m[32m        }[m
[32m+[m[32m      }[m
[32m+[m[32m    }[m
[32m+[m[41m    [m
[32m+[m[32m    //double loghr = std::log(hrhat);[m
[32m+[m[32m    double tstat;[m
[32m+[m
[32m+[m[32m    /*TEST PART 3 PART 2[m
[32m+[m[32m    sdloghr = std::sqrt(M2 / (n_ok - 1));[m
[32m+[m[32m    tcrit = R::qt(1.0 - alpha / 2.0, n_ok - 1, 1, 0);[m
[32m+[m[32m    hrlower = std::exp(loghr - tcrit * sdloghr);[m
[32m+[m[32m    hrupper = std::exp(loghr + tcrit * sdloghr);[m
[32m+[m[32m    hr_CI_type = "bootstrap";[m
[32m+[m[32m    pvalue = 2.0 * (1.0 - R::pt(std::fabs(loghr/sdloghr), n_ok - 1, 1, 0));[m
[32m+[m[32m    */[m
[32m+[m
[32m+[m
[32m+[m
[32m+[m[32m    if(n_ok >=2){[m
[32m+[m[32m      sdloghr = std::sqrt(M2/ (n_ok -1));[m
[32m+[m[32m      tcrit = R::qt(1.0 - alpha/ 2.0, n_ok - 1, 1 , 0);[m
[32m+[m[32m      hrlower = std::exp(loghr - tcrit * sdloghr);[m
[32m+[m[32m      hrupper = std:: exp(loghr + tcrit * sdloghr);[m
[32m+[m[32m      hr_CI_type = "bootstrap";[m
[32m+[m[32m      tstat = (sdloghr > 0.0) ? (loghr / sdloghr) : 0.0;[m
[32m+[m[32m      pvalue = (sdloghr > 0.0)[m
[32m+[m[32m              ? 2.0 * (1.0 - R::pt(std::fabs(tstat), n_ok - 1, 1, 0))[m
[32m+[m[32m              : NA_REAL;[m
[32m+[m[32m    } else{[m
[32m+[m[32m      hrlower = hrhat;[m
[32m+[m[32m      hrupper = hrhat;[m
[32m+[m[32m      hr_CI_type = "bootstrap (n_ok < 2)";[m
[32m+[m[32m    }[m
[32m+[m[41m  [m
[32m+[m[32m    /*TEST[m
[32m+[m[32m    NumericVector psihats1 = psihats[ok];[m
[32m+[m[32m    int m_ok = 0;[m
[32m+[m[32m    double m_mean = 0.0, m_M2 = 0.0;[m
[32m+[m[32m    for(int i = 0; i < n_boot; ++i){[m
[32m+[m[32m      if(!fails[i]){[m
[32m+[m[32m        double x = psihats[i];[m
[32m+[m[32m        if(std::isfinite(x)){[m
[32m+[m[32m          m_ok++;[m
[32m+[m
[32m+[m[32m        }[m
[32m+[m[32m      }[m
[32m+[m[32m    }[m
[32m+[m[32m    */[m
[32m+[m
[32m+[m
[32m+[m[32m    /* ORIGINAL CODE[m
     LogicalVector ok = !fails;[m
     int n_ok = sum(ok);[m
     NumericVector loghrs = log(hrhats[ok]);[m
[36m@@ -903,10 +987,23 @@[m [mList tsesimpcpp_mt(const DataFrame data,[m
     double tcrit = R::qt(1-alpha/2, n_ok-1, 1, 0);[m
     hrlower = exp(loghr - tcrit*sdloghr);[m
     hrupper = exp(loghr + tcrit*sdloghr);[m
[32m+[m[41m    [m
[32m+[m[41m    [m
[32m+[m[32m    //loghr -> log(hrhat) -> from initial -1 building[m[41m [m
[32m+[m[32m    //tcrit[m[41m [m
[32m+[m[32m    //sdloghr[m
     hr_CI_type = "bootstrap";[m
[32m+[m[41m    [m
[32m+[m
[32m+[m
     pvalue = 2*(1 - R::pt(fabs(loghr/sdloghr), n_ok-1, 1, 0));[m
[32m+[m[32m    */[m
[32m+[m[32m    //Rprintf("[boot] n_ok=%d\n", n_ok);[m
[32m+[m[32m    //Rprintf("[boot] sdloghr=%.6f, tcrit=%.6f\n", sdloghr, tcrit);[m
[32m+[m[32m    //Rprintf("[boot] CI=(%.6f, %.6f)\n", hrlower, hrupper);[m
     [m
     // obtain bootstrap confidence interval for psi[m
[32m+[m[32m    //LogicalVector ok = !fails;[m
     NumericVector psihats1 = psihats[ok];[m
     double sdpsi = sd(psihats1);[m
     psilower = psihat - tcrit*sdpsi;[m
