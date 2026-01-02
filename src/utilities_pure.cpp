#include <Rcpp.h>
//#include <R_ext/Applic.h>
#include "utilities_pure.h"

#include <vector>
#include <string>
#include <cmath>
//these bottom three are introduced to replace the r method calls that use more complicated math functions
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/chi_squared.hpp>

#include <limits>

#include <Rmath.h>

namespace bm = boost::math;
/*
A majority of these functions are rewritten from Rcpp to use std::vector or stl structs instead
to avoid Rcpp API calls and make them usable in pure C++.

The main idea is to keep the pipeline of data processing similar to the original Rcpp code
*/

//TO DO:
//double check code to have size_t i = 0, as loop because we can't compare a integer with a size_t for example
//vector length is unsigned integers so large vectors can cause issues with signed integers
//also will probably rename xxxxx_cpp to xxxxx_ in the future, but I want to avoid re writing the same name for functions


// helpers specifically for rcpp to pure cpp conversion
namespace details {
  constexpr double INV_SQRT_2PI = 0.39894228040143267793994605993438; // 1/sqrt(2*pi)
  constexpr double LOG_SQRT_2PI = 0.91893853320467274178032973640562; // log(sqrt(2*pi))
  inline double NaN() { return std::numeric_limits<double>::quiet_NaN(); }
}




typedef std::vector<std::vector<double>> matrix;


std::vector<int> seq_cpp(int start, int end) {
  if(end < start) return {};
  int len = end - start + 1;
  std::vector<int> v;
  v.reserve(len);
  for (int i = 0; i < len; ++i) {
    v.push_back(start + i);
  }
  return v;
}
/*
template<typename T>
void reorder(std::vector<T>& vec, const std::vector<int>& order) {
  std::vector<T> old = vec;
  for (size_t i = 0; i < order.size(); ++i) {
    vec[i] = old[order[i]];
  }
}
*/

matrix subset_matrix_by_row_cpp(const matrix& mat, const std::vector<int>& order) {
  if(order.empty()){
    if(mat.empty()) return {};
    return std::vector<std::vector<double>>{};
  }
  if(mat.empty()){
    return {};
  }

  int n = static_cast<int>(order.size());
  int p = static_cast<int>(mat[0].size());
  
  #ifndef NDEBUG
  for (int r = 0; r < n; ++r) {
    if (static_cast<int>(mat[r].size()) != p) {
      throw std::runtime_error("subset_matrix_by_row_cpp: ragged matrix");
    }
  }
  for (int idx : order) {
    if (idx < 0 || idx >= n) {
      throw std::out_of_range("subset_matrix_by_row_cpp: row index OOB");
    }
  }
  #endif

  matrix result(n, std::vector<double>(p));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < p; ++j) {
      result[i][j] = mat[order[i]][j];
    }
  }
  return result;
}

// sum_bool_cpp for boolean vectors
int sum_bool_cpp(const std::vector<bool>& vec) {
  int total = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    if(vec[i]) {
      total++;
    }
  }
  return total;
}

double mean_cpp(const std::vector<double>& vec) {
    if(vec.empty()) {
        throw std::invalid_argument("Mean is not defined for empty vector");
    }
    double total = 0.0;
    for (const auto& val : vec) {
        total += val;
    }
    return total / static_cast<double>(vec.size());
}

double sum_cpp(const std::vector<double>& vec) {
    double total = 0.0;
    for (const auto& val : vec) {
        total += val;
    }
    return total;
}

double sd_cpp(const std::vector<double>& vec) {
    if(vec.empty()) {
        throw std::invalid_argument("Standard deviation is not defined for empty vector");
    }
    if (vec.size() < 2) {
        throw std::invalid_argument("Standard deviation is not defined for size < 2");
    }
    if(vec.size() == 2) {
        return std::abs(vec[0] - vec[1]) / std::sqrt(2.0); // special case for two elements
    }
    double m = mean_cpp(vec);
    double sqsum = 0.0;
    for (const auto& val : vec) {
        double d = val - m;
        sqsum += d * d;
    }
    return std::sqrt(sqsum / static_cast<double>(vec.size() - 1));
}

//bottom few functions are rewritten from R::(FUNCTION_NAME) to use vectors instead but R supports it's own math library so we add Boost math equivalents to help us with certain tougher functions to recreate like
//added type swapping to boost accuracy in the lower decimal places for extreme values
double dnorm_cpp(double x, bool logp){
  return R::dnorm(x, 0.0, 1.0,logp ? 1:0);
  /*double u = -0.5 * x * x;
  if (logp) {
      return u - details::LOG_SQRT_2PI;
  } else {
      return std::exp(u) * details::INV_SQRT_2PI;
  }
      */
}

double pnorm_cpp(double x, bool lower_tail, bool log_p){
  return R::pnorm(x,0.0,1.0,lower_tail ? 1:0, log_p ? 1:0);
  
  /*
  using ld = long double;
  if (std::isnan(x)) return details::NaN();

  static const bm::normal_distribution<ld> dist(0.0L,1.0L);

  ld p;
  if(lower_tail){
    p = bm::cdf(dist, (ld)x);
  } else {
    p = bm::cdf(bm::complement(dist, (ld)x));
  }

  if(log_p){
    return (double)std::log(p);
  }
  return (double)p;
  */
}

double qnorm_cpp(double p, bool lower_tail, bool log_p){
  return R::qnorm(p, 0.0, 1.0, lower_tail ? 1 : 0, log_p ? 1 : 0);
  /*
  using ld = long double;
  if(std::isnan(p)) return details::NaN();
  
  ld pp;
  static const bm::normal_distribution<ld> dist(0.0L,1.0L);

  if (log_p) {
    pp = std::exp((ld)p);
  } 
  else{
    pp = (ld)p;
  }

  if ( !lower_tail ) {
    pp = 1.0L - pp; // upper tail
  }
  
  if(pp <= 0.0L) return -std::numeric_limits<double>::infinity();
  if(pp >= 1.0L) return std::numeric_limits<double>::infinity();

  ld quantile = bm::quantile(dist, pp);
  return (double)quantile;
  */
}

double dlogis_cpp(double x, bool log_p){
  return R::dlogis(x, 0.0, 1.0, log_p ? 1 :0);
  /*
  if(log_p){
    if(x >= 0.0){
      return -x - 2.0 * std::log1p(std::exp(-x));
    }
    else{
      return x - 2.0 * std::log1p(std::exp(x));
    }
  }
  else{
    if(x >= 0.0){
      double e = std::exp(-x);
      double denom = (1 + e);
      return e / (denom*denom);
    }
    else{
      double e = std::exp(x);
      double denom = (1 + e);
      return e / (denom*denom);
    }
  }
    */
}

double plogis_cpp(double x, bool lower_tail, bool log_p){
  return R::plogis(x, 0.0, 1.0, lower_tail ? 1:0, log_p ? 1: 0);
  /*
  if(log_p){
    if(lower_tail){
      if(x >= 0.0){
        return -std::log1p(std::exp(-x));
      }
      else{
        return x - std::log1p(std::exp(x));
      }
    }
    else{
      if(x <= 0.0){
        return -std::log1p(std::exp(x));
      }
      else{
        return -x - std::log1p(std::exp(-x));
      }
    }
  }
  else{
    if(lower_tail){
      if(x >= 0.0){
        double e = std::exp(-x);
        return 1.0 / (1.0 + e);
      }
      else{
        double e = std::exp(x);
        return e / (1.0 + e);
      }
    }
    else{
      if( x >= 0.0){
        double e = std::exp(-x);
        return e / (1.0 + e);
      }
      else{
        double e = std::exp(x);
        return 1.0 / (1.0 + e);
      }
    }
  }
    */
}

double pchisq_cpp(double x, double df, bool lower_tail, bool log_p){
  return R::pchisq(x,df,lower_tail ? 1:0, log_p ? 1:0);
  /*
  using ld = long double;
  if(!std::isfinite(df) || df <= 0){
    //Rcpp::Rcout << "[pchisq]_cpp invalid df" << df << "x= " << x << std::endl;
    return details::NaN();
  }
  if(std::isnan(x)) return details::NaN();

  if(x < 0.0){
    double prob = lower_tail ? 0.0 : 1.0;
    return log_p ? std::log(prob) : prob;
  }

  if(!std::isfinite(x)){
    double prob = lower_tail ? 1.0 : 0.0;
    return log_p ? std::log(prob) : prob;
  }

  bm::chi_squared_distribution<ld> dist((ld)df);

  ld p;
  if(lower_tail){
    p = bm::cdf(dist, (ld)x);
  } else {
    p = bm::cdf(bm::complement(dist, (ld)x));
  }

  if(log_p){
    return (double)std::log(p);
  }
  else{
    return (double)p;
  }
  */
}

double qchisq_cpp(double p, double df, bool lower_tail, bool log_p){
  return R::qchisq(p,df,lower_tail ? 1 : 0, log_p ? 1:0);
  /*
  using ld = long double;
  if(!std::isfinite(df) || df <= 0){
    //Rcpp::Rcout << "[qchisq]_cpp invalid df" << df << std::endl;
    return details::NaN();
  }

  if(std::isnan(p)) return details::NaN();
  
  ld pp;
  if (log_p) {
    pp = std::exp((ld)p);
  } 
  else{
    pp = (ld)p;
  }

  if ( !lower_tail ) {
    pp = 1.0L - pp; // upper tail
  }
  
  if(pp <= 0.0L) return 0.0;
  if(pp >= 1.0L) return std::numeric_limits<double>::infinity();

  // Inverse of the chi-squared distribution
  bm::chi_squared_distribution dist((ld)df);
  ld quantile = bm::quantile(dist,pp);
  return (double)quantile;
  */
}




/*
// iteration 1
// qchisq
// standard normal density (log-pdf if logp=true)
double dnorm_cpp(double x, bool logp) {
  static const double SQRT2PI = std::sqrt(2 * M_PI);
  double u = -0.5 * x * x;
  if (logp) {
      return u - std::log(SQRT2PI);
  } else {
      return std::exp(u) / SQRT2PI;
  }
}

// standard normal CDF Φ(x) 
double pnorm_cpp(double x, bool lower_tail, bool log_p) {
  double c = 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
  if (!lower_tail) {
    c = 1.0 - c; // upper tail
  }

  if(log_p) {
    return std::log(c); // log probability
  }

  //else return probability
  return c;
}
double pnorm_cpp(double x, bool lower_tail, bool log_p) {
  const double y = 0.5 * std::erfc(-x * M_SQRT1_2); // M_SQRT1_2 is 1/sqrt(2)
  double p = lower_tail ? y : 1.0 - y;
  if(log_p){ return std::log(p); }
  return p;
}

//uses inverse of the standard normal distribution (quantile stuff so we introduce boost::math::quantile)
double qnorm_cpp(double p, bool lower_tail, bool log_p) {
  bm::normal_distribution<double> dist(0.0,1.0);
  if (log_p) {
    p = std::exp(p);
  } 
  if ( !lower_tail ) {
    p = 1.0 - p; // upper tail
  }
  double quantile = bm::quantile(dist, p);
  return quantile;
}

//these two below don't have log_p
// standard logistic density (log-pdf if logp=true)
double dlogis_cpp(double x, bool log_p) {
  // f(x) = exp(-x)/(1+exp(-x))^2
  double e = std::exp(-x);
  if (log_p) {
      return -x - 2.0 * std::log1p(e);
  } else {
      double denom = (1 + e);
      return e / (denom*denom);
  }
}

// standard logistic CDF F(x) = 1/(1+exp(-x))
double plogis_cpp(double x, bool lower_tail, bool log_p){
  double c = 1.0 / (1.0 + std::exp(-x));
  if(!lower_tail){
    c = 1.0 - c; // upper tail
  }
  if (log_p) {
    return std::log(c); // log probability branch
  }

  return c;
}

double pchisq_cpp(double x, double df, bool lower_tail, bool log_p) {
  //few base case checks to deal with special cases and to not crash the program incase infinite or zero values
  if(!std::isfinite(df) || df <= 0){
    Rcpp::Rcout << "[pchisq]_cpp invalid df" << df << "x= " << x << std::endl;
    return NA_REAL;
  }

  if(std::isnan(x)) return NA_REAL;

  if(x < 0.0){
    double prob = lower_tail ? 0.0 : 1.0;
    return log_p ? std::log(prob) : prob;
  }
  
  if(!std::isfinite(x)){
    double prob = lower_tail ? 1.0 : 0.0;
    return log_p ? std::log(prob) : prob;
  }

  bm::chi_squared dist(df);
  
  double c = bm::cdf(dist, x);
  if (!lower_tail) {
    c = 1.0 - c; // upper tail
  }
  if (log_p) {
    return std::log(c); // log probability branch
  }
  return c;
}

double qchisq_cpp(double p, double df, bool lower_tail, bool log_p) {
  //base case checks
  if(!std::isfinite(df) || df <= 0){
    Rcpp::Rcout << "[qchisq]_cpp invalid df" << df << std::endl;
    return NA_REAL;
  }

  if(std::isnan(p)) return NA_REAL;
  
  double q = log_p ? std::exp(p) : p;
  if(!lower_tail) q = 1.0 - q;

  if(q <= 0) return 0.0;
  if(q >= 1.0) return std::numeric_limits<double>::infinity();

  // Inverse of the chi-squared distribution
  bm::chi_squared dist(df);
  return bm::quantile(dist,q);
}
*/
//bottom three are copies of the copies of the survival package from utilties.cpp but made to work with std::vectors
// The following three utilities functions are from the survival package
int cholesky2_cpp(std::vector<std::vector<double>>& matrix, int n, double toler) {
  double temp;
  int i, j, k;
  double eps, pivot;
  int rank;
  int nonneg;

  nonneg = 1;
  eps = 0;
  for (i=0; i<n; i++) {
    if (matrix[i][i] > eps) eps = matrix[i][i];
  }
  if (eps==0) eps = toler; // no positive diagonals!
  else eps *= toler;

  rank = 0;
  for (i=0; i<n; i++) {
    pivot = matrix[i][i];
    if (std::isinf(pivot) == 1 || pivot < eps) {
      matrix[i][i] = 0;
      if (pivot < -8*eps) nonneg = -1;
    }
    else  {
      rank++;
      for (j=i+1; j<n; j++) {
        temp = matrix[i][j]/pivot;
        matrix[i][j] = temp;
        matrix[j][j] -= temp*temp*pivot;
        for (k=j+1; k<n; k++) matrix[j][k] -= temp*matrix[i][k];
      }
    }
  }

  return(rank*nonneg);
}


void chsolve2_cpp(std::vector<std::vector<double>>& matrix, int n, std::vector<double>& y) {
  int i, j;
  double temp;

  for (i=0; i<n; i++) {
    temp = y[i];
    for (j=0; j<i; j++)
      temp -= y[j]*matrix[j][i];
    y[i] = temp;
  }

  for (i=n-1; i>=0; i--) {
    if (matrix[i][i] == 0) y[i] = 0;
    else {
      temp = y[i]/matrix[i][i];
      for (j=i+1; j<n; j++)
        temp -= y[j]*matrix[i][j];
      y[i] = temp;
    }
  }
}


void chinv2_cpp(std::vector<std::vector<double>>& matrix, int n) {
  double temp;
  int i, j, k;

  for (i=0; i<n; i++){
    if (matrix[i][i] > 0) {
      matrix[i][i] = 1/matrix[i][i];   // this line inverts D
      for (j=i+1; j<n; j++) {
        matrix[i][j] = -matrix[i][j];
        for (k=0; k<i; k++)     // sweep operator
          matrix[k][j] += matrix[i][j]*matrix[k][i];
      }
    }
  }

  for (i=0; i<n; i++) {
    if (matrix[i][i] == 0) {  // singular row
      for (j=0; j<i; j++) matrix[i][j] = 0;
      for (j=i; j<n; j++) matrix[j][i] = 0;
    }
    else {
      for (j=i+1; j<n; j++) {
        temp = matrix[i][j]*matrix[j][j];
        matrix[j][i] = temp;
        for (k=i; k<j; k++)
          matrix[k][i] += temp*matrix[k][j];
      }
    }
  }
}


std::vector<std::vector<double>> invsympd_cpp(std::vector<std::vector<double>>& matrix, int n, double toler) {
  int i, j;
  std::vector<std::vector<double>> v = matrix;
  i = cholesky2_cpp(v, n, toler);
  chinv2_cpp(v, n);
  for (i=1; i<n; i++) {
    for (j=0; j<i; j++) {
      v[j][i] = v[i][j];
    }
  }

  return v;
}