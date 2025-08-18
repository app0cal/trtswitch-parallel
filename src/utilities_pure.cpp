#include <Rcpp.h>
#include <R_ext/Applic.h>
#include "utilities_pure.h"

#include <vector>
#include <string>
#include <cmath>
//these bottom three are introduced to replace the r method calls that use more complicated math functions
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/chi_squared.hpp>

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

typedef std::vector<std::vector<double>> matrix;


std::vector<int> seq_cpp(int start, int end) {
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

matrix subset_matrix_by_row_cpp(
  const matrix& mat, 
  const std::vector<int>& order
  ) {
  int n = static_cast<int>(order.size());
  int p = static_cast<int>(mat[0].size());
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
  // Inverse of the chi-squared distribution
  bm::chi_squared dist(df);

  if (log_p) {
    p = std::exp(p);
  } 

  if ( !lower_tail ) {
    p = 1.0 - p; // upper tail
  }
  double quantile = bm::quantile(dist, p);
  return quantile;
}

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
