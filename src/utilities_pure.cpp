#include <Rcpp.h>
//#include <R_ext/Applic.h>
#include "utilities_pure.h"

#include <vector>
#include <string>
#include <cmath>
#include <limits>
//these bottom three are introduced to replace the r method calls that use more complicated math functions
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/special_functions/erf.hpp>    // erfc, erfc_inv
#include <boost/math/special_functions/gamma.hpp>  // gamma_p/q and inverses


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

std::vector<int>
compose_block_rows(const std::vector<int>& row_idx,
                   int start,                 // idx[h]
                   const std::vector<int>& local_order) // order2 / order2x
{
  std::vector<int> out(local_order.size());
  for (size_t k = 0; k < local_order.size(); ++k) {
    out[k] = row_idx[start + local_order[k]];
  }
  return out;
}

inline double NaN() { return std::numeric_limits<double>::quiet_NaN(); }
// details to help out with the R math calls (needed replacing bc it caused issues in mt)
namespace details_pure {
  //inline double NaN() { return std::numeric_limits<double>::quiet_NaN(); }
  inline bool is_nan(double x) { return std::isnan(x); }
  inline bool is_finite(double x) { return std::isfinite(x); }

  // constants (long double for internal accuracy)
  constexpr long double LOG_SQRT_2PI = 0.9189385332046727417803297364056176398614L; // 0.5*log(2*pi)
  constexpr long double SQRT2 = 1.4142135623730950488016887242096980785697L;

  // Convert an input (p, lower_tail, log_p) into a tail probability without cancellation.
  // Returns probability in [0,1]. If invalid -> NaN.
  inline long double prob_from(double p, bool lower_tail, bool log_p) {
    if (is_nan(p)) return std::numeric_limits<long double>::quiet_NaN();

    long double pp;
    if (log_p) {
      // R semantics: log_p means p is log(probability)
      // exp(p) is safe; but if we need 1-exp(p), use -expm1(p).
      long double lp = static_cast<long double>(p);
      if (lp > 0.0L) return std::numeric_limits<long double>::quiet_NaN(); // log(prob) cannot be > 0
      if (lower_tail) {
        pp = std::expl(lp);
      } else {
        // upper tail prob = exp(lp) (already upper tail), but for conversions where we want lower tail
        // caller should decide; so here we interpret p as the requested tail prob already.
        pp = std::expl(lp);
      }
    } else {
      pp = static_cast<long double>(p);
    }

    if (!(pp >= 0.0L && pp <= 1.0L)) return std::numeric_limits<long double>::quiet_NaN();
    return pp;
  }

  // Compute 1 - exp(x) accurately for x near 0
  inline long double one_minus_exp(long double x) {
    // 1 - exp(x) = -expm1(x)
    return -std::expm1(x);
  }

  inline double log_prob_01(long double p) {
    if (p <= 0.0L) return -std::numeric_limits<double>::infinity();
    if (p >= 1.0L) return 0.0;
    return static_cast<double>(std::log(p));
  }
} // namespace details_pure
using namespace details_pure;

// creates a vector that sequentially adds from int start to end
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

// grabs however many are true inside of a bool vector
int sum_bool_cpp(const std::vector<bool>& vec) {
  int total = 0;
  for (size_t i = 0; i < vec.size(); ++i) {
    if(vec[i]) {
      total++;
    }
  }
  return total;
}

// calculates the mean in a vector
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

// sums the values in a double vector
double sumdouble_cpp(const std::vector<double>& vec) {
    double total = 0.0;
    for (const auto& val : vec) {
        total += val;
    }
    return total;
}

// sums the values in a int vector
int sumint_cpp(const std::vector<int>& vec) {
    double total = 0.0;
    for (const auto& val : vec) {
        total += val;
    }
    return total;
}

// gets the standard deviation in a vector
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
  
  if (is_nan(x)) return NaN();
  if (!is_finite(x)) return logp ? -std::numeric_limits<double>::infinity() : 0.0;

  long double xx = static_cast<long double>(x);
  long double logdens = -0.5L * xx * xx - LOG_SQRT_2PI;
  if (logp) return static_cast<double>(logdens);
  return static_cast<double>(std::expl(logdens));
}

double pnorm_cpp(double x, bool lower_tail, bool log_p){
  if (is_nan(x)) return NaN();

  if (x == -std::numeric_limits<double>::infinity()) {
    double p = lower_tail ? 0.0 : 1.0;
    return log_p ? (p == 0.0 ? -std::numeric_limits<double>::infinity() : 0.0) : p;
  }
  if (x == std::numeric_limits<double>::infinity()) {
    double p = lower_tail ? 1.0 : 0.0;
    return log_p ? (p == 0.0 ? -std::numeric_limits<double>::infinity() : 0.0) : p;
  }

  long double xx = static_cast<long double>(x);

  // lower tail: 0.5*erfc(-x/sqrt2)
  long double lt = 0.5L * boost::math::erfc(-xx / SQRT2);
  // upper tail: 0.5*erfc(+x/sqrt2)
  long double ut = 0.5L * boost::math::erfc( xx / SQRT2);

  if (!log_p) {
    return static_cast<double>(lower_tail ? lt : ut);
  }

  // log_p case: use complement for values near 1
  if (lower_tail) {
    // log(lt) stable when lt small; when lt ~ 1 use log1p(-ut)
    if (xx > 0.0L) {
      return static_cast<double>(std::log1p(-ut));
    } else {
      return static_cast<double>(std::log(lt));
    }
  } else {
    // log(ut) stable when ut small; when ut ~ 1 use log1p(-lt)
    if (xx < 0.0L) {
      return static_cast<double>(std::log1p(-lt));
    } else {
      return static_cast<double>(std::log(ut));
    }
  }
}

double qnorm_cpp(double p, bool lower_tail, bool log_p){
  if (is_nan(p)) return NaN();

  // Interpret p as the requested tail probability directly.
  // If log_p, p is log(tail_prob). If !log_p, p is tail_prob.
  long double tp;
  if (log_p) {
    long double lp = static_cast<long double>(p);
    if (lp > 0.0L) return NaN();
    tp = std::expl(lp);
  } else {
    tp = static_cast<long double>(p);
  }

  if (!(tp >= 0.0L && tp <= 1.0L)) return NaN();
  if (tp == 0.0L) return lower_tail ? -std::numeric_limits<double>::infinity()
                                    :  std::numeric_limits<double>::infinity();
  if (tp == 1.0L) return lower_tail ?  std::numeric_limits<double>::infinity()
                                    : -std::numeric_limits<double>::infinity();

  // Use tail-specific inversion via erfc_inv:
  // lower tail prob = 0.5*erfc(-x/sqrt2) => x = -sqrt2 * erfc_inv(2p)
  // upper tail prob = 0.5*erfc( x/sqrt2) => x =  sqrt2 * erfc_inv(2p)
  long double x;
  try {
    long double arg = 2.0L * tp;
    long double inv = boost::math::erfc_inv(arg);
    x = (lower_tail ? -SQRT2 : SQRT2) * inv;
  } catch (...) {
    return NaN();
  }
  return static_cast<double>(x);
}

double dlogis_cpp(double x, bool log_p){
  if (is_nan(x)) return NaN();
  if (!is_finite(x)) return log_p ? -std::numeric_limits<double>::infinity() : 0.0;

  // Stable versions (your commented formulas were good)
  if (log_p) {
    if (x >= 0.0) {
      return -x - 2.0 * std::log1p(std::exp(-x));
    } else {
      return  x - 2.0 * std::log1p(std::exp( x));
    }
  } else {
    if (x >= 0.0) {
      double e = std::exp(-x);
      double denom = 1.0 + e;
      return e / (denom * denom);
    } else {
      double e = std::exp(x);
      double denom = 1.0 + e;
      return e / (denom * denom);
    }
  }
}

double plogis_cpp(double x, bool lower_tail, bool log_p){
  if (is_nan(x)) return NaN();

  if (x == -std::numeric_limits<double>::infinity()) {
    double p = lower_tail ? 0.0 : 1.0;
    return log_p ? (p == 0.0 ? -std::numeric_limits<double>::infinity() : 0.0) : p;
  }
  if (x == std::numeric_limits<double>::infinity()) {
    double p = lower_tail ? 1.0 : 0.0;
    return log_p ? (p == 0.0 ? -std::numeric_limits<double>::infinity() : 0.0) : p;
  }

  if (log_p) {
    // Stable log(cdf) / log(1-cdf)
    if (lower_tail) {
      if (x >= 0.0) return -std::log1p(std::exp(-x));
      else          return  x - std::log1p(std::exp( x));
    } else {
      if (x <= 0.0) return -std::log1p(std::exp( x));
      else          return -x - std::log1p(std::exp(-x));
    }
  } else {
    // Stable cdf / ccdf
    if (lower_tail) {
      if (x >= 0.0) {
        double e = std::exp(-x);
        return 1.0 / (1.0 + e);
      } else {
        double e = std::exp(x);
        return e / (1.0 + e);
      }
    } else {
      if (x >= 0.0) {
        double e = std::exp(-x);
        return e / (1.0 + e);
      } else {
        double e = std::exp(x);
        return 1.0 / (1.0 + e);
      }
    }
  }
}

double pchisq_cpp(double x, double df, bool lower_tail, bool log_p){
  if (is_nan(x) || is_nan(df)) return NaN();
  if (!is_finite(df) || df <= 0.0) return NaN();

  if (x < 0.0) {
    double p = lower_tail ? 0.0 : 1.0;
    return log_p ? (p == 0.0 ? -std::numeric_limits<double>::infinity() : 0.0) : p;
  }
  if (x == std::numeric_limits<double>::infinity()) {
    double p = lower_tail ? 1.0 : 0.0;
    return log_p ? (p == 0.0 ? -std::numeric_limits<double>::infinity() : 0.0) : p;
  }

  long double a  = static_cast<long double>(df) / 2.0L;
  long double xx = static_cast<long double>(x)  / 2.0L;

  try {
    // Use gamma_q for complement stability when needed
    long double q = boost::math::gamma_q(a, xx); // upper tail
    if (!log_p) {
      long double p = lower_tail ? (1.0L - q) : q;
      return static_cast<double>(p);
    } else {
      if (lower_tail) return static_cast<double>(std::log1p(-q)); // log(1-q)
      else            return static_cast<double>(std::log(q));
    }
  } catch (...) {
    return NaN();
  }
}

double qchisq_cpp(double p, double df, bool lower_tail, bool log_p){
  if (is_nan(p) || is_nan(df)) return NaN();
  if (!is_finite(df) || df <= 0.0) return NaN();

  // p is the requested tail prob (lower or upper depending on lower_tail flag)
  long double tp;
  if (log_p) {
    long double lp = static_cast<long double>(p);
    if (lp > 0.0L) return NaN();
    tp = std::expl(lp);
  } else {
    tp = static_cast<long double>(p);
  }

  if (!(tp >= 0.0L && tp <= 1.0L)) return NaN();
  if (tp == 0.0L) return 0.0; // matches typical R behavior for qchisq(0, df)
  if (tp == 1.0L) return std::numeric_limits<double>::infinity();

  long double a = static_cast<long double>(df) / 2.0L;

  try {
    // Avoid cancellation: use tail-specific inverse
    long double xhalf = lower_tail
      ? boost::math::gamma_p_inv(a, tp)  // solves gamma_p(a, x) = tp
      : boost::math::gamma_q_inv(a, tp); // solves gamma_q(a, x) = tp (upper tail)

    long double x = 2.0L * xhalf;
    return static_cast<double>(x);
  } catch (...) {
    return NaN();
  }
}


// bottom three are copies of the copies of the survival package from utilties.cpp but made to work with std::vectors
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

//overloaded version using MatrixRM
int cholesky2_cpp(MatrixRM& m, int n, double toler) {
  double temp;
  int i, j, k;
  double eps, pivot;
  int rank;
  int nonneg;

  nonneg = 1;
  eps = 0.0;
  for (i = 0; i < n; i++) {
    if (m(i,i) > eps) eps = m(i,i);
  }
  if (eps == 0.0) eps = toler;
  else eps *= toler;

  rank = 0;
  for (i = 0; i < n; i++) {
    pivot = m(i,i);
    if (std::isinf(pivot) || pivot < eps) {
      m(i,i) = 0.0;
      if (pivot < -8*eps) nonneg = -1;
    } else {
      rank++;
      for (j = i+1; j < n; j++) {
        temp = m(i,j) / pivot;
        m(i,j) = temp;
        m(j,j) -= temp*temp*pivot;
        for (k = j+1; k < n; k++) {
          m(j,k) -= temp * m(i,k);
        }
      }
    }
  }

  return rank * nonneg;
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

// overloaded version with MatrixRM
void chsolve2_cpp(const MatrixRM& m, int n, std::vector<double>& y) {
  int i, j;
  double temp;

  // forward solve
  for (i = 0; i < n; i++) {
    temp = y[i];
    for (j = 0; j < i; j++) {
      temp -= y[j] * m(j,i);
    }
    y[i] = temp;
  }

  // backward solve
  for (i = n-1; i >= 0; i--) {
    if (m(i,i) == 0.0) y[i] = 0.0;
    else {
      temp = y[i] / m(i,i);
      for (j = i+1; j < n; j++) {
        temp -= y[j] * m(i,j);
      }
      y[i] = temp;
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

// overloaded version 
MatrixRM invsympd_cpp(const MatrixRM& matrix, int n, double toler) {
  MatrixRM v = matrix;              // copy
  (void)cholesky2_cpp(v, n, toler); // factorize in-place

  MatrixRM inv;
  inv.rows = n;
  inv.cols = n;
  inv.a.assign((size_t)n*(size_t)n, 0.0);

  std::vector<double> col(n);

  for (int c = 0; c < n; ++c) {
    std::fill(col.begin(), col.end(), 0.0);
    col[c] = 1.0;

    chsolve2_cpp(v, n, col);        // solves A x = e_c using factorization

    for (int r = 0; r < n; ++r) {
      inv(r,c) = col[r];
    }
  }

  // force symmetry (like your old code)
  for (int i = 1; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      inv(j,i) = inv(i,j);
    }
  }

  return inv;
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