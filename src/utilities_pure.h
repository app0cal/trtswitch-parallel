#pragma once
#include <Rcpp.h>
#include <R_ext/Applic.h>

#include <vector>
#include <string>
#include <cmath>
//these bottom three are introduced to replace the r method calls that use more complicated math functions
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/chi_squared.hpp>

//all functions here are just the signatures so no code is needed
typedef std::vector<std::vector<double>> matrix;

namespace bm = boost::math;
std::vector<int> seq_cpp(
    int start, 
    int end
);
/*
template<typename T>
void reorder(
    std::vector<T>& vec, 
    const std::vector<int>& order
);
*/

matrix subset_matrix_by_row_cpp(
    const matrix& mat,
    const std::vector<int>& order
);

template<typename T>
void reorder(std::vector<T>& vec, const std::vector<int>& order) {
  std::vector<T> old = vec;
  for (size_t i = 0; i < order.size(); ++i) {
    vec[i] = old[order[i]];
  }
}
/*
template<typename Container>
std::vector<int> which_true(
    const Container& flags
);

template<class T1>
std::vector<int> which_single_cpp(
    const std::vector<T1>& a, 
    const T1& val1
);

template<class T1, class T2>
std::vector<int> which_dual_cpp(
    const std::vector<T1>& a,
    const T1& val1,
    const std::vector<T2>& b, 
    const T2& val2
);


template<typename T>
std::vector<T> subset_by_idx(
    const std::vector<T>& vec,
    const std::vector<int>& idx
);
*/
//introduced single and dual which functions to avoid Rcpp API calls
// single is for single vector, dual is for two vectors
// this introduces a more generic way to find indices of elements in vectors
// if three or more vectors are needed, we can extend this further 

//use this for bool vectors and integer vectors w only 1's and 0's representing true/false
template<typename Container>
std::vector<int> which_true(const Container& flags) {
  std::vector<int> idx;
  idx.reserve(flags.size());
  for (std::size_t i = 0; i < flags.size(); ++i) {
    if ( static_cast<bool>(flags[i]) ) {
      idx.push_back(static_cast<int>(i));
    }
  }
  return idx;
}

template<class T1>
std::vector<int> which_single_cpp(const std::vector<T1>& a, const T1& val1) {
  std::vector<int> idx;
  idx.reserve(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] == val1) {
      idx.push_back(static_cast<int>(i));
    }
  }
  return idx;
}

template<class T1, class T2>
std::vector<int> which_dual_cpp(const std::vector<T1>& a, const T1& val1,
                           const std::vector<T2>& b, const T2& val2) {
  std::vector<int> idx;
  idx.reserve(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i] == val1 && b[i] == val2) {
      idx.push_back(static_cast<int>(i));
    }
  }
  return idx;
}

//subsets a vector by indices
template<typename T>
std::vector<T> subset_by_idx(const std::vector<T>& vec,
                             const std::vector<int>& idx) {
  std::vector<T> out;
  out.reserve(idx.size());
  for (int i : idx) {
    out.push_back(vec[i]);
  }
  return out;
}

int sum_bool_cpp(
    const std::vector<bool>& vec
);

double mean_cpp(
    const std::vector<double>& vec
);

double sum_cpp(
    const std::vector<double>& vec
);

double sd_cpp(
    const std::vector<double>& vec
);

double dnorm_cpp(
    double x, 
    bool logp=false
);

double pnorm_cpp(
    double x, 
    bool lower_tail=true, 
    bool log_p=false
);

double qnorm_cpp(
    double p, 
    bool lower_tail=true, 
    bool log_p=false
);

double dlogis_cpp(
    double x, 
    bool log_p=false
);

double plogis_cpp(
    double x, 
    bool lower_tail=true, 
    bool log_p=false
);

double pchisq_cpp(
    double x, 
    double df, 
    bool lower_tail=true, 
    bool log_p=false
);

double qchisq_cpp(
    double p, 
    double df, 
    bool lower_tail=true, 
    bool log_p=false
);

int cholesky2_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n, 
    double toler
);

void chsolve2_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n, 
    std::vector<double>& y
);

void chinv2_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n
);

std::vector<std::vector<double>> invsympd_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n, 
    double toler
);