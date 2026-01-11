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

struct MatrixRM {
  //Row Major Matrix (R uses Column major but Row Major works pretty well for our code)
  int rows=0, cols=0;
  std::vector<double> a;
  int nrows() const { return rows; }
  int ncols() const { return cols; }
  double& operator()(int r,int c)       { return a[(size_t)r*cols + c]; }
  double  operator()(int r,int c) const { return a[(size_t)r*cols + c]; }

  MatrixRM() = default;

  MatrixRM(int r, int c){
    if(r < 0|| c < 0){
      throw std::invalid_argument("Matrix creatoin failed, invalid dimensions.");
    }
    const size_t rr = (size_t)r;
    const size_t cc = (size_t)c;
    rows = r;
    cols = c;
    a.assign(rr*cc,0.0); //allocates to 0
  }

  MatrixRM(int r, int c, double init_value){
    if(r < 0|| c < 0){
      throw std::invalid_argument("Matrix creatoin failed, invalid dimensions.");
    }
    const size_t rr = (size_t)r;
    const size_t cc = (size_t)c;
    rows = r;
    cols = c;
    a.assign(rr*cc,init_value); //allocates to whatever value is passed
  }
};

// expensive way to permute a matrix based off a order vector, but used for correctness
MatrixRM permute_rows(const MatrixRM& M, const std::vector<int>& order);

std::vector<int>
compose_block_rows(const std::vector<int>& row_idx,
                   int start,                 // idx[h]
                   const std::vector<int>& local_order); // order2 / order2x

struct MatrixRMView {
  int rows = 0, cols = 0;
  const double* a = nullptr; // non-owning

  int nrows() const { return rows; } // optional
  int ncols() const { return cols; } // optional

  double operator()(int r, int c) const {
    return a[(size_t)r * cols + c];
  }
};

// meant to be non owning version of row index view for matrix rm view
// this will be changed to be w templates instead later
struct RowIndexViewView {
  const MatrixRMView* base = nullptr;
  const std::vector<int>* idx = nullptr;

  int rows() const { return (int)idx->size(); }
  int cols() const { return base->cols; }

  double operator()(int r, int c) const {
    return (*base)((*idx)[r], c);
  }
};

struct RowIndexView {
  const MatrixRM* base = nullptr;
  const std::vector<int>* idx = nullptr;

  int rows() const { return (int)idx->size(); }
  int cols() const { return base->ncols(); }

  double operator()(int r,int c) const {
    return (*base)((*idx)[r], c);
  }
};


template <class ZMat>
struct RowRangeViewT {
  const ZMat* base = nullptr;
  int start = 0, len = 0;

  int rows() const { return len; }
  int cols() const { return base->cols(); }

  double operator()(int r, int c) const {
    return (*base)(start + r, c);
  }
};


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
/*
std::vector<std::vector<double>> subset_matrix_by_row_cpp(
    const std::vector<std::vector<double>>& mat,
    const std::vector<int>& order
); */

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

double sumdouble_cpp(
    const std::vector<double>& vec
);

int sumint_cpp(
    const std::vector<int>& vec
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

int cholesky2_cpp(MatrixRM& m, int n, double toler);

void chsolve2_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n, 
    std::vector<double>& y
);

void chsolve2_cpp(const MatrixRM& m, int n, std::vector<double>& y);

void chinv2_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n
);

std::vector<std::vector<double>> invsympd_cpp(
    std::vector<std::vector<double>>& matrix, 
    int n, 
    double toler
);

MatrixRM invsympd_cpp(
  const MatrixRM& matrix, 
  int n, 
  double toler
);