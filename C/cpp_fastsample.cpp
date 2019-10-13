#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

struct Comp{
  Comp(const Rcpp::NumericVector& v ) : _v(v) {}
  bool operator ()(int a, int b) { return _v[a] < _v[b]; }
  const Rcpp::NumericVector& _v;
};

// [[Rcpp::export]]
Rcpp::NumericVector RcppSample(NumericVector n, IntegerVector size, NumericVector prob){
  int num = as<int>(size), x = as<int>(n);
  Rcpp::NumericVector vx = Rcpp::clone<Rcpp::NumericVector>(x);
  Rcpp::NumericVector pr = Rcpp::clone<Rcpp::NumericVector>(prob);
  Rcpp::NumericVector rnd = rexp(x) / pr;
  for(int i= 0; i<vx.size(); ++i) vx[i] = i;
  std::partial_sort(vx.begin(), vx.begin() + num, vx.end(), Comp(rnd));
  vx = vx[seq(0, num - 1)] + 1;
  return vx;
}

/*** R
#done compiling sampling algorithm
*/