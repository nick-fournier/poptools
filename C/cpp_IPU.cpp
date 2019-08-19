//##
//# Iterative Proportional Updating written in C++ for R
//#
//# NOTICE:  All information, intellectual and technical concepts contained herein is,
//# and remains the property of Nicholas Marc Fournier. Dissemination of this 
//# information or reproduction of this material is strictly forbidden unless
//# prior written permission is obtained from Nicholas Marc Fournier
//##

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

vec colSums(const mat & X){
  int nCols = X.n_cols;
  vec out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i));
  }
  return(out);
}

//dense
double deltacalc(NumericMatrix Amat, NumericVector ccons, NumericVector Wts){
  int c = Amat.ncol();
  int r = Amat.nrow();
  double sums = 0, Wsums = 0;
  for(int j=0; j<c; j++){
    Wsums=0;
    for(int i=0; i<r; i++){
      Wsums = Wsums + Amat(i,j) * Wts[i];
    }
    Wsums = fabs(Wsums - ccons[j])/ccons[j];
    sums = sums + Wsums;
  }
  sums = sums/c;
  return(sums);
}
double errormarg(NumericMatrix Amat, NumericVector Wts, NumericVector ccons){
  int c = Amat.ncol();
  int r = Amat.nrow();
  double max = 0;
  double error = 0;
  NumericVector colsums(c);
  for(int j=0; j<c; j++){
    for(int i=0; i<r; i++){
      colsums[j] = colsums[j] + Amat(i,j) * Wts[i];
    }
    error = fabs(ccons[j] - colsums[j]);
    if(error > max) {
      max = error;
    }
  }
  return(max);
}
double errorperc(NumericMatrix Amat, NumericVector Wts, NumericVector ccons){
  int c = Amat.ncol();
  int r = Amat.nrow();
  double max = 0;
  double error = 0;
  NumericVector colsums(c);
  for(int j=0; j<c; j++){
    for(int i=0; i<r; i++){
      colsums[j] = colsums[j] + Amat(i,j) * Wts[i];
    }
    error = fabs(ccons[j] - colsums[j]) / ccons[j];
    if(error > max) {
      max = error;
    }
  }
  return(max);
}
double srmsecalc(NumericMatrix Amat, NumericVector Wts, NumericVector ccons){
  int c = Amat.ncol();
  int r = Amat.nrow();
  double error = 0;
  double esum = 0;
  double srmse = 0;
  double ybar = 0;
  NumericVector colsums(c);
  for(int j=0; j<c; j++){
    for(int i=0; i<r; i++){
      colsums[j] = colsums[j] + Amat(i,j) * Wts[i];
    }
    error = fabs(ccons[j] - colsums[j]); //raw error
    error = pow(error, 2); //square error
    esum = esum + error; //sum of error
    ybar = ybar + ccons[j];
  }
  ybar = ybar / c;
  srmse = sqrt((esum / c)) / ybar;
  return(srmse);
}
NumericVector newWeightCalc(NumericVector Acol, double con, NumericVector Wts){
  int N = Acol.size();
  double sum=0;
  //Adjust column and gather sum
  for(int i=0; i<N; i++){
    Acol[i] = Acol[i] * Wts[i];
    sum = sum + Acol[i];
  }
  //resulting new adjustment ratio
  double ratio = con/sum;
  //Update the weights for non-empty cells in column
  for(int i=0; i<N; i++){
    if(Acol[i] > 0){
      Wts[i] = ratio * Wts[i];
    }
  }
  return(Wts);
}

//sparse
double deltacalc_s(const sp_mat& Amat, NumericVector ccons, NumericVector Wts, const umat& nz){
  int r = nz.n_cols;
  int C = Amat.n_cols, i = 0, j = 0, jprev = 0;
  double sums = 0, wsums = 0, ratio = 0;
  for(int k=0; k<r; k++){
    i = nz(0,k);
    j = nz(1,k);
    if(jprev == j) {
      wsums = wsums + Amat(i,j) * Wts[i];
    } else {
      ratio = fabs(wsums - ccons[jprev])/ccons[jprev];
      jprev = j;
      wsums = Amat(i,j) * Wts[i];
      sums = ratio + sums;
    }
  }
  // Final calculation
  ratio = fabs(wsums - ccons[j])/ccons[j];
  sums = ratio + sums;
  sums = sums/C;
  return(sums);
}
double errormarg_s(const sp_mat& Amat, NumericVector Wts, NumericVector ccons, const umat& nz){
  int r = nz.n_cols, i = 0, j = 0, jprev = 0;
  double max = 0, error = 0, csums = 0;
  for(int k=0; k<r; k++){
    i = nz(0,k);
    j = nz(1,k);
    if(jprev == j) {
      csums = csums + Amat(i,j) * Wts[i];
    } else {
      error = fabs(ccons[jprev] - csums);
      jprev = j;
      csums = Amat(i,j) * Wts[i];
      if(error > max) {
        max = error;
      }
    }
  }
  // Final calculation
  error = fabs(ccons[j] - csums);
  if(error > max) {
    max = error;
  }
  return(max);
}
double errorperc_s(const sp_mat& Amat, NumericVector Wts, NumericVector ccons, const umat& nz){
  int r = nz.n_cols, i = 0, j = 0, jprev = 0;
  double max = 0, error = 0, csums = 0;
  for(int k=0; k<r; k++){
    i = nz(0,k);
    j = nz(1,k);
    if(jprev == j) {
      csums = csums + Amat(i,j) * Wts[i];
    } else {
      error = fabs(ccons[jprev] - csums) / ccons[jprev];
      jprev = j;
      csums = Amat(i,j) * Wts[i];
      if(error > max) {
        max = error;
      }
    }
  }
  // Final calculation
  error = fabs(ccons[j] - csums) / ccons[j];
  if(error > max) {
    max = error;
  }
  return(max);
}
double srmsecalc_s(const sp_mat& Amat, NumericVector Wts, NumericVector ccons, const umat& nz){
  int c = Amat.n_cols, r = nz.n_cols, i = 0, j = 0, jprev = 0;
  double error = 0, esum = 0, srmse = 0, ybar = 0, csums = 0;
  for(int k=0; k<r; k++){
    i = nz(0,k);
    j = nz(1,k);
    if(jprev == j) {
      csums = csums + Amat(i,j) * Wts[i];
    } else {
      error = fabs(ccons[jprev] - csums); //raw error
      error = pow(error, 2);          //square error
      esum = esum + error;            //sum of error
      ybar = ybar + ccons[jprev];
      
      jprev = j;
      csums = Amat(i,j) * Wts[i];
    }
  }
  // Final calculation
  error = fabs(ccons[jprev] - csums); //raw error
  error = pow(error, 2);          //square error
  esum = esum + error;            //sum of error
  ybar = ybar + ccons[jprev];
  
  ybar = ybar / c;
  srmse = sqrt((esum / c)) / ybar;
  return(srmse);
}
NumericVector newWeightCalc_s(const sp_mat& Amat, double con, NumericVector Wts, int colidx, const umat& nz){
  umat nzcol = nz.row(1);
  uvec nzrows = find(nzcol == colidx);
  int N = nzrows.n_elem, idx = 0;
  NumericVector Acol(N);
  double sum=0;
  //Gather sum of adjusted column
  for(int i=0; i<N; i++){
    idx = nz(0,nzrows[i]);
    Acol[i] = Amat(idx,colidx) * Wts[idx];
    sum = sum + Acol[i];
  }
  //resulting new adjustment ratio
  double ratio = con/sum;
  //Update the weights for non-empty cells in column
  for(int i=0; i<N; i++){
    idx = nz(0,nzrows[i]);
    Wts[idx] = ratio * Wts[idx];
  }
  return(Wts);
}
arma::umat get_locations(const sp_mat& B) {
  // Make const iterator
  arma::sp_mat::const_iterator start = B.begin();
  arma::sp_mat::const_iterator end   = B.end();
  // Calculate number of points
  int n = std::distance(start, end);
  // Kill process if no values are found (very sparse matrix)
  if (n <= 0) { Rcpp::stop("No values found!"); }
  // Build a location storage matrix
  arma::umat locs(2, n);
  // Create a vector to store each row information in. (Row, Col)
  arma::uvec temp(2);
  // Start collecting locations
  arma::sp_mat::const_iterator it = start; 
  for(int i = 0; i < n; ++i)
  {
    temp(0) = it.row();
    temp(1) = it.col();
    locs.col(i) = temp;
    ++it; // increment
  }
  return locs;
}

// [[Rcpp::export]]
NumericVector IPU_cpp(NumericMatrix A,
                      NumericVector Cons,
                      NumericVector HHidx,
                      double corncrit = 1,
                      double crit = 1e-8,
                      int maxit = 1000,
                      bool print = false) {
  int nR = A.nrow();
  int nC = A.ncol();
  int nH = HHidx.size();
  int idx = 0;
  int iters = 0;
  //initial weights vector
  NumericVector W(nR, 1.0);
  //new weights vector
  NumericVector newW(nR);
  double Diff = 1;
  double delta = 0;
  double deltaprev = deltacalc(A, Cons, W); //initial delta value
  double deltamin = deltaprev;
  
  for(int i=0; i<=maxit; i++) {
    iters = i;
    for(int j=0; j<nC; j++) {
      newW = newWeightCalc(A(_,j), Cons[j], W);
    }
    deltaprev = delta ;
    //Calcualte new goodness-of-fit using new weights
    delta = deltacalc(A, Cons, newW);
    //Improvement of goodness-of-fit
    Diff = fabs(delta - deltaprev);
    //if the new weights are better, then we replace the weights and update min
    if(delta < deltamin) {
      deltamin = delta;
      W = newW;
    }
    if(Diff < crit) break;
    if(print) Rcout << "Iteration: " << i << " Delta diff: " << Diff << " Delta: " << deltamin << std::endl;
  }
  //if a perfect solution is not found, then find a corner solution for HH's by running once more for just HH constraints
  if(delta > corncrit) {
    if(print) Rcout << "Corner solution not found, finding corner for households" << std::endl;
    for(int j=0; j<nH; j++) {
      idx = HHidx[j];
      W = newWeightCalc(A(_,idx), Cons[idx], W);
    }
  }
  
  double merror = errormarg(A, W, Cons);
  double perror = errorperc(A, W, Cons);
  double rmse = srmsecalc(A, W, Cons);
  if(print) {
    Rcout << "Finished after " << iters << " iterations" << std::endl;
    Rcout << "Final Delta difference criterion: " << Diff << std::endl;
    Rcout << "Final small Delta: " << deltamin << std::endl;
    Rcout << "Max error margin: " << merror << std::endl;
    Rcout << "Max percent error margin: " << perror << std::endl;
    Rcout << "RMSE: " << rmse << std::endl;
  }
  
  return(W);
}
// [[Rcpp::export]]
NumericVector IPU_cpp_sparse(const sp_mat& A,
                             NumericVector Cons,
                             NumericVector HHidx,
                             double corncrit = 1,
                             double crit = 1e-8,
                             int maxit = 1000,
                             bool print = false) {

  int nR = A.n_rows, nC = A.n_cols, nH = HHidx.size(), idx = 0, iters = 0;
  //get nonzero elements
  arma::umat nonzero = get_locations(A);
  //initial weights vector
  NumericVector W(nR, 1.0);

  //new weights vector
  NumericVector newW(nR);
  double Diff = 1;
  double delta = 0;
  double deltaprev = deltacalc_s(A, Cons, W, nonzero); //initial delta value
  double deltamin = deltaprev;
  for(int i=0; i<=maxit; i++) {
    iters = i;
    for(int j=0; j<nC; j++) {
      newW = newWeightCalc_s(A, Cons[j], W, j, nonzero);
    }
    deltaprev = delta ;
    //Calculate new goodness-of-fit using new weights
    delta = deltacalc_s(A, Cons, newW, nonzero);
    //Improvement of goodness-of-fit
    Diff = fabs(delta - deltaprev);
    //if the new weights are better, then we replace the weights and update min
    if(delta < deltamin) {
      deltamin = delta;
      W = newW;
    }
    if(Diff < crit) break;
    if(print) Rcout << "Iteration: " << i << " Delta diff: " << Diff << " Delta: " << deltamin << std::endl;
  }
  //if a perfect solution is not found, then find a corner solution for HH's by running once more for just HH constraints
  if(delta > corncrit) {
    if(print) Rcout << "Corner solution not found, finding corner for households" << std::endl;
     for(int j=0; j<nH; j++) {
       idx = HHidx[j];
       W = newWeightCalc_s(A, Cons[idx], W, idx, nonzero);
     }
  }

  double merror = errormarg_s(A, W, Cons, nonzero);
  double perror = errorperc_s(A, W, Cons, nonzero);
  double rmse = srmsecalc_s(A, W, Cons, nonzero);
  if(print) {
    Rcout << "Finished after " << iters << " iterations" << std::endl;
    Rcout << "Final Delta difference criterion: " << Diff << std::endl;
    Rcout << "Final small Delta: " << deltamin << std::endl;
    Rcout << "Max error margin: " << merror << std::endl;
    Rcout << "Max percent error margin: " << perror << std::endl;
    Rcout << "RMSE: " << rmse << std::endl;
  }
  return(W);
}
/*** R
#done compiling
#Parameters for simple test
# library(Matrix)
# A <- matrix( c(1, 0, 1, 1, 1,  1, 0, 1, 0, 1,  1, 0, 2, 1, 0,  0, 1, 1, 0, 2,
#               0, 1, 0, 2, 1,  0, 1, 1, 1, 0,  0, 1, 2, 1, 2,  0, 1, 1, 1, 0), ncol=5, byrow = T)
# spA <- as(A, 'sparseMatrix')
# C <- c(35.00, 65.00, 91.00, 65.00, 104.00)
#system.time(IPU_cpp_sparse(NA, spA, C, HHidx=1:2, corncrit = 0.1, crit = 1e-8, print=T))
#system.time(IPU_cpp(NA, A, C, HHidx=1:2, crit = 1e-8, print=F))
*/