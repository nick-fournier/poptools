#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
IntegerVector dupes(const sp_mat& A) {
  
  int nc = A.n_cols, nr = A.n_rows;
  IntegerVector x(nr);
  for(int i=0; i < nr; i++) {
    Rcout << i+1 << " of " << nr << std::endl;
    if(x(i) == 1) break;
    for(int j=0; j < nr; j++) {
      for(int k=0; k < nc; k++){
        if(A(i,k) != A(j,k)) break;
        if(k == nc-1 && i != j) x(j) = 1;
      }
    }
  }
  return x;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#compiled
*/
