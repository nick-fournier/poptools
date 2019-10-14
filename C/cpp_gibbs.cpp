//##
//# Gibbs sampler in C++ for R
//#
//# NOTICE:  All information, intellectual and technical concepts contained herein is,
//# and remains the property of Nicholas Marc Fournier. Dissemination of this 
//# information or reproduction of this material is strictly forbidden unless
//# prior written permission is obtained from Nicholas Marc Fournier
//##

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <progress.hpp>
#include <iostream>
#include <vector>
#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// Function to map they multidimensional matrix at vector: index = i + j*I + k*IJ ...
// [[Rcpp::export]]
int matkey(IntegerVector coordinates, IntegerVector dim, int base=1){
  int n = dim.size(), prod, idx = 0;
  IntegerVector coord = clone(coordinates);
  //adjusting coord vector to start at 0
  if(base != 0)
    for(int a=0; a < n; a++) {
      coord[a] = coord[a]-1;
    }
  //main loop
  for(int a=0; a < n; a++){
    prod = coord[a];
    for(int b=0; b < a; b++){
      prod = prod * dim[b];
    }
    idx = idx + prod;
  }
  return idx;
}

// Approximate equals
bool double_equals(double a, double b, double epsilon)  {
  return std::abs(a-b) < epsilon;
}

// Progress bar
class MinimalProgressBar: public ProgressBar{
public:
  MinimalProgressBar()  {
    _finalized = false;
  }
  
  ~MinimalProgressBar() {}
  
  void display() {
    REprintf("Progress: ");
  }
  
  void update(float progress) {
    if (_finalized) return;
    REprintf("+");
  }
  
  void end_display() {
    if (_finalized) return;
    REprintf("\n");
    _finalized = true;
  }
  
private:
  
  bool _finalized;
  
};

//Simple discrete Gibbs sampler
// [[Rcpp::export]]
IntegerVector RcppGibbs(int n, int burnin, int thn,
                       IntegerVector start,
                       List condtargets,
                       List condvars,
                       List condex,
                       List condprobs,
                       bool display_progress=true) {
  
  int h,i,k,c,t, it=0, maxit=n;
  IntegerMatrix mat(n, start.size());
  IntegerVector current = start, vardims(start.size());
  
  //MinimalProgressBar pb;
  Progress p(n+burnin, display_progress);
  
  //getting the dimensions of the variables
  for(i=0; i<condvars.size(); i++) {
    IntegerVector tmp = condvars[i];
    vardims[i] = tmp.size();
  }
  
  //Looping over i to n + burn in
  for (i=0; i<n+burnin; i++) {
    bool stop = false;
    //Thinning every thn'd
    for (h=0; h<thn; h++) {
      // Generating a candidate by looping over the conditionals for the variables
      for (c=0; c<condex.size(); c++) {
        // extracting the conditional probs, and making temporary prob vector
        NumericVector cprobs = condprobs[c];
        // Extracting the targets, indices, and dimensions.
        IntegerVector targets = condtargets[c], vardex = condex[c];
        // Initializing dimensions, and spreads.
        IntegerVector dims(vardex.size()), tmpspread(vardex.size()), spread(vardex.size());
        // Going through the targets, recycling conditionals if possible
        for (t=0; t<targets.size(); t++) {
          // Index of target overall
          int targ = targets[t], sprdtarg=0;
          IntegerVector varvals = condvars[targ];
          NumericVector tmp(vardims[targ]);
          // Getting index of target within the spread
          for (k=0; k<vardex.size(); k++) if ( targets[t] == vardex[k] ) sprdtarg = k;
          // Making a copy of current draw, and getting dimensions for matrix key finder
          for (k=0; k<vardex.size(); k++) {
            int idx = vardex[k];
            tmpspread[k] = current[idx];
            dims[k] = vardims[idx];
          }
          // Making deep copy
          spread = clone(tmpspread);
          //Varying target variable in spread, getting the probability index each time.
          for (k=0; k<vardims[targ]; k++) {
            spread[sprdtarg] = k+1; //+1 because indexed at 1, which is -1 to 0 in matkey()
            int key = matkey(spread, dims);
            tmp[k] = cprobs[key];
          }
          // Making deep copy of probability vector because "sample" alters it in memory
          NumericVector probs = clone(tmp);
          double sum = 0;
          //Check if dead end, break to last good draw
          for (k=0; k<probs.size(); k++) sum += probs[k];
          if (sum == 0 ) {
            c = 0;
            i -= 2;
            it++;
            stop = true;
            break;
          } else {
            it=0;
            current[targ] = Rcpp::RcppArmadillo::sample(varvals, 1, TRUE, probs)[0];
          }
        }
        //Break conditional loop
        if (stop == true) break;
      }
      //Break from thinning loop too
      if (stop == true) break;
    }
    if (stop != true) p.increment();
    if (it > maxit ) {
      Rcout << "Stuck in dead end" << std::endl;
      break;
      }
    if(i>=burnin) mat(i-burnin,_) = current;
  }
  return(mat);
}


