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
#include "multinom.h"
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

// Getting the log likelihood
// [[Rcpp::export]]
double MLR(IntegerVector draws, NumericVector probs, IntegerVector vars) {
  
  int i,j;
  double LR = 0;
  NumericVector pi;
  IntegerVector v;
  //First remove zero probabilities
  for(j=0; j<vars.size(); j++) {
    if( probs[j] > 0 ) {
      pi.push_back(probs[j]);
      v.push_back(vars[j]);
    }
  }
  
  NumericVector x(v.size(), 1e-16);
  
  for(j=0; j<v.size(); j++) {
    for(i=0; i<draws.size(); i++) {
      if(draws[i] == v[j] ) x[j] += 1; //add one if matched
    }
  }
  
  //LR test
  for(j=0; j<v.size(); j++) LR += x[j] * log(draws.size() * pi[j] / x[j] ) ;
  LR = -2*LR;
  
  //if ( isfinite(LR) == 0 ) LR = -1e16;
  return LR;
}

// Checking the change in proportions of current size
NumericVector check_prop(IntegerMatrix mtx, List varlist)  {
  int nvars=0, idx=0, i,j,k,h;
  double sum=0;
  IntegerVector vardims(varlist.size());

  //getting the dimensions of the variable
  for(j=0; j<varlist.size(); j++) {
    IntegerVector tmp = varlist[j];
    vardims[j] = tmp.size();
    nvars = nvars + vardims[j];
  }
  //Setting up sums vector
  NumericVector varsums(nvars);
  //Going through matrix and summing up variable frequencies
  for(j=0; j<mtx.ncol(); j++) {
    IntegerVector vars = varlist[j];
      for(k=0; k<vars.size(); k++) {
        for(i=0; i<mtx.nrow(); i++) {
          if(mtx(i,j) == vars[k] ) {
            //If matrix element matches var, add one.
            varsums[idx] = varsums[idx]+1;
            sum = sum + 1;
          }
        }
        idx = idx + 1;
      }
  }
  //taking the proportion
  for(h=0; h<varsums.size(); h++) {
    varsums[h] = varsums[h]/sum;
  }
  return varsums;
}

// Checking the change in proportions for an expected distribution size N
NumericVector check_nprop(IntegerMatrix mtx, List varlist, double n)  {
  int nvars=0, idx=0, i,j,k,h;
  IntegerVector vardims(varlist.size());
  
  //getting the dimensions of the variable
  for(j=0; j<varlist.size(); j++) {
    IntegerVector tmp = varlist[j];
    vardims[j] = tmp.size();
    nvars = nvars + vardims[j];
  }
  //Setting up sums vector
  NumericVector varsums(nvars);
  //Going through matrix and summing up variable frequencies
  for(j=0; j<mtx.ncol(); j++) {
    IntegerVector vars = varlist[j];
    for(k=0; k<vars.size(); k++) {
      for(i=0; i<mtx.nrow(); i++) {
        if(mtx(i,j) == vars[k] ) {
          //If matrix element matches var, add one.
          varsums[idx] = varsums[idx]+1;
        }
      }
      idx = idx + 1;
    }
  }
  //taking the proportion
  for(h=0; h<varsums.size(); h++) {
    varsums[h] = varsums[h]/n;
  }
  return varsums;
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

//Conditional Gibbs sampler
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

// Metropolis-Hastings
// [[Rcpp::export]]
IntegerVector RcppMH(int n, int burnin, int thn,
                        IntegerVector start,
                        List condtargets,
                        List condvars,
                        List condex,
                        List condprobs,
                        bool display_progress=true) {
  
  int h,i,k,c,t, it=0, maxit=burnin;
  IntegerMatrix mat(n+burnin, start.size()), outmat(n, start.size());
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
            it++;
            stop = true;
            break;
          } else {
            it=0;
            //the draw
            current[targ] = Rcpp::RcppArmadillo::sample(varvals, 1, TRUE, probs)[0];
            
            //Metroplis Hasting acceptance criteria
            if(i>=burnin) {
              //Adding to vector and checking likelihood
              IntegerVector vardraws(i);
              for (k=0; k<(vardraws.size() - 1); k++) vardraws[k] = mat(k,targ);
              vardraws[vardraws.size()] = current[targ];
              
              double lr_old = MLR(mat(_,targ), probs, varvals);
              double lr_new = MLR(vardraws, probs, varvals);
              double a = lr_new / lr_old;
              
              //If no improvement, decide whether to keep...
              if( a < 1 ) {
                IntegerVector coin = IntegerVector::create(0,1); // {0=reject, 1=keep}
                NumericVector accept = NumericVector::create(1-a,a);
                int cointoss = Rcpp::RcppArmadillo::sample(coin, 1, TRUE, accept)[0];
                if (cointoss == 0) {
                  it++;
                  stop = true;
                  break;
                }
              }
            }
          }
        }
        //Break conditional loop
        if (stop == true) break;
      }
      //Break from thinning loop too
      if (stop == true) break;
    }
    if (stop != true) {
      p.increment();
      mat(i,_) = current;
    } else {
      current = mat(i-1,_);
      i -= 2;
    }
    if (it > maxit ) {
      Rcout << "Stuck in dead end" << std::endl;
      break;
    }
  }
  
  // for (i=0; i < n; i++) {
  //   outmat(i,_) = mat(i+burnin,_);
  // }
  
  return(mat);
}


