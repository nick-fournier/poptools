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

//Binomial coefficient
int binomial(int n, int k)
{
  int num, den ;
  if ( n < k ) 
  {
    return(0) ; 
  }
  else 
  {
    den = 1;
    num = 1 ; 
    for (int i =  1  ; i <= k   ; i = i+1)
      den =    den * i;
    for (int j = n-k+1; j<=n; j=j+1)	
      num = num * j;
    return(num/den);
  } 
}

// multinomial coefficient
double multinomCoeff(NumericVector x){
  int nx = x.size();
  multinomial::SVI v(nx);
  int i;
  
  for(i = 0; i < nx; i++){
    v.at(i) = x[i];
  }
  
  double u = multinomial::multi<double>(v);

  return u;
}

// Getting the log likelihood
double loglikFunc(IntegerMatrix mtx, List varlist, NumericVector probs)  {
  int nvars=0, idx=0, i,j,k,h;
  double loglik;
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
        if(mtx(i,j) == vars[k] ) varsums[idx] = varsums[idx]+1; //add one if matched
      }
      idx = idx + 1; //next index key
    }
  }
  
  //Calculate multinom coefficient
  loglik = multinomCoeff(varsums);
  loglik = log(loglik);
  
  //Calculate log likelihood summation
  for(h=0; h<varsums.size(); h++) {
    loglik = loglik + varsums[h]*log(probs[h]);
  }
  return loglik;
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

//Single draw Gibbs sampling generator
IntegerVector GibbsGen(IntegerVector start,
                       List condvars,
                       List condex,
                       List condprobs) {
  
  int i,j,k,l,c, key;
  IntegerVector current = start, seq(condex.size()), vardims(condex.size()), bin = {0,1};
  NumericVector errorvec(1), jumps(1), P(2), acceptprobs(2);

  //getting the dimensions of the variable
  for(i=0; i<condex.size(); i++) {
    seq[i] = i;
    IntegerVector tmp = condvars[i];
    vardims[i] = tmp.size();
  }
  
  // Generating a candidate by looping over the conditionals for the variables
  for (j=0; j<condex.size(); j++) {
    c = seq[j]; //this is incase we want to randomize the conditional order
    // extracting the conditional probs, variables, and dimensions vectors
    NumericVector cprobs = condprobs[c], tmp(vardims[c]);
    IntegerVector vardex = condex[c], vars = condvars[c], conddims(condex.size()), tmpspread(condex.size()), spread(condex.size());
    
    // Making a copy of the last drawn variables to make a spread from
    for(k=0; k<vardex.size(); k++) {
      l = vardex[k];
      tmpspread[k] = current[l];
      conddims[k] = vardims[l];
    }
    
    // Making deep copy
    spread = clone(tmpspread);
    // Making a temporary vector of variable to clone the probabilities fro
    // looping within a variable to get the probability range
    for(k=0; k<vardims[c]; k++) {
      // varying the target variable
      spread[0] = k+1; //+1 because indexed at 1, which is -1 to 0 in matkey()
      key = matkey(spread, conddims);
      tmp[k] = cprobs[key];
    }
    
    // Making copy of probability vector b/c sample alters it into conditional probs
    NumericVector probs = clone(tmp);
    // Checking if all probabilities are 0
    double sum=0;
    for(k=0; k<vardims[c]; k++)
      sum = sum + probs[k];
    
    // if there are no probabilities, reject and reset to entirely new start condition
    if ( double_equals(sum, 0, 1e-8) ) {
      //restart from from beginning with different sequence
      current = start;
      NumericVector probs(seq.size(),1);
      seq = Rcpp::RcppArmadillo::sample(seq, 1, TRUE, probs)[0];
      j=0;
    } else {
      current[c] = Rcpp::RcppArmadillo::sample(vars, 1, TRUE, probs)[0];
    }
  }
  return(current);
}

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


//Simulated annealing (not working)
List RcppSA(int n, double tol, int maxiter,
                            IntegerMatrix startmat,
                            NumericVector margprops,
                            List condvars,
                            List condex,
                            List condprobs,
                            bool display_progress=true) {
  
  int i,j, limiter=0,  nvars = condex.size(), startid, accept=0;
  double rmse=1, delta;
  
  IntegerMatrix mat(n, nvars);
  IntegerVector vardims(nvars), vars(nvars), current(nvars), previous(nvars), starts(startmat.nrow()), bin = {0,1};
  NumericVector seq(nvars), currentmargs(margprops.size()), startprobs(starts.size(),1), errorvec(1), jumps(1), P(2), acceptprobs(2);
  List outlist = List::create(Named("Simulation")=0, Named("RMSE")=0, Named("Jumps")=0);
  
  //MinimalProgressBar pb;
  Progress p(n, display_progress);
  
  //Was a start matrix given? If not, generate one at random
  if ( startmat.nrow() > 1 ) { 
    IntegerMatrix startmat(n, nvars);
    for(i=0; i<n; i++) {
      for(j=0; j<nvars; j++){
        vars = condvars[j];
        NumericVector probs(vars.size(),1);
        startmat(i,j) = Rcpp::RcppArmadillo::sample(vars, 1, TRUE, probs)[0];
      }
    }
  }
  
  //Set random start point
  for(i=0; i<starts.size(); i++)
    starts[i] = i;
  
  startid = Rcpp::RcppArmadillo::sample(starts, 1, TRUE, startprobs)[0];
  current = startmat(startid,_);
  
  // Start of algorithm
  i = 0;
  // looping for n samples
  while ( i<n ) {
    
    // Generating a candidate by looping over the conditionals for the variables
    mat(i,_) = GibbsGen(current, condvars, condex, condprobs);
    
    // Getting the RMSE
    currentmargs = check_nprop(mat, condvars, n);
    double thisrmse = 0;
    for(j=0; j<currentmargs.size(); j++)
      thisrmse = thisrmse + (currentmargs[j] - margprops[j])*(currentmargs[j] - margprops[j]);
    
    thisrmse = sqrt(thisrmse);
    delta = (thisrmse - rmse);
    
    //Check if marginals are met
    bool checkprops = true;
    for(j=0; j<currentmargs.size(); j++)
      if( currentmargs[j] > margprops[j] ) checkprops = false;
    
    // If good candidate, keep it.
    if( checkprops == true || thisrmse < rmse ) {
      accept = 1; 
    } else {
      // Else randomly decide if accepting worse candidate, gets smaller as i approaches n
      P[1] = exp(-tol*std::abs(delta)*n*n/(n-i) ) ;
      P[0] = 1 - P[1];
      acceptprobs = clone(P);
      accept = Rcpp::RcppArmadillo::sample(bin, 1, TRUE, acceptprobs)[0];
      Rcout << i << "|" << thisrmse << "|" << rmse << "|"  << delta << "|" << accept << "|" << P << std::endl;
    }
    
    errorvec.push_back(rmse);
    //Check if it improves the fit
    if( accept == 1 ) {
      rmse = thisrmse;
      jumps.push_back(0);
      //Next iter
      i++;
      limiter=0;
      p.increment();
    } else {
      //Checking if no improvement is being made, if not then find new start point
      startid = Rcpp::RcppArmadillo::sample(starts, 1, TRUE, startprobs)[0];
      current = startmat(startid,_);
      jumps.push_back(1);
      limiter++;
    }
    
    if (limiter>maxiter )
      break;
  }
  
  Rcout << "RMSE: " << rmse << std::endl;
  outlist[0] = mat;
  outlist[1] = errorvec;
  outlist[2] = jumps;
  
  return outlist; // Return to R
}
// Metropolis-Hastings (not working)
List RcppMH(int n, int burn,
            IntegerMatrix startmat,
            NumericVector margprops,
            List condvars,
            List condex,
            List condprobs,
            bool display_progress=true) {
  
  int i,j, nvars = condex.size(), startid;
  double a, u, loglik_new, loglik_prev=1, rmse=1, delta;
  
  IntegerMatrix mat(n, nvars);
  IntegerVector vardims(nvars), vars(nvars), current(nvars), previous(nvars), starts(startmat.nrow()), bin = {0,1};
  NumericVector seq(nvars), currentmargs(margprops.size()), startprobs(starts.size(),1), errorvec(1), jumps(1), coin = {1,1};
  List outlist = List::create(Named("Simulation")=0, Named("RMSE")=0, Named("Jumps")=0);
  
  //MinimalProgressBar pb;
  Progress p(n, display_progress);
  
  //Was a start matrix given? If not, generate one at random
  if ( startmat.nrow() > 1 ) { 
    IntegerMatrix startmat(n, nvars);
    for(i=0; i<n; i++) {
      for(j=0; j<nvars; j++){
        vars = condvars[j];
        NumericVector probs(vars.size(),1);
        startmat(i,j) = Rcpp::RcppArmadillo::sample(vars, 1, TRUE, probs)[0];
      }
    }
  }

  //Set random start point
  for(i=0; i<starts.size(); i++) starts[i] = i;
  startid = Rcpp::RcppArmadillo::sample(starts, 1, TRUE, startprobs)[0];
  mat(0,_) = startmat(startid,_);

  // looping for n samples
  for(i=0; i<n; i++) {
    // Generating a candidate by looping over the conditionals for the variables
    previous = mat(i,_);
    mat(i,_) = GibbsGen(previous, condvars, condex, condprobs);
    
    // Getting the RMSE
    currentmargs = check_nprop(mat, condvars, n);
    double thisrmse = 0;
    for(j=0; j<currentmargs.size(); j++) thisrmse = thisrmse + (currentmargs[j] - margprops[j])*(currentmargs[j] - margprops[j]);
    thisrmse = sqrt(thisrmse);
    delta = thisrmse - rmse;
    
    // Calculate current likelihood
    loglik_new = loglikFunc(mat, condvars, margprops);

    a = loglik_new / loglik_prev;
    u = Rcpp::runif(1,0,1)[0];
    Rcout << a << "<=" << u << "|" << loglik_new << "|" << loglik_prev << std::endl;
    if( u <= a ) {
      mat(i,_) = current;
      loglik_prev = loglik_new;
      rmse = thisrmse;
      jumps.push_back(0);
    }
    
    Rcout << a << "|" << delta << std::endl;
    
    // //Checking if no improvement is being made, if not then find new start point
    // startid = Rcpp::RcppArmadillo::sample(starts, 1, TRUE, startprobs)[0];
    // current = startmat(startid,_);
    // jumps.push_back(1);
    errorvec.push_back(rmse);
    
    p.increment();
  }
  
  Rcout << "RMSE: " << rmse << std::endl;
  outlist[0] = mat;
  outlist[1] = errorvec;
  outlist[2] = jumps;
  
  return outlist; // Return to R
}

// Tabu search (not working)
List RcppTabu(int n, double threshold, int maxloop, 
              IntegerMatrix startmat,
              NumericVector margprops,
              List condvars,
              List condex,
              List condprobs,
              bool display_progress=true) {
  
  int h=0,i,j,k, nvars = condex.size(), startid;
  double rmse=1, delta=1;
  
  IntegerMatrix mat(n, nvars);
  IntegerVector vardims(nvars), vars(nvars), current(nvars), previous(nvars), starts(startmat.nrow()), bin = {0,1}, tmp(nvars);
  NumericVector seq(nvars), currentmargs(margprops.size()), startprobs(starts.size(),1), errorvec(1), jumps(1), P(2), acceptprobs(2);
  List outlist = List::create(Named("Simulation")=0, Named("RMSE")=0, Named("Jumps")=0);
  
  //MinimalProgressBar pb;
  Progress p(n, display_progress);
  
  //Was a start matrix given? If not, generate one at random
  if ( startmat.nrow() > 1 ) { 
    IntegerMatrix startmat(n, nvars);
    for(i=0; i<n; i++) {
      for(j=0; j<nvars; j++){
        vars = condvars[j];
        NumericVector probs(vars.size(),1);
        startmat(i,j) = Rcpp::RcppArmadillo::sample(vars, 1, TRUE, probs)[0];
      }
    }
  }
  
  //Set random start point
  for(i=0; i<starts.size(); i++)
    starts[i] = i;
  
  startid = Rcpp::RcppArmadillo::sample(starts, 1, TRUE, startprobs)[0];
  current = startmat(startid,_);
  
  
  // Generating starter matrix to curate
  for(i=0; i<n; i++) {
    mat(i,_) = GibbsGen(current, condvars, condex, condprobs);
    current = mat(i,_);
  }
  
  // Check the RMSE against marginals
  currentmargs = check_nprop(mat, condvars, n);
  double thisrmse = 0;
  for(k=0; k<currentmargs.size(); k++)
    thisrmse = thisrmse + (currentmargs[k] - margprops[k])*(currentmargs[k] - margprops[k]);
  
  rmse = sqrt(thisrmse);
  
  // looping for n samples
  while ( delta > threshold && h < maxloop ) {
    // Goes through each candidate, checking if a replacement is better.
    for(i=0; i<n; i++){
      //Saving last entry
      tmp = mat(i,_);
      previous = clone(tmp);
      mat(i,_) = current;
      
      // Getting the RMSE
      currentmargs = check_nprop(mat, condvars, n);
      double thisrmse = 0;
      for(k=0; k<currentmargs.size(); k++)
        thisrmse = thisrmse + (currentmargs[k] - margprops[k])*(currentmargs[k] - margprops[k]);
      
      thisrmse = sqrt(thisrmse);
      
      //Check if marginals are met
      bool checkprops = true;
      for(k=0; k<currentmargs.size(); k++) 
        if( currentmargs[k] > 1.1*margprops[k]) 
          checkprops = false;
          // Rcout << checkprops << "|" << current << "|" << thisrmse << "|" << rmse << "|"  << delta << std::endl;
      // If good candidate, keep it and continue exploring this locality
      // If not, then reset the last one and draw a new random place
      if( checkprops == true || thisrmse < rmse ) { 
        current = GibbsGen(current, condvars, condex, condprobs);
        jumps.push_back(0);
        rmse = thisrmse;
      } else {
        mat(i,_) = previous;
        startid = Rcpp::RcppArmadillo::sample(starts, 1, TRUE, startprobs)[0];
        current = startmat(startid,_);
        current = GibbsGen(current, condvars, condex, condprobs);
        jumps.push_back(1);
      }
      // Adding RMSE to output
      delta = (thisrmse - rmse);
      errorvec.push_back(rmse);
    }
    
    h++;
    p.increment();
  }
  
  Rcout << "RMSE: " << rmse << std::endl;
  outlist[0] = mat;
  outlist[1] = errorvec;
  outlist[2] = jumps;
  
  return outlist; // Return to R
}

