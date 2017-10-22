// [[Rcpp::depends(RcppArmadillo)]]
// Usage: <in Rstudio>
// library(Rcpp)
// library(Rcpp11)
// library(RcppArmadillo)
// sourceCpp('~/R_dir/balm.cpp')
//#include <Rcpp.h>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Enable C++11
// [[Rcpp::plugins("cpp11")]]

/*************************************************************************
Utility functions IMPLEMENTATION
**************************************************************************/

//utility function: column vector
//y <- x[, col]
// [[Rcpp::export]]
NumericVector col_vec(NumericMatrix x, int col) {
  arma::mat m1 = as<arma::mat>(x);
  arma::vec v1 = m1.col(col-1);
  return (wrap(v1));
}


//utility funciton: row vector
//y <- x[row, ]
//object returned contains all element in the row,
//but has dimension of a column vector (transpose required if needed)
// [[Rcpp::export]]
NumericVector row_vec(NumericMatrix x, int row) {
  arma::mat m1 = as<arma::mat>(x);
  arma::vec v1 = (m1.row(row-1)).t();
  return (wrap(v1));
}
//utility funcion which2
//INPUT : one logical vector. OUTPUT: numeric vector
// [[Rcpp::export]]
NumericVector which2(LogicalVector x) {
  if (x.size() == 0) {
    printf("WARNING: which2 input vector size is zero\n");
    NumericVector ret(0);
    return ret;
  }
  //printf("size is normal for input vector which2\n");
  int nx = x.size();
  std::vector<double> y;
  y.reserve(nx);
  for (int i = 0; i < nx; i++) {
    if (x[i] == true) y.push_back(i+1);
  }
  NumericVector ret(y.size());
  for (int i = 0; i < y.size(); i++)  {
    ret[i] = y[i];
  }
  return ret;
}

//utility function: removeNA
//INPUT and OUTPUT: NumericVector
// [[Rcpp::export]]
NumericVector naomit(NumericVector x) {
  if (x.size() == 0) {
    printf("INPUT vector has size zero\n");
    return x;
  }
  return (wrap(as<arma::vec>(Rcpp::na_omit(x))));
}


//utility function: row_erase
// [[Rcpp::export]]
NumericMatrix row_erase (NumericMatrix x, IntegerVector rowID) {

  if (x.size() == 0 || rowID.size() == 0) {
  //printf("row_erase: input vector size NULL\n");
    return x;
  }
  rowID = rowID.sort();
  NumericMatrix x2(x.nrow()- rowID.size(), x.ncol());
  int iter = 0;
  int del = 0; // to count deleted elements
  for (int i = 0; i < x.nrow(); i++) {
    if (i != rowID[del]-1) {
      x2(iter, _ ) = x(i, _ );
      iter++;
    } else {
      del++;
    }
  }
  return x2;
}



//utility function: naomit_m
// [[Rcpp::export]]
NumericMatrix naomit_m(NumericMatrix x) {
  if (x.nrow() * x.ncol() == 0) {
    //printf("INPUT matrix has size zero \n");
    return x;
  }
  std::vector<int>idx_vec;
  for (int i = 0; i < x.nrow(); i++) {
    if (is_true(any(is_na(row_vec(x,i+1))))) {
      idx_vec.push_back(i+1);
    }
  }
  IntegerVector idx = wrap(idx_vec);
  NumericMatrix ret = row_erase(x, idx);
  return ret;
}


//utiliy function: concatenation
//combine two vectors together to return a new vector
// [[Rcpp::export]]
// returns a column vector with dimensin n x 1
NumericVector concat(NumericVector a, NumericVector b) {
  arma::vec v1 = as<arma::vec>(a);
  arma::vec v2 = as<arma::vec>(b);
  arma::vec ret = arma::join_cols(v1,v2);
  return (wrap(ret));
}

//utility function: subset contiguous
//a <- c[in:out]
// [[Rcpp::export]]
NumericVector subset_vec(NumericVector a, int fromin, int toin) {
  int from = fromin - 1;
  int to = toin - 1;
  arma::vec v1 = as<arma::vec>(a);
  arma::vec subset = v1.subvec(from, to);
  return (wrap(subset));
}

//utility function: subset
//c <- a[b] where a,b,c are all NumericVectors.
// [[Rcpp::export]]
NumericVector subset_vec2(NumericVector a, NumericVector b) {
  arma::uvec idx = as<arma::uvec>(b) - 1;
  arma::vec v1 = as<arma::vec>(a);
  arma::vec subset = v1.elem(idx);
  return (wrap(subset));
}



//utility function: repNA
// x <- rep(NA, a)
// TEST status: passed
// [[Rcpp::export]]
NumericVector repNA(int a) {
  NumericVector ret(a);
  std::fill(ret.begin(), ret.end(), NumericVector::get_na());
  return ret;
}

//utility function: Csample_sum
//sample function in R. wrapper for armadillo::sample
// exception handling: for empty sample set, return numeric(0) as result.
// [[Rcpp::export]]
NumericVector Csample(NumericVector x, int size, bool replace,
                          NumericVector prob= NumericVector::create())
{
  if (x.size() == 0) {
    printf("WARNING Csample: input vector has size zero\n");
    NumericVector ret(0);
    return ret;
  }
  //printf("input vector has nonzero size\n");
  NumericVector ret = wrap(RcppArmadillo::sample(x, size, replace, prob));
  return ret;
}

//utility function: cbind2, bind columns of two vector to form a matrix
// [[Rcpp::export]]
NumericMatrix cbind_vec(NumericVector x, NumericVector y) {
  arma::vec v1 = as<arma::vec>(x);
  arma::vec v2 = as<arma::vec>(y);
  arma::mat m2(v1.n_elem, 2);
  m2.col(0) = v1;
  m2.col(1) = v2;
  return (wrap(m2));
}

//utility funciton: bind two vector as column matrix to form a matrix
//z <- rbind(x,y)
// [[Rcpp::export]]
NumericMatrix rbind_vec(NumericVector x, NumericVector y) {
  arma::vec v1 = as<arma::vec>(x);
  arma::vec v2 = as<arma::vec>(y);
  arma::mat m2(2, v1.n_elem);
  m2.row(0) = v1.t();
  m2.row(1) = v2.t();
  return (wrap(m2));
}

//utility function: bind a matrix and a vector by column
// [[Rcpp::export]]
NumericMatrix col_add(NumericMatrix x, NumericVector y){
  arma::mat m1 = as<arma::mat>(x);
  arma::vec v1 = as<arma::vec>(y);
  m1.insert_cols(m1.n_cols, v1);
  return (wrap(m1));
}


//utility function: bind a matrix and a vecotr by row
// [[Rcpp::export]]
NumericMatrix row_add(NumericMatrix x, NumericVector y){
  arma::mat m1 = as<arma::mat>(x);
  arma::vec v1 = as<arma::vec>(y);
  m1.insert_rows(m1.n_rows, v1.t());
  return (wrap(m1));
}


//utility function: bind two matrices by row
// [[Rcpp::export]]
NumericMatrix rbind_mat(NumericMatrix x, NumericMatrix y){
  arma::mat m1 = as<arma::mat>(x);
  arma::mat m2 = as<arma::mat>(y);
  arma::mat ret = arma::join_cols(m1,m2);
  return (wrap(ret));
}

//utility function: bind two matrices by column
// [[Rcpp::export]]
NumericMatrix cbind_mat(NumericMatrix x, NumericMatrix y){
  arma::mat m1 = as<arma::mat>(x);
  arma::mat m2 = as<arma::mat>(y);
  arma::mat ret = arma::join_rows(m1,m2);
  return (wrap(ret));
}


//utility function: stlsort in descnding order
// [[Rcpp::export]]
//returns a row vector
NumericVector stl_sort(NumericVector x) {
  if (x.size() == 0) {
    printf("WARNING: stl_sort vector has zero length\n");
    return x;
  }
  NumericVector y =  clone(x);
  std::sort(y.begin(), y.end(), std::greater<double>());
  return y;
}

//utility function: subset a matrix by row, contiguously
// [[Rcpp::export]]
NumericMatrix r_subset_mat(NumericMatrix x, int fromin, int toin) {
  int from = fromin - 1;
  int to = toin - 1;
  arma::mat m1 = as<arma::mat>(x);
  arma::mat s1 = m1.rows(from, to);
  return (wrap(s1));
}

//utility function: row subsetting a matrix object to get a subset matrix
//z <- x[y, ]
// [[Rcpp::export]]
NumericMatrix r_subset_mat2(NumericMatrix x, NumericVector y) {
  arma::mat m1 = as<arma::mat>(x);
  arma::uvec idx = as<arma::uvec>(y) - 1;
  arma::mat s1 = m1.rows(idx);
  return (wrap(s1));
}

//utility funciton: col subsetting a matrix contiguously
// [[Rcpp::export]]
NumericMatrix c_subset_mat(NumericMatrix x, int fromin, int toin) {
  int from = fromin - 1;
  int to = toin - 1;
  arma::mat m1 = as<arma::mat>(x);
  arma::mat s1 = m1.cols(from, to);
  return (wrap(s1));
}

//utility function: col subsetting matrix non-contiguously
// [[Rcpp::export]]
NumericMatrix c_subset_mat2(NumericMatrix x, NumericVector y){
  arma::mat m1 = as<arma::mat>(x);
  arma::uvec idx = as<arma::uvec>(y) -1 ;
  arma::mat s1 = m1.cols(idx);
  return (wrap(s1));
}

//utility: double to string
// [[Rcpp::export]]
std::string d_to_s(double x) {
  std::ostringstream strs;
  strs << x;
  std::string str = strs.str();
  return str;
}

//utility funciotn: order
// [[Rcpp::export]]
//returns a column vector
NumericVector order_(NumericVector x) {
  if (x.size() == 0) {
    printf("order_: input vector size is NULL\n");
    return x;
  }
  int len = x.size();
  //one vector to store 1s, one for 0s and one for index
  std::vector<double> result;
  std::vector<double> ones;
  std::vector<double> zeros;

  for (int i = 0; i < len; i++) {
    if (x[i] == 0) {
        //0 index goes to zeros vector
        zeros.push_back(i+1);
    }
    else if (x[i] == 1) {
      ones.push_back(i+1);
    }
    else {
      printf("ERROR: not ones or zeros\n");
    }
  }
  //after updating a and b, add them to the result vector
  for (int i = 0; i < ones.size(); i++) {
    result.push_back(ones[i]);
  }
  for (int i = 0; i < zeros.size(); i++) {
    result.push_back(zeros[i]);
  }
  //printf("before the wrap function call\n");
  return (wrap(result));

}



//function used for debugging purpose only
void vector_dump(NumericVector x) {
  if (x.size() == 0) {printf("empty vector\n");}
  std::cout<<"vector length: "<<x.size()<<std::endl;
  for (int i = 0; i < x.size(); i++) {
    std::cout<<x[i]<<" ";
  }
  printf("\n");
}

//function used for debugging purpose only
void mat_dump(NumericMatrix x) {
  std::cout<<"Matrix dimensions: "<<x.nrow()<<" x "<<x.ncol()<<std::endl;
  for (int i = 0; i < x.nrow(); i++) {
    for (int j = 0; j < x.ncol(); j++) {
      std::cout<<x(i,j)<<" ";
    }
    printf("\n");
  }
}

//fuction used for debugging purpose only
void logi_dump(LogicalVector x) {
  for (int i = 0;  i < x.size(); i++) {
    std::cout<<x[i]<<" ";
  }
  printf("\n");
}


//function wrapper used for debugging rnorm
// [[Rcpp::export]]
NumericVector Rnorm(int num, double mean, double std) {
  return (Rcpp::rnorm(num, mean, std));
}

//funciton wrapper used for testing rbeta
// [[Rcpp::export]]
double Rbeta(double mean, double std) {
  return (R::rbeta(mean, std));
}

// [[Rcpp::export]]
double Pnorm(double quantile, double std) {
  //default: given value in R function is quantile
  return  (R::pnorm5(quantile, 0, std, 1, 0));

}

/*------------------------------------------------------------------
Main
*-------------------------------------------------------------------*/
//test: using it as the main function
// [[Rcpp::export]]
void balmR(NumericVector r_obs, NumericVector n_obs, int type, bool random=true,
                     int side =1,int pn = 1,
                     NumericVector cpoint = NumericVector::create(0.5), int M = 15000,
                     int burn = 15000,
                     int thin = 35,
                     NumericVector w_prior = NumericVector::create(0.1,0.1),
                     NumericVector mu_prior = NumericVector::create(0, 10000),
                     NumericVector tau_prior = NumericVector::create(0,100),
                     int mu_start = 0, double tau_start = 0.2)
{
  NumericVector n = n_obs;
  NumericVector z_obs(r_obs.size());
  for (int i = 0; i < r_obs.size(); i++) {
    z_obs[i] = 0.5 * log((1+r_obs[i])/(1-r_obs[i]));
  }

  int nk = cpoint.size() + 1;
  NumericMatrix index(nk, 2);
  NumericMatrix j(z_obs.size(), nk);  //length(z_obs) by nk matrix with NAs
  std::fill(j.begin(), j.end(), NumericVector::get_na() );

  index(0,1) = cpoint[0];
  int len = z_obs.size();
  NumericVector p(len);

  if (type == 1) {
      if (side == 1) {
        for(int i =0; i < len; i++) {
          p[i] = 2*(1-Pnorm( std::abs(z_obs[i]), std::sqrt(1/(n[i]-3))));
        }
      }
      else if (side == 2){
        for (int i =0; i < len; i++) {
          p[i] = 1 - Pnorm(z_obs[i], std::sqrt(1/(n[i]-3)));
        }
      }
      else if (side == 3) {
        for (int i = 0; i < len; i++) {
          p[i] = Pnorm(z_obs[i], std::sqrt(1/(n[i]-3)));
        }
      }
  }
  else if (type == 2) {
    if (pn == 1) {
      for (int i = 0; i < len; i++) {
        p[i] = 1 - Pnorm(z_obs[i], std::sqrt(1/(n[i]-3)));
      }
    }
    else if (pn == 2) {
      for (int i = 0; i < len; i++) {
        p[i] = Pnorm(z_obs[i], std::sqrt(1/(n[i]-3)));
      }
    }
  }

  NumericVector logicvec = which2(p > index(0,0) & p <= index(0,1));
  if(logicvec.size() ==0) {
    for (int i = 0; i < z_obs.size(); i++) {
      j(i,0) = NumericVector::get_na();
    }
  }
  else {
    for (int i = 0; i < logicvec.size(); i++) {
      j(i,0) = logicvec[i];
    }
  }
  //modified
  if (nk > 2) {
    for (int k = 2; k < nk; k++) {
      index(k-1,0) = cpoint[k-2];
      index(k-1,1) = cpoint[k-1];
      NumericVector logicvec2 = which2(p > index(k-1,0) & p <= index(k-1,1));
      if (logicvec2.size() == 0) {
         for (int i = 0; i < z_obs.size(); i++) {
         j(i,k-1) = NumericVector::get_na();
         }
      }
      else {
         for (int i =0; i < logicvec2.size(); i++) {
           j(i,k-1) = logicvec2[i];
         }
      }
    }
  } //end of nk>2
  index(nk-1, 0) = cpoint[nk-2];
  index(nk-1, 1) = 1;
  NumericVector logicvec3 = which2(p > index(nk-1,0) & p <= index(nk-1,1));
  if (logicvec3.size() == 0) {
    for (int i = 0; i < z_obs.size(); i++) {
      j(i, nk-1) = NumericVector::get_na();
    }
  }
  else {
    for (int i = 0; i < logicvec3.size(); i++) {
      j(i, nk-1) = logicvec3[i];
    }
  }
  int a = z_obs.size();
  NumericMatrix data(a, 2*nk);
  std::fill(data.begin(), data.end(), NumericVector::get_na());

  for (int k = 1; k < nk+1; k++) {//iteration: same as R due to function declaration
    //column vector
    NumericVector colvec = col_vec(j, k);
    NumericVector NaNcolvec = naomit(colvec);
    NumericVector NAvec = repNA(a - NaNcolvec.size());
    NumericVector subvec1 = subset_vec2(n, NaNcolvec);
    NumericVector subvec2 = subset_vec2(z_obs, NaNcolvec);
    data(_, k-1) = concat(subvec1, NAvec);
    data(_, nk+k-1) = concat(subvec2, NAvec);
  }

  //string vector as names
  NumericMatrix w(M, nk);

  std::vector<std::string> names1(nk);
  for (int i = 0; i < nk; i++) {
   names1[i] = "w " + d_to_s(i+1);
  }
  colnames(w) = wrap(names1);

  NumericMatrix m(M,nk);
  std::vector<std::string> names2(nk);
  for(int i = 0; i < nk; i++) {
    names2[i] = "m " + d_to_s(i+1);
  }
  colnames(m) = wrap(names2);

  double alpha = w_prior[0];
  double beta = w_prior[1];
  a = mu_prior[0];
  double b = mu_prior[1];
  double c = tau_prior[0];
  double d = tau_prior[1];

  for (int i = 0; i < nk; i++) {
    w(0, i) = 0.5;
  }

  NumericVector mu(M);
  mu[0] = mu_start;
  NumericVector tau(M);
  tau[0] = tau_start;

  NumericVector normalizedw = row_vec(w, 1)/max(row_vec(w,1));
  for (int i = 0; i < w.ncol(); i++) {
    w(0, i) = normalizedw[i];
  }


     //NOTE: Testing for the part before for loop
/*
     printf("testing elements before the for loop\n");
     printf("--------------------------------------\n");
     printf("data matrix info: \n");
     mat_dump(data);
     printf("---------------------------------------\n");
     printf("index matrix info: \n");
     mat_dump(index);
     printf("---------------------------------------\n");
     printf(" j matrix info: \n");
     mat_dump(j);
     printf("---------------------------------------\n");
     std::cout<<"m matrix dimensions: "<<m.nrow()<<" x "<<m.ncol()<<std::endl;
     printf("---------------------------------------\n");
     printf(" w matrix info: \n");
     std::cout<<"w matrix dimensions: "<<w.nrow()<<" x "<<w.ncol()<<std::endl;
     vector_dump(w(0, _));
     vector_dump(w(1, _));
     printf("---------------------------------------\n");
     printf(" n vector dump\n");
     vector_dump(n);
     printf("---------------------------------------\n");
     printf(" z_obs dump\n");
     vector_dump(z_obs);
     printf("---------------------------------------\n");
     printf("mu_prior dump\n");
     vector_dump(mu_prior);
     printf("---------------------------------------\n");
     printf(" p vector dump\n");
     vector_dump(p);
*/


}
