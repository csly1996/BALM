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
/*
VERY IMPORTANT: index in R starts at 1 but in C++ it is 0
USE  armadillo objects instead of Rcpp object to faciliate the function calls
*/

// Enable C++11
// [[Rcpp::plugins("cpp11")]]

/**************************************************
Utility functions Declaration
****************************************************/
// NumericVector which2(LogicalVector x);
// NumericVector naomit(NumericVector x);
// NumericMatrix naomit_m(NumericMatrix x);
// NumericMatrix row_erase (NumericMatrix x, IntegerVector rowID);
// NumericVector concat(NumericVector a, NumericVector b); //c(a,b)
// NumericVector subset_vec(NumericVector a, int fromin, int toin); //c[in:out]
// NumericVector subset_vec2(NumericVector a, NumericVector b); //a[b]
// NumericVector col_vec(NumericMatrix x, int col); // x[, col]
// NumericVector row_vec(NumericMatrix x, int row); //x[row, ]
// NumericVector repNA(int a);
// NumericVector Csample(NumericVector x, int size, bool replace, NumericVector prob = NumericVector::create());
// NumericMatrix cbind_vec(NumericVector x, NumericVector y); // cbind(x,y)
// NumericMatrix rbind_vec(NumericVector x, NumericVector y); // rbind(x,y)
// NumericMatrix col_add(NumericMatrix x, NumericVector y); //cbind(x,y)
// NumericMatrix row_add(NumericMatrix x, NumericVector y); //rbind(x,y)
// NumericMatrix rbind_mat(NumericMatrix x, NumericMatrix y); //rbind(x,y)
// NumericMatrix cbind_mat(NumericMatrix x, NumericMatrix y); //cbind(x,y)
// NumericVector stl_sort(NumericVector x);
// NumericMatrix r_subset_mat(NumericMatrix x, int fromin, int toin); //x[a:b, ]
// NumericMatrix r_subset_mat2(NumericMatrix x, NumericVector y); //x[y, ]
// NumericMatrix c_subset_mat(NumericMatrix x, int fromin, int toin); //x[, a:b]
// NumericMatrix c_subset_mat2(NumericMatrix x, NumericVector y); //x[, y]
// std::string d_to_s(double x);
// NumericVector order_(NumericVector x);



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



/**************************************************************************
utility function end
***************************************************************************/
//test: using it as the main function
// [[Rcpp::export]]
List balm_func1(NumericVector d_obs, NumericMatrix n_obs, int type, bool random=true,
                     int side =1,int pn = 1,
                     NumericVector cpoint = NumericVector::create(0.5), int M = 15000,
                     int burn = 15000,
                     int thin = 35, bool adjust = true,
                     NumericVector w_prior = NumericVector::create(0.1,0.1),
                     NumericVector mu_prior = NumericVector::create(0, 10000),
                     NumericVector tau_prior = NumericVector::create(0,100),
                     int mu_start = 0, double tau_start = 0.2)
{
  int nrow = n_obs.nrow();

  NumericVector n1 = col_vec(n_obs, 1);
  NumericVector n2 = col_vec(n_obs, 2);
  NumericVector cm(nrow);
  cm = (1-3/(4*(n1+n2-2)-1));

  if (adjust == true) {
    d_obs = cm * d_obs;
  }

  int nk = cpoint.size() + 1;
  NumericMatrix index(nk, 2); //nk by 2 matrix with zeros
  std::fill(index.begin(), index.end(), 0);
  NumericMatrix j(d_obs.size(), nk);  //length(d_obs) by nk matrix with NAs
  std::fill(j.begin(), j.end(), NumericVector::get_na() );

  //the equivalently behaving code:
  index(0,1) = 0.5;

  NumericVector p(nrow);
  if (type == 1) {
    if (side == 1) {
      //R::pt()
      // p =  2*(1-R::pt(abs(d_obs/cm)*sqrt(n1/(n1+n2) * n2), df=n1+n2-2))
      for (int i = 0; i < p.size(); i++) {
	       p[i] = 2*(1-R::pt(std::abs(d_obs[i]/cm[i]) * std::sqrt(n1[i]/(n1[i]+n2[i])*n2[i]),
			             n1[i]+n2[i]-2, true, false));
      }
    }
    else if (side == 2) {
      for (int i = 0; i < p.size(); i++) {
	       //p = 1 - R::pt(d_obs/cm*sqrt(n1/(n1+n2)*n2),df=n1+n2-2);
	        p[i] = 1 - R::pt(d_obs[i]/cm[i] * std::sqrt(n1[i]/(n1[i]+n2[i])*n2[i]),
		                       n1[i] + n2[i] -2, true, false);
      }
    }
    else if (side == 3) {
      for (int i = 0; i < p.size(); i++) {
	       p[i] = R::pt(d_obs[i]/cm[i]*std::sqrt(n1[i]/(n1[i]+n2[i])*n2[i]),
		                  n1[i] + n2[i] -2, true, false);
      }
    }
  }
  else if (type == 2) {
    if (pn==1) {
      for (int i = 0; i < p.size(); i++) {
	       p[i] = 1-R::pt(std::abs(d_obs[i]/cm[i])*std::sqrt(n1[i]/(n1[i]+n2[i])*n2[i]),
		                    n1[i]+n2[i]-2, true, false);
      }
    }
    else if (pn==2) {
      for (int i = 0; i < p.size(); i++) {
        p[i] = R::pt(d_obs[i]/cm[i]*std::sqrt(n1[i]/(n1[i]+n2[i])*n2[i]),
		                 n1[i] + n2[i] -2, true, false);
      }
    }
  }


  NumericVector logicvec = which2(p > index(0,0) & p <= index(0,1));
  if(logicvec.size() ==0) {
    for (int i = 0; i < d_obs.size(); i++) {
      j(i,0) = NumericVector::get_na();
    }
  }
  else {
    for (int i = 0; i < logicvec.size(); i++) {
      j(i,0) = logicvec[i];
    }
  }

  //NOTE: matrix access use () instead of []
  //if nk

  //modified
  if (nk > 2) {
    for (int k = 2; k < nk; k++) {
      index(k-1,0) = cpoint[k-2];
      index(k-1,1) = cpoint[k-1];
      NumericVector logicvec2 = which2(p > index(k-1,0) & p <= index(k-1,1));
      if (logicvec2.size() == 0) {
	       for (int i = 0; i < d_obs.size(); i++) {
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

  //NOTE: test right after nk > 2 condition



  index(nk-1, 0) = cpoint[nk-2];
  index(nk-1, 1) = 1;

  NumericVector logicvec3 = which2(p > index(nk-1,0) & p <= index(nk-1,1));
  if (logicvec3.size() == 0) {
    for (int i = 0; i < d_obs.size(); i++) {
      j(i, nk-1) = NumericVector::get_na();
    }
  }
  else {
    for (int i = 0; i < logicvec3.size(); i++) {
      j(i, nk-1) = logicvec3[i];
    }
  }

  int a = d_obs.size();
  NumericMatrix data(a, 3*nk);
  std::fill(data.begin(), data.end(), NumericVector::get_na());

  for (int k = 1; k < nk+1; k++) {//iteration: same as R due to function declaration
    //column vector
    NumericVector colvec = col_vec(j, k);
    NumericVector NaNcolvec = naomit(colvec);
    NumericVector NAvec = repNA(a - NaNcolvec.size());
    NumericVector subvec1 = subset_vec2(n1, NaNcolvec);
    NumericVector subvec2 = subset_vec2(n2, NaNcolvec);
    NumericVector subvec3 = subset_vec2(d_obs, NaNcolvec);
    data(_, k-1) = concat(subvec1, NAvec);
    data(_, nk+k-1) = concat(subvec2, NAvec);
    data(_, 2*nk+k-1) = concat(subvec3, NAvec);
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
   printf(" n_obs matrix dump\n");
   mat_dump(n_obs);
   printf("---------------------------------------\n");
   printf(" vector cm dump\n");
   vector_dump(cm);
   printf("---------------------------------------\n");
   printf("mu_prior dump\n");
   vector_dump(mu_prior);
   printf("---------------------------------------\n");
   printf(" p vector dump\n");
   vector_dump(p);
   */
    //
    // NOTE: EVERYTHING GOOD BEFORE THE FOR LOOP
    //


    //NOTE: no handling of colnames in helper funtions (might need lets see)

  for (int t = 2; t <= M; t++) {
     //int t = 2;
     //printf("iteration: %d\n", t-1);
     NumericVector logi(1000);

     /*
     Behavior of k in R: if a vector of size > 1, set k to be (int):1
     Else, if the vector has only one element, set k to be that number.
     */

     NumericVector temp_k = which2(row_vec(w, t-1) == 1);
     double k;
     //default behavior
     k = temp_k[0];
     //exceptional behavior.
     if (temp_k.size() > 1 || temp_k.size() == 0) { k =  1; }
     m(t-1, k-1) = 0;

     //printf("m should now have second row to be zero\n");
     //vector_dump(r_subset_mat(m, 1, 3));

     //sample function is in RcppArmadillo package
     NumericVector lenvec(d_obs.size());
     for (int i = 0; i < d_obs.size(); i++) {
       lenvec[i] = i+1;
     } // c(1:length(d_obs))

     //use of sample function in RcppArmadillo
     NumericVector sample_ret = Csample(lenvec, 1000-d_obs.size(), true);
     //std::cout<<"n_obs dimension: "<<n_obs.nrow()<<" x "<<n_obs.ncol()<<std::endl;
     NumericMatrix n_mis = r_subset_mat2(n_obs, sample_ret);
     //printf("rbind_mat object specification\n");
     //std::cout<<"n_mis dimension: "<<n_mis.nrow()<<" x "<<n_mis.ncol()<<std::endl;
     //std::cout<<"n_obs dimension: "<<n_obs.nrow()<<" x "<<n_obs.ncol()<<std::endl;
     NumericMatrix n_up = rbind_mat(n_obs, n_mis);

     //NOTE: test dump for first sample call

     /*
     std::cout<<"value of k is: "<<k<<std::endl;
     printf("----------------------------------\n");
     printf("dump for n_mis matrix dimension\n");
     std::cout<<n_mis.nrow()<<" x "<<n_mis.ncol()<<std::endl;
     printf("----------------------------------\n");
     printf("dump for n_up matrix dimension\n");
     std::cout<<n_up.nrow()<<" x "<<n_mis.ncol()<<std::endl;
     */

     // use of rnorm function from Rcpp::rnorm
     //-2 due to the index t used same as in

     //printf("rnorm section\n");
     NumericVector td  = Rcpp::rnorm(1000, mu[t-2], tau[t-2]);
     NumericVector col1 = col_vec(n_up, 1);
     NumericVector col2 = col_vec(n_up, 2);
     NumericVector v = (col1 + col2)/col1/col2 + td*td/(2*(col1 + col2));
     NumericVector d_up(1000);
     // d_up <- rnorm(1000, mean=td, sd=sqrt(v).
     for (int i = 0; i < 1000; i++) {
        d_up[i] = R::rnorm(td[i], std::sqrt(v[i]));
     }


     //NOTE: test dump for first rnorm call
     /*
     printf("first rnorm call\n");
     printf("dump for td vector:\n");
     std::cout<<td.size()<<std::endl;
     vector_dump(subset_vec(td, 1, 20));
     printf("---------------------------------\n");
     printf("dump for v vector\n");
     std::cout<<v.size()<<std::endl;
     vector_dump(subset_vec(v, 1, 20));
     printf("---------------------------------\n");
     printf("dump for d_up vector\n");
     vector_dump(subset_vec(d_up, 1, 20));
     */

     //printf("p_up section and pt\n");
     n1 = col1;
     n2 = col2;
     cm = (1-3/(4*(n1+n2-2)-1));
     NumericVector p_up(d_up.size());
     NumericVector df = n1+n2-2;
     NumericVector temp = (n1/(n1+n2)*n2);
     NumericVector p_dis_abs(d_up.size());
     NumericVector p_dis(d_up.size());

     for (int i = 0; i < d_up.size(); i++) {
       p_dis_abs[i] = R::pt(std::abs(d_up[i]/cm[i]) * std::sqrt(temp[i]), df[i], true, false);
       p_dis[i] = R::pt(d_up[i]/cm[i] * std::sqrt(temp[i]), df[i], true, false);
     }
     if (type == 1) {
       if (side == 1) {
         p_up = 2*(1-p_dis_abs);
       }
       else if (side == 2) {
         p_up = 1 - p_dis;
       }
       else if (side == 3) {
         p_up = p_dis;
       }
     }
     else if (type == 2) {
       if (pn == 1) {
         p_up = 1-p_dis_abs;
       }
       else if (pn == 2) {
         p_up = p_dis;
       }
     }



     //NOTE: dump for the pt part
     //printf(" vector n1 dump\n");
     /*
     std::cout<<n1.size()<<std::endl;
     vector_dump(subset_vec(n1, 1, 20));
     printf(" vector n2 dump\n");
     std::cout<<n2.size()<<std::endl;
     vector_dump(subset_vec(n2, 1, 20));
     printf("vector cm dump\n");
     std::cout<<n2.size()<<std::endl;
     vector_dump(subset_vec(cm, 1, 20));
     printf("vector p_up dump\n");
     std::cout<<p_up.size()<<std::endl;
     vector_dump(p_up);
     */


     //printf("logi seciton\n");
     //logi part
     //k is a vector right now.
     //r_subset_mat(index,k) == index[k,]
     NumericVector whicharr = which2((p_up > index[k-1,0]) & (p_up <= index[k-1,1]));
     for (int i = 0; i < whicharr.size(); i++) {
       logi[whicharr[i]-1] = 1;
     }

     //NOTE: dump for logi
     /*
     printf("----------------------------\n");
     printf(" whicharr vector: \n");
     vector_dump(whicharr);
     printf("----------------------------\n");
     printf("dump for the logi vector\n");
     vector_dump(logi);
     */
     //order(logi, decreasing=T)
     //NOTE: order need further improvement, since it now gives only 1 and 0 cases.

    // printf("order_ section\n");
     //printf("function input check: size of the input: %d\n", logi.size());

     logi = order_(logi);
     /*
     printf("dump for logi after\n");
     printf("------------------------------------\n");
     vector_dump(logi);
     */
     NumericVector logiNaN = naomit(col_vec(j,k));

     /*
     printf("dump for logiNaN\n");
     vector_dump(logiNaN);
     printf("---------------------------------------\n");
     */
     double et;
     et = logi[logiNaN.size()-1];
     if (!et) {et = 0;}
     double m2 = et - d_obs.size();
     if (m2 < 0)  {m2 = 0;}
     //std::cout<<"et is: "<<et<<std::endl;
     //std::cout<<"m2 is: "<<m2<<std::endl;
     double bound = m2 + d_obs.size();
     n_up = r_subset_mat(n_up, 1, bound);
     d_up = subset_vec(d_up, 1, bound);
     p_up = subset_vec(p_up, 1, bound);


     //NOTE: dump before the inner for loop
     //NOTE: order_ function behaves differently in vectors with duplicate
     /*
     printf("before the q loop\n");
     printf("------------------------------\n");
     printf("n_up dump:\n");
     mat_dump(n_up);
     printf("------------------------------\n");
     printf("d_up dump\n");
     vector_dump(d_up);
     printf("------------------------------\n");
     printf("p_up dump\n");
     vector_dump(p_up);
     */

     //NOTE: loop break

    for (int q = 1; q < nk+1; q++) {
       double bound1 = index(q-1, 0);
       double bound2 = index(q-1, 1);
       LogicalVector which_idx2 = ((p_up > bound1) & (p_up <= bound2));
       NumericVector tempwhich = which2(which_idx2);
       m(t-1, q-1) = tempwhich.size() - (naomit(col_vec(j, q))).size();
       if (m(t-1,q-1) < 0) {
         m(t-1,q-1) = 0;
         //sample function use: sample array
         NumericVector id1 = which2((p_up > bound1) & (p_up <= bound2));
         NumericVector samparr((naomit(col_vec(j, q))).size());
         for (int i = 0; i < (naomit(col_vec(j,q)).size()); i++) {
           samparr[i] = i+1;
         }
         //default input as false: need clarification
         NumericVector nid = Csample(samparr, id1.size(), true);
         NumericVector omitset1 = subset_vec2(naomit(col_vec(data, 2*nk +q)),nid);
         for (int i = 0; i < omitset1.size(); i++) {
           d_up[id1[i]-1] = omitset1[i];
         }
         NumericMatrix omitset3 = naomit_m(r_subset_mat2(cbind_vec(col_vec(data, q), col_vec(data, nk+q)),nid));
         for (int i  = 0; i < omitset3.nrow(); i++) {
           n_up(id1[i]-1, _ ) = omitset3( i, _ ); //matrix row substitution
         }
       }
       else  {
         NumericVector id2 = Csample(which2((p_up > bound1) & (p_up <= bound2)), (naomit(col_vec(j, q))).size(), true);
         NumericVector omitset2 = naomit(col_vec(data, 2*nk+q));
         for (int i = 0; i < omitset2.size(); i++) {
           d_up[id2[i]-1] = omitset2[i];
         }
         //n_up to be dealt
         NumericMatrix omitset4 = naomit_m(cbind_vec(col_vec(data,q),col_vec(data,nk+q)));
         for (int i  = 0; i < omitset4.nrow(); i++) {
           n_up(id2[i]-1, _ ) = omitset4( i, _ ); //matrix row substitution
         }
       }


       w(t-1,q-1) = R::rbeta(naomit(col_vec(j, q)).size() + alpha,
                             m(t-1,q-1) + beta);
     }

    //NOTE: dump for rbeta and loop part
    /*
    printf("w matrix dump:\n");
    std::cout<<"w matrix size: "<<w.nrow()<<" x "<<w.ncol()<<std::endl;
    mat_dump(r_subset_mat(w, 1, 20));
    printf("------------------------\n");
    printf("dump for n_up again:\n");
    vector_dump(n_up);
    printf("-------------------------\n");
    printf("dump for d_up again\n");
    vector_dump(d_up);
    printf("--------------------------\n");
    */
     //rescale part
     w(t-1, _) = w(t-1, _) /Rcpp::max(w(t-1, _ ));

    //  printf("after scaling\n");
    //printf("dump for the w matrix with specific rows:\n");
    //mat_dump(r_subset_mat(w, 1, t-1));
     NumericVector dve = d_up;
     NumericMatrix nve = n_up;

     //NOTE: dve nve dump
     /*
     printf("dump for dve:\n");
     vector_dump(dve);
     printf("-----------------\n");
     printf("dump for nve:\n");
     mat_dump(nve);
     */

     NumericVector temp2 = (col_vec(nve, 1) + col_vec(nve, 2))/col_vec(nve,1)/col_vec(nve,2) + dve*dve/2/(col_vec(nve,1) + col_vec(nve,2)) + tau[t-2]*tau[t-2];
     NumericVector numervec = (dve/temp2);
     double numer = std::accumulate(numervec.begin(), numervec.end(), 0.0);
     //std::cout<<"numer: "<<numer<<std::endl;
     NumericVector denovec = (1/temp2);
     double deno = std::accumulate(denovec.begin(), denovec.end(), 0.0)+1/b;
     //std::cout<<"deno: "<<deno<<std::endl;
     mu[t-1] = R::rnorm(numer/deno, std::sqrt(1/deno));

     //NOTE: mu tau dump
     /*
     printf("dump for vector mu\n");
     vector_dump(mu);
     printf("---------------------\n");
     printf("tau vector dump before random condition\n");
     vector_dump(tau);
     */

     if (random == true) {
       tau[t-1] = tau[t-2] + runif(1, -0.1, 0.1)[0];
       if (tau[t-1] > d || tau[t-1] < c) {
         tau[t-1] = tau[t-2];
       }

       NumericVector temp1 = (col_vec(nve, 1) + col_vec(nve, 2))/col_vec(nve,1)/col_vec(nve,2) + dve*dve/2/(col_vec(nve,1) + col_vec(nve,2)) + tau[t-1]*tau[t-1];
       temp2 = (col_vec(nve, 1) + col_vec(nve, 2))/col_vec(nve,1)/col_vec(nve,2) + dve*dve/2/(col_vec(nve,1) + col_vec(nve,2)) + tau[t-2]*tau[t-2];
       double u = runif(1)[0];
       //std::cout<<"value of u : "<<u<<std::endl;
       //NOTE: dump for temp vectors
       /*
       printf("dump temp1: \n");
       vector_dump(temp1);
       printf("--------------------\n");
       printf("dump for temp2: \n");
       vector_dump(temp2);
       printf("-----------------------\n");
       */

       NumericVector logitau_t_vec1(temp1.size());
       for (int i = 0; i < temp1.size(); i++) {
         logitau_t_vec1[i] = (-0.5)* log(temp1[i]);
       }
       NumericVector logitau_t_vec2 = -0.5*(dve-mu[t-1])*(dve-mu[t-1])/temp1;
       NumericVector logitau_t_1_vec1(temp2.size());
       for (int i = 0; i < temp2.size(); i++) {
         logitau_t_1_vec1[i]= -0.5* log(temp2[i]);
       }
       NumericVector logitau_t_1_vec2 = -0.5*(dve-mu[t-1])*(dve-mu[t-1])/temp2;

       /*
       printf("dump for the summing vectors");
       printf("vec1: \n");
       vector_dump(logitau_t_vec1);
       printf("---------------------------\n");
       printf("vec2: \n");
       vector_dump(logitau_t_vec2);
       printf("-----------------------------\n");
       printf("vec3: \n");
       vector_dump(logitau_t_1_vec1);
       printf("--------------------------------\n");
       printf("vec4: \n");
       vector_dump(logitau_t_1_vec2);
       printf("-------------------------------\n");
       */

       double logtau_t = std::accumulate(logitau_t_vec1.begin(), logitau_t_vec1.end(), 0.0) +
                         std::accumulate(logitau_t_vec2.begin(), logitau_t_vec2.end(), 0.0);

       double logtau_t_1 = std::accumulate(logitau_t_1_vec1.begin(), logitau_t_1_vec1.end(), 0.0) +
                           std::accumulate(logitau_t_1_vec2.begin(), logitau_t_1_vec2.end(), 0.0);
       if (std::log(u) > (logtau_t - logtau_t_1)) {
         tau[t-1] = tau[t-2];
       }
      //std::cout<<"logtau_t: "<<logtau_t<<std::endl;
      //std::cout<<"logtau_t_1: "<<logtau_t_1<<std::endl;
     }
     else {
         tau[t-1] = 0;
      }

      //NOTE: dump for the sum part
      /*
      printf("tau vector dump\n");
      NumericVector print_tau = subset_vec(tau, 1, t);
      vector_dump(print_tau);
      printf("---------------------------------------\n");
      */
  }



  //NOTE: PARTS after the main for loop

  NumericMatrix mutau = cbind_vec(mu, tau);
  NumericMatrix results = cbind_mat(w, cbind_mat(m, mutau));

  //seq(burn, M. thin)
  //size: compiler always truncate.
  NumericVector pick((M-burn)/thin + 1);
  pick[0] = burn;
  for (int i = 0; i < pick.size(); i++) {
     pick[i] = burn + i * thin;
  }

  //printf("dump for pick vector\n");
  //vector_dump(pick);
  //printf("------------------------\n");

  //pick is a NumericVector
   w = r_subset_mat2(w, pick);
   m = r_subset_mat2(m, pick);
   mu = subset_vec2(mu, pick);
   tau = subset_vec2(tau, pick);

   //std::cout<<"dump for m and w: "<<std::endl;
   //mat_dump(w);
   //printf("----------------------\n");
   //mat_dump(m);

   std::string note = "";
   for (int i = 0; i < nk; i++) {
     std::string note_tem = " " + names1[i] + " and " + names2[i] + " are corresponding to p value interval: "
                            + d_to_s(index[i,0]) + " to "+ d_to_s(index[i,1]) +" ; ";
     note += note_tem;
   }
   List results1 = List::create(note, Rcpp::Named("w")=w, Rcpp::Named("m")=m, Rcpp::Named("mu")=mu, Rcpp::Named("tau")=tau);
   List results2 = List::create(note, Rcpp::Named("w")=w, Rcpp::Named("m")=m, Rcpp::Named("mu")=mu);
   if (random) {
      return results1;
   }
   else {
       return results2;
   }
}
