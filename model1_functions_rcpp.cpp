#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <cmath>
#include <algorithm>
#include <vector>    
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
//betabinomial function
// [[Rcpp::export]]
double dbbinom(int x, int size, double alpha, double beta, bool log) {
  double ldense = R::lchoose(size, x) +
    R::lbeta((double)x + alpha, (double)size - (double)x + beta) -
    R::lbeta(alpha, beta);
  
  if (log) {
    return(ldense);
  }
  else {
    return(std::exp(ldense));
  }
}

//mu(t)
// [[Rcpp::export]]
double mut_c(double t, double a1, double a2, double a3){
  return a3*(R::dlnorm(t, a1, a2,0));
}


//new tau function
// [[Rcpp::export]]
double tau_c(List paras){
  double et=0;
  double st=1;
  NumericVector a1=paras[0];
  NumericVector a2=paras[1];
  NumericVector a3=paras[2];
  NumericVector theta=paras[3];
  for (int t=0; t<101; t++){
    double at=mut_c(t+1, a1(0), a2(0), a3(0))*((1-mut_c(t+1, a1(0), a2(0), a3(0)))/theta(0)-1);
    double bt=(1-mut_c(t+1, a1(0), a2(0), a3(0)))*((1-mut_c(t+1, a1(0), a2(0), a3(0)))/theta(0)-1);
    et+=st*(1-R::beta(at,7+bt)/R::beta(at,bt))*(t+1);
    st=st*R::beta(at,7+bt)/R::beta(at,bt);
  }
  return et;
}

// [[Rcpp::export]] 
vec lambda_new_cal(vec lambda_new, mat distance, DataFrame inf_list,double b1, double b2, double b3, int j, int n){
 NumericVector place=inf_list[0];
 NumericVector inf=inf_list[2];
 j=j-1;
 for (int t=1; t<n;t++){
   int i=place[t-1]-1;
   if(inf[t-1]>0){
     if(distance(i,j)==1){
       lambda_new[t]=lambda_new[t-1]+b1;
     }else if(distance(i,j)==2){
       lambda_new[t]=lambda_new[t-1]+b1+b2;
     }else if(distance(i,j)==3){
       lambda_new[t]=lambda_new[t-1]+b1+b2+b3;
     }else{
       lambda_new[t]=lambda_new[t-1];
     }
   }else{
       if(distance(i,j)==1){
         lambda_new[t]=lambda_new[t-1]-b1;
       }else if(distance(i,j)==2){
         lambda_new[t]=lambda_new[t-1]-(b1+b2);
       }else if(distance(i,j)==3){
         lambda_new[t]=lambda_new[t-1]-(b1+b2+b3);
       }else{
         lambda_new[t]=lambda_new[t-1];
       }
   }
 }
 return lambda_new;
}

// [[Rcpp::export]] 
double loglikelihood_c(vec lambda_new , vec deltat, int n){
  double loglikelihood=0;
  for(int i=0;i<n;i++){
    loglikelihood-=lambda_new[i]*deltat[i];
  }
  return loglikelihood;
}


//loglikelihood of logposterior
// [[Rcpp::export]]
vec logposterior_c(List paras, mat y, vec onset_c, vec ret, vec first_inf, vec issym, DataFrame data, vec end_place, vec inf_place, vec inf_order, vec deltat, mat lambda_i, int n, int start, int end){ // C version
  //double ret = 0;
  NumericVector a1=paras[0];
  NumericVector a2=paras[1];
  NumericVector a3=paras[2];
  NumericVector theta=paras[3];
  NumericVector m1=paras[4];
  NumericVector block=data["block"];
  NumericVector inf = data["t_inf"];
  NumericVector S1=data["S1"];
  NumericVector S0=data["S0"];
  NumericVector T1=data["T1"];
  NumericVector T0=data["T0"];
  ret[4] += R::dlnorm(m1(0), 1.88, 0.06, 1);
  ret[1] += R::dlnorm(m1(0), 1.88, 0.06, 1);
  for(int i = 0; i < n; i++){
    if(i==first_inf[0]-1||i==first_inf[1]-1||i==first_inf[2]-1||i==first_inf[3]-1||i==first_inf[4]-1) continue;
    if(inf[i]>start&&inf[i]<=end){
      int t= inf_place[i]-1;
      if(lambda_i(i,t-1)<=0){
        ret[0]=ret[1]=-999999;
        return ret;}
      double loglikelihood_i=log(lambda_i(i,t-1));
      for(int j=0;j<t;j++){
        loglikelihood_i=loglikelihood_i-lambda_i(i,j)*deltat[j];
      }
      if(issym[i]>0){
        for(int j=(ceil(inf(i))>0?(ceil(inf(i))-1):0); j<60; j++){
          if (std::isnan(y(issym[i]-1,j))==0)
          {
            double t=(double)j-inf(i);
            double at=mut_c(t+2, a1(0), a2(0), a3(0))*((1-mut_c(t+2, a1(0), a2(0), a3(0)))/theta(0)-1);
            double bt=(1-mut_c(t+2, a1(0), a2(0), a3(0)))*((1-mut_c(t+2, a1(0), a2(0), a3(0)))/theta(0)-1);
            loglikelihood_i += dbbinom(y(issym[i]-1, j), 7, at, bt, 1);
          }
        }
        loglikelihood_i +=R::dlnorm(onset_c(issym[i]-1)-inf(i), log(m1(0))-0.125,0.5, 1);
        ret[4]=ret[4]+loglikelihood_i;
        ret[block[i]+1]=ret[block[i]+1]+loglikelihood_i;
      }
      ret[1]=ret[1]+loglikelihood_i;
      ret[block[i]]=ret[block[i]]+loglikelihood_i;
      
    }else if(S1[i]==0){
      int m=end_place[T1[i]-139];
      for(int j=0; j<m-1; j++){
        ret[2]=ret[2]-lambda_i(i,j)*deltat[j];
      }
      ret[2]=ret[2]-lambda_i(i,m-1)*(T1[i]-139-inf_order[m-1]);
    }else if(S0[i]==0){  // isnan(S1[i])==1&& deleted from condition, not necessary
      int m=end_place[T0[i]-139];
      for(int j=0; j<m-1; j++){
        ret[3]=ret[3]-lambda_i(i,j)*deltat[j];
      }
      ret[3]=ret[3]-lambda_i(i,m-1)*(T0[i]-139-inf_order[m-1]);
    }
  }
  ret[0]=ret[1]+ret[2]+ret[3];
  return ret;
}



//copy lambda_i function (for lambda_matrix)
// [[Rcpp::export]]
mat copy_c(mat lambda_i, mat lambda_small, vec matching, int n, int m){
  for (int i=0;i<n;i++){
    for (int t=0;t<m;t++){
      lambda_i(i,t)=lambda_small(matching[i]-1,t);
    }
  }
  return lambda_i;
}

//incubation likelihood changes
//[[Rcpp::export]]
double incubation_c(List current, double inf_current, double inf_old, mat y, vec onset_c, int k){ // C version
  NumericVector a1=current[0];
  NumericVector a2=current[1];
  NumericVector a3=current[2];
  NumericVector theta=current[3];
  NumericVector m1=current[4];
  double logdelta = R::dlnorm(onset_c(k)-inf_current, log(m1(0))-0.125, 0.5, 1)-R::dlnorm((double)onset_c(k)-inf_old, log(m1(0))-0.125, 0.5, 1);
  for(int j=(ceil(inf_current)>0?(ceil(inf_current-1)):0); j<60; j++){
    if (std::isnan(y(k,j))==0)
    {
      double t=(double)j-inf_current;
      double at=mut_c(t+2, a1(0), a2(0), a3(0))*((1-mut_c(t+2, a1(0), a2(0), a3(0)))/theta(0)-1);
      double bt=(1-mut_c(t+2, a1(0), a2(0), a3(0)))*((1-mut_c(t+2, a1(0), a2(0), a3(0)))/theta(0)-1);
      logdelta += dbbinom(y(k, j), 7, at, bt, 1);
    }
  }
  for(int j=(ceil(inf_old)>0?(ceil(inf_old-1)):0); j<60; j++){
    if (std::isnan(y(k,j))==0)
    {
      double t=(double)j-inf_old;
      double at=mut_c(t+2, a1(0), a2(0), a3(0))*((1-mut_c(t+2, a1(0), a2(0), a3(0)))/theta(0)-1);
      double bt=(1-mut_c(t+2, a1(0), a2(0), a3(0)))*((1-mut_c(t+2, a1(0), a2(0), a3(0)))/theta(0)-1);
      logdelta -= dbbinom(y(k, j), 7, at, bt, 1);
    }
  }
  return logdelta;
}  

