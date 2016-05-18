#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
#include <math.h> 

// SDistanceMat must be a symmetric matrix
// [[Rcpp::export]]
Rcpp::NumericMatrix  TDistanceMat_cpp(Rcpp::NumericMatrix& SDistanceMat,
  float& A,
  float& b){
    
    int nr=SDistanceMat.nrow(),nc=SDistanceMat.ncol();
    Rcpp::NumericMatrix tmat(nr,nc);
    tmat.attr("dimnames")=SDistanceMat.attr("dimnames");
    for(int i=0;i<nr;i++){
      for(int j=i;j<nc;j++){
        float val=A*exp(-b*SDistanceMat(i,j));
        tmat(i,j)=val;
        tmat(j,i)=val;
      }
    }
    
    return(tmat);
}

// [[Rcpp::export]]
float ICod_onepair_cpp(Rcpp::NumericMatrix& SDis,
  Rcpp::NumericMatrix& TDis,
  float& C){
    
  float numerator=0,denominator=0;
  for(int i=0; i<SDis.nrow();i++){
    for(int j=0;j<SDis.ncol();j++){
      if(SDis(i,j)<=C){
        numerator +=TDis(i,j);
      }
      denominator += TDis(i,j);
    }
  }
  return(numerator/denominator);
}







