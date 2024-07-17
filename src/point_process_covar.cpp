#include <Rcpp.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector newfun (NumericVector c){
  return c;
}

// [[Rcpp::export]]
SEXP test(SEXP h){
  return CDR(h);
}

class point{
  
};



typedef  std::vector<point> pointprocess;
  


