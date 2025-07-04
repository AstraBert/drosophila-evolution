#include <Rcpp.h>

using namespace Rcpp;

double dist1 (NumericVector x, NumericVector y){
  int n = y.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += (x(i)-y(i))*(x(i)-y(i));
  }
  total = sqrt(total);
  return total;
}

// [[Rcpp::export]]
NumericVector calDIST (NumericVector y, NumericMatrix x){
  int n = x.nrow();
  NumericVector out(n);

  for (int i = 0 ; i < n ; i++){
      NumericVector v = x.row(i);
      double d = dist1(y, v);
      out(i) = d;
	}
  return (out);
}
