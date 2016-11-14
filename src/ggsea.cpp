#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix s2n_C(NumericMatrix x, LogicalMatrix y) {
  int n_response = y.ncol();
  int n_gene = x.ncol();
  int n_sample = y.nrow();
  double mean1, mean2, sum1, sum2, var1, var2, m;
  int n1, n2;

  NumericMatrix coef(n_gene, n_response);

  if (n_sample != x.nrow()) {
    return coef;
  }

  for (int i_response=0; i_response<n_response; i_response++) {
    for (int i_gene=0; i_gene<n_gene; i_gene++) {
      sum1 = 0.0;
      sum2 = 0.0;
      var1 = 0.0;
      var2 = 0.0;
      n1 = 0;
      n2 = 0;
      // First pass to compute mean
      for (int i_sample=0; i_sample<n_sample; i_sample++) {
          if (y(i_sample, i_response) == NA_LOGICAL) {
            ;
          } else if (y(i_sample, i_response)) {
            sum1 += x(i_sample, i_gene);
            n1++;
          } else {
            sum2 += x(i_sample, i_gene);
            n2++;
          }
      }
      mean1 = sum1 / n1;
      mean2 = sum1 / n2;
      sum1 = 0.0;
      sum2 = 0.0;
      // Second pass to compute more accurate mean
      for (int i_sample=0; i_sample<n_sample; i_sample++) {
          if (y(i_sample, i_response) == NA_LOGICAL) {
            ;
          } else if (y(i_sample, i_response)) {
            sum1 += x(i_sample, i_gene) - mean1;
          } else {
            sum2 += x(i_sample, i_gene) - mean2;
          }
      }
      mean1 += sum1 / n1;
      mean2 += sum2 / n2;
      // Third pass to compute variance
      sum1 = 0.0;
      sum2 = 0.0;
      for (int i_sample=0; i_sample<n_sample; i_sample++) {
          if (y(i_sample, i_response) == NA_LOGICAL) {
            ;
          } else if (y(i_sample, i_response)) {
            m = x(i_sample, i_gene) - mean1;
            sum1 += m*m;
          } else {
            m = x(i_sample, i_gene) - mean2;
            sum2 += m*m;
          }
      }
      var1 = sum1 / n1;
      var2 = sum2 / n2;
      coef(i_gene, i_response) = (mean1 - mean2) / (sqrt(var1) + sqrt(var2));
    }
  }

  return coef;
}
