#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix dpidtheta_calc_cpp(
  NumericMatrix detadtheta,    // (k-1, k-1+nbeta)
  NumericVector dpideta // k-1
)
{
  int    km1=dpideta.length();
  int    ncol1=detadtheta.ncol();
  int    i,j;        // i denotes row, j denotes column,
  
  NumericMatrix dpidtheta(km1, ncol1);
    for(j=0; j<ncol1; ++j)  // j is the column, outer-loop due to column-wise nature
    {
      for(i=0; i<km1; ++i) // i is the row
      {
        dpidtheta[i+j*km1] = detadtheta[i+j*km1]*dpideta[i];
      }
    }
  return(dpidtheta);
}

