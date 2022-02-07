#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix dhd_gamma_calc(
    NumericMatrix hmat,    // IN (k-1, k) 
    NumericVector mm_lag,  // IN (k)
    NumericMatrix tmp_mat) // IN (k-1, (k-1)^2)
{
  int km1 = hmat.nrow();
  int i, j, k;
  int loc;
  
  NumericMatrix dhd_gamma(tmp_mat.nrow(), tmp_mat.ncol());
  
  for(i=0; i<km1; ++i)
  {
    for(j=0; j<km1; ++j)
    {
      for(k=0; k<km1; ++k)
      {
        // dhdgamma.mm.lag.2[, (i-1)+(j-1)*Km1+1] <- -hmat[j, i]*mm.lag[i]
        loc = k+(i+j*km1)*km1;
        dhd_gamma[loc] = 
          -1.0          *
          hmat[j+i*km1] *
          mm_lag[i]     *
          (tmp_mat[loc] - hmat[k+ i*km1]);
      }
    }
  }
  

  return(dhd_gamma);
}
