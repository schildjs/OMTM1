#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix dfdDelta_calc_cpp(
  NumericVector mmlag,    // k
  NumericMatrix hmat // (k-1, k)
)
{
  int    k=mmlag.length();
  int    km1=k-1;      // K - 1
    int    i,j,l;        // i denotes row, j denotes column,
  //double temp;
  
  NumericMatrix dfdDelta(km1, km1);
  
    // df.dDelta   <- matrix(0, K1, K1)
    // memset is a very fast shortcut that needs the number
    // of bytes to overwrite with a byte (0). 
    // IEEE754 all zero bytes is a zero for a double.
    // sizeof returns the number of bytes of a double,
    memset(&dfdDelta[0], 0, km1*km1*sizeof(double));
    
    // for (l in 1:K1)
    //   df.dDelta[l,l] <- sum(hmat[l,]*(1-hmat[l,])*mm.lag)
    for(i=0; i<km1; ++i)
    {
      for(j=0; j<k; ++j)
      {
        // Use packed "U" col major storage, 
        // Bytes are stored using triangle number offsets
        // https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/lapack-routines/matrix-storage-schemes-for-lapack-routines.html
        //dfdDelta[i+i*(i+1)/2] += hmat[i+j*km1] * (1-hmat[i+j*km1])*mmlag[j];
        dfdDelta[k*i] += hmat[i+j*km1] * (1-hmat[i+j*km1])*mmlag[j];
          
      }
    }

    // ## off-diagonals for df.dDelta
    // R code only shows for upper triangle
    // for (l in 1:(K1-1))
    // {
    //   for (m in (l+1):K1)
    //   {
    //     df.dDelta[l,m] <- -sum(hmat[l,]*hmat[m,]*mm.lag)
    //   }
    // }
    for(i=0; i<(km1-1); ++i)   // for (l in 1:(K1-1))
    {
      for(j=i+1; j<km1; ++j)   //   for (m in (l+1):K1)
      {
        for(l=0; l<k; ++l)  // Elements to sum
        {
          dfdDelta[i+j*km1] -= hmat[i+l*km1]*hmat[j+l*km1]*mmlag[l];
          dfdDelta[j+i*km1] -= hmat[i+l*km1]*hmat[j+l*km1]*mmlag[l];
        }
      }
    }
  return(dfdDelta);
}
