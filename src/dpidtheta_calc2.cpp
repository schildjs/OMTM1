#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix dpidtheta_calc2_cpp(
  NumericMatrix diagmtx,    // (k-1, k-1)
  NumericVector dpideta, // k-1
  NumericVector xit, // ncol in non-ppo part of design matrix
  NumericMatrix ppoMat    // (k-1, ncol) ppo part of design matrix    
)
{
  int    km1=dpideta.length();
  int    nbeta=xit.length();
  int    nbetappo = ppoMat.ncol();
  int    ncol1=nbeta+km1;
  int    ncol2=ncol1+nbetappo;
  int    i,j;        // i denotes row, j denotes column,
  
//ncol1=nbetappo+nbeta+km1;
  NumericMatrix dpidtheta(km1, ncol2);
    for(j=0; j<km1; ++j)  // j is the column, outer-loop due to column-wise nature
    {
      for(i=0; i<km1; ++i) // i is the row
      {
        dpidtheta[i+j*km1] = diagmtx[i+j*km1]*dpideta[i];
      }
    }
    for(j=km1; j<ncol1; ++j)  // j is the column, outer-loop due to column-wise nature
    {
      for(i=0; i<km1; ++i) // i is the row
      {
          dpidtheta[i+j*km1] = xit[(j-km1)]*dpideta[i];
      }
    }
    for(j=ncol1; j<ncol2; ++j)  // j is the column, outer-loop due to column-wise nature
    {
        for(i=0; i<km1; ++i) // i is the row
        {
            dpidtheta[i+j*km1] = ppoMat[i+(j-ncol1)*km1]*dpideta[i];
        }
    }
  return(dpidtheta);
}

