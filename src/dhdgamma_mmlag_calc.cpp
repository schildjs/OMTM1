#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix dhdgamma_mmlag_calc_cpp(
  NumericMatrix hmat1,
  NumericVector mmlag1
)
{
    int    km1=hmat1.nrow();      // K - 1
    int    ncol=km1*km1;
    int    row,col;        // i denotes row, j denotes column,
  //double temp;
    
    NumericMatrix hmat1t(km1, ncol);
    for(col=0; col<ncol; ++col)
    {
        for(row=0; row<km1; ++row)
        {
            hmat1t[row+km1*col] = hmat1[col]*mmlag1[col];
            
        }
    }
//  NumericMatrix dhdgamma_mmlag(km1, ncol);
  //  for(col=0; col<ncol; ++col)
    //{
 //     for(row=0; row<km1; ++row)
   //   {
     //   dhdgamma_mmlag[row+col*km1] = hmat1t[col];
 //   }
   // }
//  NumericMatrix mmlagmat(km1, ncol);
//    for(row=0; row<km1; ++row)
//    {
  //  for(col=0; col<ncol; ++col)
    //{
//        mmlagmat[row+col*km1] = mmlag1[row];
  //  }
    //}
  return(hmat1t);
}
