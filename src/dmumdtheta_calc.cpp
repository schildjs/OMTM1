#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix dmumdtheta_calc_cpp(
  NumericMatrix dpidtheta   // (k-1, k)
)
{
    int    ncol=dpidtheta.ncol(); // K-1+nbeta
    int    km1=dpidtheta.nrow(); // K-1
    int    k=km1+1; // K
    int    row,col;
  
    NumericMatrix dmumdtheta(k, ncol);
    // dpidtheta is a km1 by ncol matrix where pis are cumulative probabilities
    // and we need to calculate derivative of the state probability (mu) with respect
    // to parameter theta.  Notice that each dpidtheta is basically pi*(1-pi).
    //dpi.k.dtheta      <- rbind(dpidtheta, 0)
    // dpi.k1.dtheta     <- rbind(0,dpidtheta)
    // dmum.k.dtheta.tmp <- dpi.k.dtheta - dpi.k1.dtheta
  for(col=0; col<ncol; ++col)  // outer-loop due to column-wise nature of calc
  {
      dmumdtheta[0+col*k] = dpidtheta[0+col*(k-1)];
      dmumdtheta[km1+col*k] = 0.0-dpidtheta[km1-1+col*(k-1)];
      for(row=1; row<km1; ++row)
    {
        dmumdtheta[row+col*k] = dpidtheta[row+col*(k-1)]-dpidtheta[row-1+col*(k-1)];
    }
  }
  return(dmumdtheta);
}
