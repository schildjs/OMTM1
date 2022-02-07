#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericMatrix hmat_calc_cpp(
  NumericVector Deltavec,  // k-1
  NumericMatrix gammamat   // (k-1, k)
)
{
  int    km1=Deltavec.length();
  int    k=km1+1;      // K - 1
  int    row,col;
  double denom;
  
  NumericMatrix hmat(km1, k);

  // Delta.mat   <- matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
  // hmat.num    <- exp(Delta.mat+gamma.mat)
  // hmat.denom  <- 1 + colSums(hmat.num)
  // hmat        <- sweep(hmat.num,2,hmat.denom,"/")
  for(col=0; col<k; ++col)  // outer-loop due to column-wise nature of calc
  { 
    denom = 1.0; // Denominator for column
    for(row=0; row<km1; ++row) 
    {
      // matrix memory is a column-wise in layout
      // vector sequential in memory.
      // i.e., all matrices appear as as.vector(matrix) in memory
      // location = row + col*nrow 
      hmat[row+col*km1] = exp(Deltavec[row]+gammamat[row+col*km1]);
      denom += hmat[row+col*km1];
    }
    for(row=0; row<km1; ++row) hmat[row+col*km1] /= denom;
  }
 
  return(hmat);
}
