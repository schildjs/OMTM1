#include <Rcpp.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

using namespace Rcpp;

// [[Rcpp::export()]]
NumericVector Delta_calc_cpp(
  NumericVector mm,       // k
  NumericVector mmlag,    // k
  NumericMatrix gammamat, // (k-1, k)
  double        tol,
  int           maxit,
  int           trace
)
{
  int    k=mm.length();
  int    km1=k-1;      // K - 1
  int    i,j,l;        // i denotes row, j denotes column, 
  int    itr;          // Track iterations
  double temp;
  double temp2;
  double maxslope;     // Finds maximum slope
  
  NumericMatrix hmat(km1, k);
  NumericMatrix dfdDelta(km1, km1);
  NumericVector del(km1);
  NumericVector fDelta(km1);
  NumericVector Deltavec(km1);
  IntegerVector ipiv(km1);
  NumericMatrix work(km1, km1);
  
  // Double check assumptions about size
  if(mm.length() != mmlag.length())
  {
    ::Rf_error("mm[%d] != mm.lag[%d]\n", mm.length(), mmlag.length());
  }
  if(gammamat.nrow() != km1)
  {
    ::Rf_error("gamma[%d, %d] rows != k-1\n", gammamat.nrow(), gammamat.ncol());
  }
  if(gammamat.ncol() != k)
  {
    ::Rf_error("gamma[%d, %d] cols != k\n", gammamat.nrow(), gammamat.ncol());
  }
  
  itr=0;
  maxslope=tol+0.1; // Make sure first iteration happens
                
  while(maxslope > tol)
  {
    if(++itr > maxit) ::Rf_error("Maximum Iterations Exceeded");
    
     if(trace > 0) Rprintf("Iteration %d, ", itr);
    
    // Delta.mat   <- matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
    // hmat.num    <- exp(Delta.mat+gamma.mat)
    // hmat.denom  <- 1+ colSums(hmat.num)
    // hmat        <- sweep(hmat.num,2,hmat.denom,"/")
    for(j=0; j<k; ++j)  // j is the column, outer-loop due to column-wise nature
    { 
      temp = 1.0; // Denominator for column
      for(i=0; i<km1; ++i) // i is the row
      {
        // matrix memory is a column-wise in layout
        // vector sequential in memory.
        // i.e., all matrices appear as as.vector(matrix) in memory
        // location = row + col*nrow 
        hmat[i+j*km1] = exp(Deltavec[i]+gammamat[i+j*km1]);
        temp += hmat[i+j*km1];
      }
      for(i=0; i<km1; ++i) hmat[i+j*km1] /= temp;
    }
 
    // fDelta <- hmat %*% mm.lag - mm
    for(i=0; i<km1; ++i) fDelta[i] = -mm[i]; // BLAS will modify this vector
    i    = 1;   // Each entry in vector is one apart, BLAS allows arbitrary spacing
    temp = 1.0; // Operations are scaled by 1.0
    dgemv_("N", &km1, &k, &temp, &hmat[0], &km1, &mmlag[0], &i, &temp, &fDelta[0], &i);
    // Ref: http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gadd421a107a488d524859b4a64c1901a9.html#gadd421a107a488d524859b4a64c1901a9
  
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
        dfdDelta[i+i*(i+1)/2] += hmat[i+j*km1] * (1-hmat[i+j*km1])*mmlag[j];
      }
    }

    // ## upper triangle for df.dDelta
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
          dfdDelta[i+j*(j+1)/2] -= hmat[i+l*km1]*hmat[j+l*km1]*mmlag[l];
        }
      }
    }
  
    //del <- solve(df.dDelta) %*% fDelta
    // Solve inverse of real symmetric indefinite matrix in packed storage
    // https://www.math.utah.edu/software/lapack/lapack-d/dsptrf.html
    dsptrf_("U", &km1, &dfdDelta[0], &ipiv[0], &i);
    if(i != 0) ::Rf_error("Bunch-Kaufman factorization failed: info %d", i);
    // https://www.math.utah.edu/software/lapack/lapack-d/dsptri.html
    dsptri_("U", &km1, &dfdDelta[0], &ipiv[0], &work[0], &i);
    if(i != 0) ::Rf_error("Inversion failed, DSPTRI Info %d", i);
    
    temp=1.0;  // 1.0 * fDelta
    temp2=0.0; // Ignore contents of del
    i=1;       // Elements are next to each other (i.e. no comb like gaps)
    // This is the  "%*% fDelta" piece with a triangular packed matrix
    // Result is left in del
    //http://www.netlib.org/lapack/explore-html/d7/d15/group__double__blas__level2_gab746575c4f7dd4eec72e8110d42cefe9.html#gab746575c4f7dd4eec72e8110d42cefe9
    dspmv_("U", &km1, &temp, &dfdDelta[0], &fDelta[0], &i, &temp2, &del[0], &i);
    
    // Modify Deltavec elementwise and find maximum absolute slope.
    maxslope = 0.0;
    if(trace > 0) Rprintf("  del: ");

    for(i=0; i<km1; ++i)
    {
      // Delta.vec  <- Delta.vec - del
      Deltavec[i] -= del[i];

      if(fabs(del[i]) > maxslope) maxslope = fabs(del[i]); 
      
      if(trace > 0) Rprintf("%le ", del[i]);
    }
    if(trace > 0) Rprintf("\n");

  }
  
  return(Deltavec);
}
