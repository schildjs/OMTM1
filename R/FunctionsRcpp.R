# Rcpp interfers with inline and crashes R
# Doing this last seems to work
library(Rcpp)
# For Local compile
Sys.setenv(PKG_LIBS="-llapack")
# For ACCRE compile
#Sys.setenv(PKG_LIBS="-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl")
#Sys.setenv(PKG_CXXFLAGS="-I${MKLROOT}/include")

ThePath <- "src/"
#######################################################################################################
sourceCpp(paste(ThePath,"Delta_calc.cpp", sep=""))

Delta_calc <- function(mm, mm.lag, gamma.mat, tol=1e-4, maxit=10000, trace=0)
    Delta_calc_cpp(mm, mm.lag, gamma.mat, tol, maxit, trace) 

#######################################################################################################
sourceCpp(paste(ThePath,"Delta_calcWith0TxProbs.cpp", sep=""))

Delta_calcWith0TxProbs <- function(mm, mm.lag, gamma.mat, CalcTxMtx, tol=1e-4, maxit=10000, trace=0)
    Delta_calcWith0TxProbs_cpp(mm, mm.lag, gamma.mat, CalcTxMtx, tol, maxit, trace) 

#######################################################################################################
sourceCpp(paste(ThePath,"dfdDelta_calc.cpp", sep=""))

dfdDelta_calc <- function(mm.lag, hmat)
    dfdDelta_calc_cpp(mm.lag, hmat) 
#######################################################################################################
sourceCpp(paste(ThePath,"dhd_gamma_calc.cpp", sep=""))

dhdgamma_mmlag_calc <- function(hmat, mm.lag, tmp.mat)
    dhd_gamma_calc(hmat, mm.lag, tmp.mat) 

#######################################################################################################
sourceCpp(paste(ThePath,"dmumdtheta_calc.cpp", sep=""))

dmumdtheta_calc <- function(dpi.dtheta)
    dmumdtheta_calc_cpp(dpi.dtheta) 

#######################################################################################################
sourceCpp(paste(ThePath,"dpidtheta_calc.cpp", sep=""))

dpidtheta_calc <- function(deta.k.dtheta, dpi.deta)
    dpidtheta_calc_cpp(deta.k.dtheta, dpi.deta) 

#######################################################################################################
sourceCpp(paste(ThePath,"dpidtheta_calc1.cpp", sep=""))

dpidtheta_calc1 <- function(diagmtx, dpi.deta, xit)
    dpidtheta_calc1_cpp(diagmtx, dpi.deta, xit) 

#######################################################################################################
sourceCpp(paste(ThePath,"dpidtheta_calc2.cpp", sep=""))

dpidtheta_calc2 <- function(diagmtx, dpi.deta, xit, ppo.Imat)
    dpidtheta_calc2_cpp(diagmtx, dpi.deta, xit, ppo.Imat) 

#######################################################################################################
sourceCpp(paste(ThePath,"hmat_calc.cpp", sep=""))

hmat_calc <- function(Delta.vec, gamma.mat)
    hmat_calc_cpp(Delta.vec, gamma.mat) 

