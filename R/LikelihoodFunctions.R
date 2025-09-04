## 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
####################################################################################
### This is the original R version of the Delta solver that was converted to cpp for speed
### Solving for the (K-1)-vector of Delta_{itk}... at time t
findDeltait = function(mm,            # row of marginal (mean) probabilities of each response category k=1...K
                        mm.lag,       # lagged value of mm
                        gamma.mat,    # Transition log odds ratios matrix ((a,b)^th element corresponds to being in state a on day j conditional on being in state b on day j-1
                        K){     # K
    K1 = K-1
    mm          = mm[-K]
    Delta.vec   = rep(0, K1)
    del         = rep(1, K1)
    while (max(abs(del))>1e-4){

        Delta.mat   = matrix(rep(Delta.vec,each=K), ncol=K, byrow=TRUE)
        hmat.num    = exp(Delta.mat+gamma.mat)
        hmat.denom  = 1+ colSums(hmat.num)
        hmat        = sweep(hmat.num,2,hmat.denom,"/")
        fDelta      = hmat %*% mm.lag - mm
        df.dDelta   = matrix(0, K1, K1)
        ## diagonal for df.dDelta
        for (l in 1:K1){
            df.dDelta[l,l] = sum(hmat[l,]*(1-hmat[l,])*mm.lag)
        }
        ## upper triangle for df.dDelta
        for (l in 1:(K1-1)){
            for (m in (l+1):K1){
                df.dDelta[l,m] = -sum(hmat[l,]*hmat[m,]*mm.lag)
            }}
        ## lower triangle for df.dDelta
        df.dDelta[lower.tri(df.dDelta)] =  t(df.dDelta)[lower.tri(df.dDelta)]
        #print(df.dDelta)
        del = solve(df.dDelta) %*% fDelta
        Delta.vec  = Delta.vec - del
    }
    Delta.vec
}
####################################################################################
expit = plogis

## Transition probability, joint probability and overall state probability
#' Calculate transition probabilities, joint probabilities and overall state probability
#' Calculate transition probabilities, joint probabilities and overall state probability
#' @param Y response variable which should be an integer value
#' @param id id variables that is the same length as Y
#' @param digs number of rounding digits
#' @return transition probabilities, joint probabilities and overall state probability
#' @export
#'
Calc.TransProbs = function(Y, id, digs){

    L = length(Y)
    Y1 = c(NA,Y[-L])
    dup = duplicated(id)
    Y = Y[dup]
    Ylag = Y1[dup]
    L = length(Y)

    tab = table(Y,Ylag)
    joint = tab/L
    trans = sweep(tab, 2, colSums(tab), "/")

    out = list(Tx.Mtx.Col2Row=round(trans,digs), Joint.Mtx=round(joint, digs), Marg.Prob=round(table(Y)/L, digs))
    out
}

## This function is involved in the likelihood function, but does get used
CreateTranIndMtx = function(Y, id){

    K    = length(unique(Y))
    L    = length(Y)
    Y1   = c(NA,Y[-L])
    dup  = duplicated(id)
    Y    = Y[dup]
    Ylag = Y1[dup]
    L    = length(Y)

    AnyTrans = 1*(table(Y,Ylag)>0)[-K,]
    AnyTrans
}

####################################################################################
####################################################################################
#### We modified the likelihoods and data generating models to permit the dependence model
#### to be non-stationary or to depend on covariates.  In limited testing this seems to
#### work, but it is not fully tested.
####################################################################################
####################################################################################
####################################################################################

GenDatOMTM1.ppo = function(id,
                            XMat, #### design matrix that should NOT include any intercepts
                            alpha, #### K-1 vectors of intercepts
                            beta, #### vector of proportional odds coefficient
                            UMat=NULL, #### Design matrix for the transitions, e.g., G1 * 1 + G2 * tx if G1 is the log relative risk ratio matrix for tx=0 and G1+G2 for tx=1
                            gamma.mat.list, ## This should be a list of (K-1) by (K-1) matrices for the dependence model.  It can be a list of length 1 for a standard transition model
                            ppo.k=NULL,     ## Dichotomizations that correspond to non-PO parts of the model e.g., ppo.k=c(2,3,2,3)
                            XMat.ppo=NULL,  ## Design matrix for the non-PO parts of the model, e.g. XMat.ppo=XMat[,c(2,2,3,3)]
                            beta.ppo=NULL){ ## Coefficients corresponding to deviations from the PO assumption

    ppo.n = length(ppo.k)

    x          =  XMat
    lp         = x %*% beta
    gamma.mat0 = cbind(gamma.mat.list, 0)
    K1         = length(alpha)
    K          = K1+1
    start      = 1

    if (is.null(UMat)) UMat = as.matrix(rep(1, length(lp)), ncol=1)
    u      = UMat
    ncol.u = ncol(u)

    l.gam           = length(gamma.mat.list)
    gamma.mat0.list = list(l.gam)
    for (mm in 1:l.gam) gamma.mat0.list[[mm]] = cbind(gamma.mat.list[[mm]], 0)
    #print(gamma.mat0.list)
    yval = rep(NA, length(lp))
    for (i in unique(id)){
        idi = id[id==i]
        lpi = lp[id==i]
        mi  = length(idi)
        xi  = x[id==i,]
        ui  = u[id==i,,drop=FALSE]
        ## marginal probability matrix and lagged values
        ## rows correspond to time
        ## columns correspond to outcome category
        lpi.mat = matrix(NA, mi, K1)

        ## po or ppo linear predictor
        if (is.null(ppo.k)){
            for (j in 1:mi){
                for( k in 1:K1){
                    lpi.mat[j,k] = lpi[j]+alpha[k] }}
        }else{
            xi.ppo = XMat.ppo[id==i,]
            for (j in 1:mi){
                xit.ppo = xi.ppo[j,]
                for( k in 1:K1){
                    tmp.1        = matrix(0, 1, ppo.n)
                    blah         = which(ppo.k==k)
                    tmp.1[blah]  = xit.ppo[blah]
                    lpi.mat[j,k] = lpi[j]+alpha[k] + tmp.1 %*% beta.ppo }}
        }

        cprobi.mat = cbind(0, expit(lpi.mat), 1) ## cumulative prob matrix
        probi.mat  = matrix(NA, mi, K)           ## multinomial probabilities
        for (k in 2:(K+1)){
            probi.mat[,(k-1)] = cprobi.mat[,k]-cprobi.mat[,(k-1)]
        }

        ### Lagged value of the marginal probability matrix
        ### The first row here has no impact on data generation.  We just need to put something there
        probi.matlag = rbind(probi.mat[1,], probi.mat[-mi,])

        ## Calculate Deltait across all timepoints.
        Deltait = NULL
        for (j in 1:mi){
            uit     = ui[j,,drop=FALSE]
            if (l.gam==1){
                gamma.mat0 = gamma.mat0.list[[1]]
            }else{
                gamma.mat0 = do.call("+", Map("*", gamma.mat0.list, uit) )
            }
            #print(gamma.mat0)
            Deltait = rbind(Deltait, c(Delta_calc(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0)))}

        # use marginal probs to simulate Yi1
        yi = yilag = c(rmultinom(1,1,probi.mat[1,]))

        ## simulate Yi2... Yimi
        for (j in 2:mi){
            tmp      = exp(Deltait[j,] + c(gamma.mat0 %*% yilag))
            tmp1     = 1+sum(tmp)
            prob.y.c = c(tmp,1)/tmp1
            yilag    = c(rmultinom(1,1,prob.y.c))
            yi       = rbind(yi, yilag)
        }
        yival                          = rep(1e6, nrow(yi))
        for(j in 1:nrow(yi)){ yival[j] = which(yi[j,]==1) }
        yval[c(start:(start+mi-1))]    = yival
        start                          = start+mi
    }
    yval
}


####################################################################################
#### Create the list of data that is used in the likelihood fitting
####################################################################################
CreateSubjectData = function(id, yval,x,x.ppo,u,ZMat,YMat, wt, cprob.m,
                              prob.m,one.min.cprob.m,phi.k, g.phi.k,
                              dphi.k.deta.k,dphi.k.deta.k1,dpi.k.deta.k, K,K1){

    id.tmp               = split(id,id)
    yival.tmp            = split(yval,id)
    xi.tmp               = split(x,id)
    ui.tmp               = split(u,id)
    ZiMat.tmp            = split(ZMat,id)
    YiMat.tmp            = split(YMat,id)
    wt.tmp               = split(wt,id)
    cprobi.m.tmp         = split(cprob.m,id)
    probi.m.tmp          = split(prob.m,id)
    one.min.cprobi.m.tmp = split(one.min.cprob.m,id)
    phii.k.tmp           = split(phi.k,id)
    g.phii.k.tmp         = split(g.phi.k,id)
    dphii.k.deta.k.tmp   = split(dphi.k.deta.k,id)
    dphii.k.deta.k1.tmp  = split(dphi.k.deta.k1,id)
    dpii.k.deta.k.tmp    = split(dpi.k.deta.k,id)
    ncolx                = ncol(x)
    ncolu                = ncol(u)

    #################################################
    if (is.null(x.ppo)){
        xi.ppo.tmp = split(x,id)
        ncolx.ppo = ncolx
    }else{
        xi.ppo.tmp = split(x.ppo,id)
        ncolx.ppo  = ncol(x.ppo)
    }
    #################################################

    subjectData = vector('list', length=length(unique(id)))
    subjectData = list()
    uid = as.character(unique(id))
    for(j in seq(along=uid)){
        i = uid[j]
        subjectData[[j]] = list(idi      = as.character(unique(id.tmp[[i]])),
                                yival    = yival.tmp[[i]],
                                mi       = length(yival.tmp[[i]]),
                                xi       = matrix(xi.tmp[[i]],ncol=ncolx),
                                xi.ppo   = matrix(xi.ppo.tmp[[i]],ncol=ncolx.ppo),
                                ui       = matrix(ui.tmp[[i]],ncol=ncolu),
                                ZiMat    = matrix(ZiMat.tmp[[i]],ncol=K),
                                YiMat    = matrix(YiMat.tmp[[i]],ncol=K),
                                wt       = unique(wt.tmp[[i]]),
                                cprobi.m = matrix(cprobi.m.tmp[[i]],ncol=K),
                                probi.m  = matrix(probi.m.tmp[[i]],ncol=K),
                                one.min.cprobi.m = matrix(one.min.cprobi.m.tmp[[i]],ncol=K),
                                phii.k           = matrix(phii.k.tmp[[i]],ncol=K1),
                                g.phii.k         = matrix(g.phii.k.tmp[[i]],ncol=K1),
                                dphii.k.deta.k   = matrix(dphii.k.deta.k.tmp[[i]],ncol=K1),
                                dphii.k.deta.k1  = matrix(dphii.k.deta.k1.tmp[[i]],ncol=K1),
                                dpii.k.deta.k    = matrix(dpii.k.deta.k.tmp[[i]],ncol=K1))
    }
    names(subjectData) = uid
    subjectData
}

logLikeCalc = function(params, ## vector of parameter values.  Gamma is row major.  When modified by covariates, input the Gamma one matrix at a time, e.g., params = c(alpha, beta, c(t(gamma.mat.list[[1]]), t(gamma.mat.list[[2]])))
                        yval, ## vector of ordered response values.  In this function, these must take on integer values from 1 to K.
                        x, ## design matrix
                        x.ppo=NULL, ## design matrix for non-proportional odds
                        ppo.k=NULL,  ## outcome dichotomizations corresponding to non-PO in the columns of x.ppo
                        u=NULL, ## design matrix for response dependence model (modifiers of the relationship between Y_i(t_{ij}) and Y_i(t_{ij-1})).  The number of rows should equal the length of yval
                        wt, ## vector of sampling weights
                        id, ## id vector, same length as yval
                        ProfileCol=NA, ## column value to fix.  This is particularly important when we have absorbing states and need transition matrix probabilities to be 1 or 0
                        ref.muc=NA, ## The reference values for the lagged response value in the transition matrix.  I believe it is important for this reference state to have observed transition for it to every other state
                        TransIndMtx=NA, ## used in case there are 0 probability transition.  It should be a K-1 by K matrix with 1s in all cells except for those that correspond to 0 probability transtions.  Do not use this right now.  Take care of 0 transition probabilities by fixing them with gamma values and using ProfileCol.
                        UseGrad=TRUE, ## should we output the gradient as an attribute to be used in the nlm fit
                        CheeseCalc=FALSE){ ## should we calculate the robust standard errors

    #
    if (is.null(u)) u = matrix(1, length(yval), 1)
    States            = sort(unique(yval))
    K                 = length(States)
    K1                = K-1
    K1sq              = K1^2
    alpha.len         = K1
    ncol.x            = ncol(x)
    ncol.x.ppo        = ifelse(is.null(x.ppo), 0, ncol(x.ppo))
    ncol.u            = ncol(u)
    npar              = length(params)
    uid               = unique(id)
    N                 = length(uid)
    Ntot              = length(yval)

    alpha.e   = params[1:K1]
    beta.e    = params[(K1+1):(K1+ncol.x+ncol.x.ppo)]
    gamma.vec = params[(K1+ncol.x+ncol.x.ppo+1):npar]
    #gamma.e   = matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
    gamma.len         = length(gamma.vec)
    gamma.len.test    = if (gamma.len != K1sq*ncol.u){stop("Something went wrong with the dimension of the parameters.")}

    # gamma is a list of matrices.
    gamma.mtx         = list(ncol.u)
    for (g in 1:ncol.u){ tmp.indx       = 1:K1sq + (g-1)*K1sq
    tmp.mtx        = matrix(gamma.vec[tmp.indx], K1, K1, byrow=TRUE)
    gamma.mtx[[g]] = cbind(tmp.mtx, 0)
    }

    alpha.len = K1
    beta.len  = length(beta.e)
    gamma.len = length(gamma.vec)

    ## matrices and indices for partial PO model
    ## numbering of utility matrices tmp* is not good and needs to be fixed
    if (!is.null(ppo.k)){
        x.ppo0 = x.ppo*0
        n.ppo = ncol(x.ppo)
        tmp3a = matrix(0, K1, n.ppo)
    }

    ## This is not being used but still appears in a c++ function, and I have not yet removed it.
    if (!is.matrix(TransIndMtx)){ TransIndMtx = matrix(1, K1, K)}

    ## Z and Y matrices
    ZMat = YMat = matrix(0, Ntot, K)
    for (k in 1:K){ YMat[,k] = ifelse(yval==k, 1, 0)
    ZMat[,k] = ifelse(yval<=k, 1, 0)
    }

    ## Marginal cumulative probs and state probs
    ## The else part regards partial PO
    cprob.m = prob.m = matrix(NA, Ntot, K)
    for (k in 1:K1){
        if (is.null(ppo.k)){cprob.m[,k]   = expit(alpha.e[k] + x %*% beta.e)
        }else{  x.ppo.k       = x.ppo0
        aaa           = which(ppo.k==k)
        x.ppo.k[,aaa] = x.ppo[,aaa]
        cprob.m[,k]   = expit(alpha.e[k] + cbind(x,x.ppo.k) %*% beta.e)
        }
    }
    cprob.m[,K] = 1
    prob.m[,1]  = cprob.m[,1]
    for (k in 2:K){
        prob.m[,k] = cprob.m[,(k)]- cprob.m[,(k-1)]
    }
    one.min.cprob.m = 1-cprob.m
    one.min.prob.m  = 1-prob.m

    ## Phi and g(Phi) (we only use the first observation for each subject)
    phi.k = matrix(NA, Ntot, K1)
    for (k in 1:K1){
        phi.k[,k] = log(cprob.m[,k]/(prob.m[,(k+1)]))
    }
    g.phi.k     = log(1+exp(phi.k))
    exp.g.phi.k = exp(g.phi.k)

    ## dPhi/deta (we only use the first observation for each subject)
    ## dpi/deta
    dphi.k.deta.k = dphi.k.deta.k1 = dpi.k.deta.k = matrix(NA, Ntot, K1)
    for (k in 1:K1){
        dphi.k.deta.k[,k]  =  exp.g.phi.k[,k]*one.min.cprob.m[,k]
        dphi.k.deta.k1[,k] =  -exp.g.phi.k[,k]*one.min.cprob.m[,(k+1)]
        dpi.k.deta.k[,k]   =  cprob.m[,k]*one.min.cprob.m[,k]
    }

    ## Some useful matrices to be used in gradient calculations
    ## Numbering of these needs to be improved
    tmp1  = diag(1,K1)
    tmp2  = rbind(tmp1[-1,], rep(0,K1))
    tmp4  = matrix(0, K1, ncol.u*(K1^2))

    ## To be used in the gradient calculation
    tmp.mat = matrix(0, K1, K1^2)
    for (ppp in 1:K1){
        tmp.mat[ppp,(c(1:K1)+(ppp-1)*K1)] = 1
    }

    li1 = li2 = rep(0,N)
    Grad1Vec = Grad2Vec = rep(0,npar)
    Grad1Mat = Grad2Mat = matrix(0,N, npar)

    subjectData = CreateSubjectData(id, yval,x,x.ppo, u, ZMat,YMat, wt, cprob.m,prob.m,one.min.cprob.m,phi.k, g.phi.k,dphi.k.deta.k,dphi.k.deta.k1,dpi.k.deta.k,K,K1)
    blah        = lapply(subjectData, LogLikeiCalc, gamma.mtx=gamma.mtx, tmp1=tmp1, tmp2=tmp2, tmp3a=tmp3a, tmp4=tmp4, tmp.mat=tmp.mat, ref.muc=ref.muc,
                         K=K, K1=K1, n.ppo=n.ppo, ppo.k=ppo.k, npar=npar, alpha.len=alpha.len, beta.len=beta.len,
                         gamma.len=gamma.len, TransIndMtx=TransIndMtx)


    li1 = lapply(blah, function(x) x[['logLi1']])
    loglike1 = -1*sum(c(unlist(li1)))
    li2 = lapply(blah, function(x) x[['logLi2']])
    loglike2 = -1*sum(c(unlist(li2)))
    gradi1 = lapply(blah, function(x) x[['Grad1iVec']])
    gradient1 = -1*Reduce('+', gradi1)
    gradi2 = lapply(blah, function(x) x[['Grad2iVec']])
    gradient2 = -1*Reduce('+', gradi2)

    gradient1[ProfileCol] = 0
    gradient2[ProfileCol] = 0

    loglikelihood  = loglike1+loglike2
    if (UseGrad==TRUE) attr(loglikelihood,"gradient") = c(gradient1+gradient2)
    if (CheeseCalc==TRUE){
        gradi = mapply(function(a,b) a+b, gradi1, gradi2, SIMPLIFY=FALSE)
        gradi.outer = mapply(outer, gradi, gradi, SIMPLIFY=FALSE)
        attr(loglikelihood,"cheese") = Reduce("+", gradi.outer)
    }
    loglikelihood
}

LogLikeiCalc = function(subjectData, gamma.mtx, tmp1, tmp2, tmp3a, tmp4, tmp.mat, ref.muc, K, K1,
                         n.ppo, ppo.k, npar, alpha.len, beta.len, gamma.len, TransIndMtx){

    yival            = subjectData[["yival"]]
    xi               = subjectData[["xi"]]
    xi.ppo           = subjectData[["xi.ppo"]]
    ui               = subjectData[["ui"]]
    mi               = subjectData[["mi"]]
    ZiMat            = subjectData[["ZiMat"]]
    YiMat            = subjectData[["YiMat"]]
    wti              = subjectData[["wt"]]
    cprobi.m         = subjectData[["cprobi.m"]]
    probi.m          = subjectData[["probi.m"]]
    one.min.cprobi.m = subjectData[["one.min.cprobi.m"]]
    phii.k           = subjectData[["phii.k"]]
    g.phii.k         = subjectData[["g.phii.k"]]
    dphii.k.deta.k   = subjectData[["dphii.k.deta.k"]]
    dphii.k.deta.k1  = subjectData[["dphii.k.deta.k1"]]
    dpii.k.deta.k    = subjectData[["dpii.k.deta.k"]]
    K1sq             = K1^2
    ncol.u           = ncol(ui)

    if (mi==1){
        xi       = matrix(xi,nrow=1)
        xi.ppo   = matrix(xi.ppo,nrow=1)
        ZiMat    = matrix(ZiMat,nrow=1)
        YiMat    = matrix(YiMat,nrow=1)
        cprobi.m = matrix(cprobi.m, nrow=1)
        probi.m  = matrix(probi.m, nrow=1)
        phii.k          = matrix(phii.k, nrow=1)
        g.phii.k        = matrix(g.phii.k, nrow=1)
        dphii.k.deta.k  = matrix(dphii.k.deta.k, nrow=1)
        dphii.k.deta.k1 = matrix(dphii.k.deta.k1, nrow=1)
        dpii.k.deta.k   = matrix(dpii.k.deta.k, nrow=1)
    }

    logLi1    = logLi2 = 0
    logLij2   = rep(0, mi)
    Grad1iVec = Grad2iVec = rep(0, npar)

    ## LogLi1 and dLogLi1/d(alpha,beta)
    ## LogLi1 only comes from the marginal portion of the model (first observation)
    for (k in 1:K1){ #print(cprobi.m[1,]);print(probi.m[1,])
        logLi1 = logLi1 +  ZiMat[1,k]*phii.k[1,k] - ZiMat[1,(k+1)]*g.phii.k[1,k]
    }

    ## LogLi1/dtheta
    ## the else term involves the partial proportional odds part of the model
    xi1     = xi[1,]
    xi1.ppo = xi.ppo[1,]
    tmp3    = matrix(rep(xi1,K1), nrow=K1, byrow=TRUE)

    if (is.null(ppo.k)){
        deta.k.dtheta  = cbind(tmp1,tmp3,tmp4)
        deta.k1.dtheta = cbind(tmp2,tmp3,tmp4)
    }else{
        tmp3b = tmp3a
        for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xi1.ppo[rrr]}
        deta.k.dtheta         = cbind(tmp1,tmp3,tmp3b, tmp4)
        deta.k1.dtheta        = cbind(tmp2,tmp3,tmp3b, tmp4)
    }

    for (k in 1:K1){
        Grad1iVec = Grad1iVec+(ZiMat[1,k]-(ZiMat[1,(k+1)]*cprobi.m[1,k]/cprobi.m[1,(k+1)]))*(dphii.k.deta.k[1,k]*deta.k.dtheta[k,] + dphii.k.deta.k1[1,k]*deta.k1.dtheta[k,])
    }

    if (mi>1){
        ## logLi2 calculated from observations 2 to mi
        ## Need to save this as a matrix to allow for the binary case
        YiMat2 = matrix(YiMat[,-K], ncol = K1)

        ############ Altered for new ref.muc (reference state in conditional model)
        if (!is.na(ref.muc)){
            probi.m   = cbind(probi.m[,-ref.muc], probi.m[,ref.muc])
            YiMat.tmp = cbind(YiMat[,-ref.muc], YiMat[,ref.muc])
            YiMat2    = YiMat.tmp[,-K]
            yival2    = ifelse(yival<ref.muc, yival,
                               ifelse(yival==ref.muc, K, yival-1))
            yival     = yival2
        }

        ### Note that I am keeping some of the R code to help myself remember more easily what each of these steps are doing
        ### without having to go to the C++ code.
        for (j in 2:mi){

            mm           = probi.m[j,]
            mm.lag       = probi.m[(j-1),]
            Yit          = YiMat2[j,]
            Yit1         = YiMat2[(j-1),]
            xit          = xi[j,]
            xit1         = xi[(j-1),]
            xit.ppo      = xi.ppo[j,]
            xit1.ppo     = xi.ppo[(j-1),]
            uit          = ui[j,]
            dpi.deta     = dpii.k.deta.k[j,]
            dpi.deta.lag = dpii.k.deta.k[(j-1),]

            ## u*gamma in matrix form.  This is the matrix that has the element that pre-mutiplies y[t-1] in the dependence model
            ## Multiply the gamma matrices by corresponding elements of uit and then sum them.
            gamma.u = Reduce("+", Map("*", gamma.mtx, uit))


            #############################################################################################################
            #Deltaij2    = Delta_calc(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx)
            Deltaij2    = Delta_calcWith0TxProbs(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.u, CalcTxMtx=TransIndMtx)
            #############################################################################################################

            ## we can see that ordered categories must take on integer values from 1 to K.
            #lpij.c     = Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
            lpij.c      = Deltaij2 + gamma.u[,yival[(j-1)]] #%*% YiMat2[(j-1),]
            probij.cK   = 1/(1+sum(exp(lpij.c)))

            logLij2[j]  = Yit %*% lpij.c + log(probij.cK)
            muit.c      = exp(lpij.c)*probij.cK

            #############################################################################################################
            #hmat = hmat_calc(Deltaij2, gamma.mtx)
            hmat = hmat_calc(Deltaij2, gamma.u)*TransIndMtx

            #############################################################################################################

            # df.dDelta   = matrix(0, K1, K1)
            # ## Left side of system of equations that solve for dDelta/dtheta
            # for (l in 1:K1){df.dDelta[l,l] = sum(hmat[l,]*(1-hmat[l,])*mm.lag)}
            # ## upper triangle for df.dDelta
            # for (l in 1:(K1-1)){
            #     for (m in (l+1):K1){
            #         df.dDelta[l,m] = -sum(hmat[l,]*hmat[m,]*mm.lag)
            #     }}
            # ## lower triangle for
            # df.dDelta[lower.tri(df.dDelta)] =  t(df.dDelta)[lower.tri(df.dDelta)]
            df.dDelta = dfdDelta_calc(mm.lag=mm.lag, hmat=hmat)

            ## right side of system of equations that solve for dDelta/dtheta (only alpha and beta)
            ## the else part is for partial PO
            if (is.null(ppo.k)){
                dpidtheta.lag = dpidtheta_calc1(tmp1, dpi.deta.lag, xit1)
                dpidtheta     = dpidtheta_calc1(tmp1, dpi.deta, xit)
            }else{
                tmp3b = tmp3a
                tmp3c = tmp3a
                for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xit1.ppo[rrr]
                tmp3c[ppo.k[rrr],rrr] = xit.ppo[rrr]
                }
                dpidtheta.lag = dpidtheta_calc2(tmp1, dpi.deta.lag, xit1, tmp3b)
                dpidtheta     = dpidtheta_calc2(tmp1, dpi.deta,     xit,  tmp3c)
            }
            ########################################################################
            dmumdtheta.lag = dmumdtheta_calc(dpidtheta.lag)
            ## for new ref.muc #############################
            if (!is.na(ref.muc)) dmumdtheta.lag = rbind(dmumdtheta.lag[-ref.muc,], dmumdtheta.lag[ref.muc,])
            h.dmumdtheta.lag     = hmat %*% dmumdtheta.lag

            ########################################################################
            dmumdtheta = dmumdtheta_calc(dpidtheta)
            ## for new ref.muc #############################
            if (is.na(ref.muc)){dmumdtheta  = dmumdtheta[-K,]
            }else{              dmumdtheta  = dmumdtheta[-ref.muc,]}

            rhs.alpha.beta  = dmumdtheta - h.dmumdtheta.lag

            ########################################################################
            #tmp.mat.tmp = list(ncol.u)
            #for (gg in 1:ncol.u) tmp.mat.tmp[[gg]] = tmp.mat
            #tmp.mat.1 = do.call("cbind", tmp.mat.tmp)
            #### NOTE: This is multiplied by -1 in the cpp function.  So it's actually -1*dhdgamma.mm.lag
            dhdgamma.mm.lag = dhdgamma_mmlag_calc(hmat, mm.lag, tmp.mat)
            #dhdgamma.mm.lag = dhdgamma_mmlag_calc(hmat, mm.lag, tmp.mat.1)
            dhdgamma.mm.lag.list = list(ncol.u)
            #for (gg in 1:ncol.u) dhdgamma.mm.lag.list[[gg]] = dhdgamma.mm.lag[,(1:K1sq + (gg-1)*K1sq)]
            for (gg in 1:ncol.u) dhdgamma.mm.lag.list[[gg]] = dhdgamma.mm.lag
            rhs.gamma = do.call("cbind", Map("*",dhdgamma.mm.lag.list, uit))

            rhs                 = cbind(rhs.alpha.beta, rhs.gamma)
            ddelta.dtheta       = solve(df.dDelta, rhs)
            y.min.muc.ylag      = c(rep(Yit-muit.c, each=K1) * rep(Yit1, K1))
            y.min.muc.ylag.list = list(ncol.u)
            for (gg in 1:ncol.u) y.min.muc.ylag.list[[gg]] = y.min.muc.ylag*uit[[gg]]
            y.min.muc.dugammadgamma.ylag = do.call("c", y.min.muc.ylag.list)

            GradiAdd = c((Yit-muit.c) %*% ddelta.dtheta + c(rep(0, alpha.len+beta.len), y.min.muc.dugammadgamma.ylag))
            Grad2iVec = Grad2iVec + GradiAdd
        }
    }

    out = list(logLi1=logLi1*wti,
               logLi2=sum(logLij2)*wti,
               Grad1iVec=Grad1iVec*wti,
               Grad2iVec=Grad2iVec*wti)
    out
}

CalcVarCov = function(MOD, epsilon, yval, x, wt, id, ProfileCol=NA, ref.muc=NA, TransIndMtx=NA, x.ppo=NULL, ppo.k=NULL, u=NULL){
    npar        = length(MOD$estimate)
    eps.mtx     = diag(rep(epsilon, npar))
    grad.at.max = MOD$gradient
    ObsInfo.tmp = ObsInfo = matrix(NA, npar, npar)
    if (is.null(u)) u = matrix(1, length(yval), 1)

    ## Observed Information
    for (j in 1:npar){
        temp            = logLikeCalc(MOD$estimate+eps.mtx[j,], yval=yval, x=x, wt=wt,
                                      x.ppo=x.ppo, ppo.k=ppo.k, u=u, id=id,
                                      ProfileCol = ProfileCol, ref.muc=ref.muc,
                                      TransIndMtx=TransIndMtx, UseGrad=TRUE)
        ObsInfo.tmp[j,] = (attr(temp,"gradient")-grad.at.max)/(epsilon)
    }
    for (m in 1:npar){ for (n in 1:npar){ ObsInfo[m,n] =  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}

    ## remove rows and columns of information matrix with no effective information.  This is usually caused by absorbing states
    rmv.noInfo = which(diag(ObsInfo)<.01)

    ## Model based covariance
    if (length(rmv.noInfo)>0){
        mod.cov = solve(ObsInfo[-rmv.noInfo,-rmv.noInfo])
    }else{                     mod.cov = solve(ObsInfo)}

    ## Robust, sandwich covariance
    cheese  = attr(logLikeCalc(params=MOD$estimate, yval=yval, x=x, wt=wt,
                               x.ppo=x.ppo, ppo.k=ppo.k, u=u, id=id, ProfileCol=ProfileCol,
                               ref.muc=ref.muc, TransIndMtx=TransIndMtx, CheeseCalc=TRUE), "cheese")
    if (length(rmv.noInfo)>0){
        cheese  = cheese[-rmv.noInfo,-rmv.noInfo]
    }
    rob.cov = mod.cov %*% cheese %*% mod.cov

    out = list(mod.cov=mod.cov, rob.cov=rob.cov, cheese=cheese, low.info=rmv.noInfo)
    out

}

cluster.summary = function( id, x, fun ){
    xlist = split( x, id )
    nj = unlist( lapply( xlist, length ) )
    xj = unlist( lapply( xlist, fun) )
    xsummary = rep( xj, nj )
    xsummary}

#############################
#############################
## Older versions of data generation and fitting functions that I do not want to throw away yet
#############################
#############################


#
# ####################################################################################
# GenDatOMTM1 = function(id,
#                         XMat, #### design matrix that should NOT include any intercepts
#                         alpha,
#                         beta,
#                         gamma.mat){
#
#     lp         = XMat %*% beta
#     gamma.mat0 = cbind(gamma.mat, 0)
#     K1         = length(alpha)
#     K          = K1+1
#     start      = 1
#
#     ## linear predictor without the intercepts
#     #y = yval = NULL
#     yval = rep(NA, length(lp))
#
#     for (i in unique(id)){
#         idi = id[id==i]
#         lpi = lp[id==i]
#         mi  = length(idi)
#
#         ## marginal mean / probability matrix and lagged values
#         ## rows correspond to time
#         ## columns correspond to outcome category
#         lpi.mat = matrix(NA, mi, K1)
#         for (j in 1:mi){ for( k in 1:K1){ lpi.mat[j,k] = lpi[j]+alpha[k] }}
#
#         cprobi.mat = cbind(0, expit(lpi.mat), 1) ## cumuulative prob matrix
#         probi.mat  = matrix(NA, mi, K) ## multinomial probabilities
#         for (k in 2:(K+1)){
#             probi.mat[,(k-1)] = cprobi.mat[,k]-cprobi.mat[,(k-1)]
#         }
#
#         ### This value has any impact based on the way I am generating the data.
#         probi.matlag = rbind(probi.mat[1,], probi.mat[-mi,])
#
#         ## Calculate Deltait across all timepoints.
#         Deltait = NULL
#         for (j in 1:mi){ Deltait = rbind(Deltait, c(Delta_calc(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0)))}
#
#         # use marginal probs to simulate Yi1
#         yi = yilag = c(rmultinom(1,1,probi.mat[1,]))
#
#         ## sim Yi2... Yimi
#         for (j in 2:mi){
#             tmp      = exp(Deltait[j,] + c(gamma.mat0 %*% yilag))
#             tmp1     = 1+sum(tmp)
#             prob.y.c = c(tmp,1)/tmp1
#             yilag    = c(rmultinom(1,1,prob.y.c))
#             yi       = rbind(yi, yilag)
#         }
#         yival = rep(10, nrow(yi))
#         for(j in 1:nrow(yi)){ yival[j] = which(yi[j,]==1) }
#         yval[c(start:(start+mi-1))] = yival
#         start = start+mi
#
#     }
#     yval
# }
#
# GenDatOMTM1.ppo = function(id,
#                         XMat, #### design matrix that should NOT include any intercepts
#                         alpha, #### K-1 vectors of intercepts
#                         beta, #### vector of proportional odds coefficient
#                         gamma.mat, ## K-1 by K-1 matrix of transition log relative risks ()
#                         ppo.k=NULL, ## Dichotomizations that correspond to non-PO parts of the model
#                         XMat.ppo=NULL, ## Design matrix for the non-PO parts of the model
#                         beta.ppo=NULL){ ## Coefficients corresponding to deviations from the PO assumption
#
#     ppo.n = length(ppo.k)
#
#     x          = XMat
#     lp         = x %*% beta
#     gamma.mat0 = cbind(gamma.mat, 0)
#     K1         = length(alpha)
#     K          = K1+1
#     start      = 1
#
#     yval = rep(NA, length(lp))
#     for (i in unique(id)){
#         idi = id[id==i]
#         lpi = lp[id==i]
#         mi  = length(idi)
#         xi  = x[id==i,]
#
#         ## marginal probability matrix and lagged values
#         ## rows correspond to time
#         ## columns correspond to outcome category
#         lpi.mat = matrix(NA, mi, K1)
#
#         ## po or ppo linear predictor
#         if (is.null(ppo.k)){
#             for (j in 1:mi){
#                 for( k in 1:K1){
#                     lpi.mat[j,k] = lpi[j]+alpha[k] }}
#         }else{
#             xi.ppo = XMat.ppo[id==i,]
#             for (j in 1:mi){
#                 xit.ppo = xi.ppo[j,]
#                 for( k in 1:K1){
#                     tmp.1        = matrix(0, 1, ppo.n)
#                     blah         = which(ppo.k==k)
#                     tmp.1[blah]  = xit.ppo[blah]
#                     lpi.mat[j,k] = lpi[j]+alpha[k] + tmp.1 %*% beta.ppo }}
#         }
#
#         cprobi.mat = cbind(0, expit(lpi.mat), 1) ## cumulative prob matrix
#         probi.mat  = matrix(NA, mi, K)           ## multinomial probabilities
#         for (k in 2:(K+1)){
#             probi.mat[,(k-1)] = cprobi.mat[,k]-cprobi.mat[,(k-1)]
#         }
#
#         ### The first row here has no impact on data generation.  We just need to put something there
#         probi.matlag = rbind(probi.mat[1,], probi.mat[-mi,])
#
#         ## Calculate Deltait across all timepoints.
#         Deltait = NULL
#         for (j in 1:mi){ Deltait = rbind(Deltait, c(Delta_calc(mm=probi.mat[j,], mm.lag=probi.matlag[j,], gamma.mat=gamma.mat0)))}
#
#         # use marginal probs to simulate Yi1
#         yi = yilag = c(rmultinom(1,1,probi.mat[1,]))
#
#         ## sim Yi2... Yimi
#         for (j in 2:mi){
#             tmp      = exp(Deltait[j,] + c(gamma.mat0 %*% yilag))
#             tmp1     = 1+sum(tmp)
#             prob.y.c = c(tmp,1)/tmp1
#             yilag    = c(rmultinom(1,1,prob.y.c))
#             yi       = rbind(yi, yilag)
#         }
#         yival                          = rep(1e6, nrow(yi))
#         for(j in 1:nrow(yi)){ yival[j] = which(yi[j,]==1) }
#         yval[c(start:(start+mi-1))]    = yival
#         start                          = start+mi
#     }
#     yval
# }

#
#
#
#
# ####################################################################################
# #### Create the list of data that is used in the likelihood fitting
# ####################################################################################
# CreateSubjectData = function(id, yval,x,x.ppo,ZMat,YMat,cprob.m,
#                              prob.m,one.min.cprob.m,phi.k, g.phi.k,
#                              dphi.k.deta.k,dphi.k.deta.k1,dpi.k.deta.k, K,K1){
#
#     id.tmp               = split(id,id)
#     yival.tmp            = split(yval,id)
#     xi.tmp               = split(x,id)
#     ZiMat.tmp            = split(ZMat,id)
#     YiMat.tmp            = split(YMat,id)
#     cprobi.m.tmp         = split(cprob.m,id)
#     probi.m.tmp          = split(prob.m,id)
#     one.min.cprobi.m.tmp = split(one.min.cprob.m,id)
#     phii.k.tmp           = split(phi.k,id)
#     g.phii.k.tmp         = split(g.phi.k,id)
#     dphii.k.deta.k.tmp   = split(dphi.k.deta.k,id)
#     dphii.k.deta.k1.tmp  = split(dphi.k.deta.k1,id)
#     dpii.k.deta.k.tmp    = split(dpi.k.deta.k,id)
#     ncolx                = ncol(x)
#
#     #################################################
#     if (is.null(x.ppo)){
#         xi.ppo.tmp = split(x,id)
#         ncolx.ppo = ncolx
#     }else{
#         xi.ppo.tmp = split(x.ppo,id)
#         ncolx.ppo  = ncol(x.ppo)
#     }
#     #################################################
#
#     subjectData = vector('list', length=length(unique(id)))
#     subjectData = list()
#     uid = as.character(unique(id))
#     for(j in seq(along=uid)){
#         i = uid[j]
#         subjectData[[j]] = list(idi = as.character(unique(id.tmp[[i]])),
#                                  yival = yival.tmp[[i]],
#                                  mi = length(yival.tmp[[i]]),
#                                  xi   = matrix(xi.tmp[[i]],ncol=ncolx),
#                                  xi.ppo = matrix(xi.ppo.tmp[[i]],ncol=ncolx.ppo),
#                                  ZiMat = matrix(ZiMat.tmp[[i]],ncol=K),
#                                  YiMat   = matrix(YiMat.tmp[[i]],ncol=K),
#                                  cprobi.m = matrix(cprobi.m.tmp[[i]],ncol=K),
#                                  probi.m   = matrix(probi.m.tmp[[i]],ncol=K),
#                                  one.min.cprobi.m = matrix(one.min.cprobi.m.tmp[[i]],ncol=K),
#                                  phii.k   = matrix(phii.k.tmp[[i]],ncol=K1),
#                                  g.phii.k = matrix(g.phii.k.tmp[[i]],ncol=K1),
#                                  dphii.k.deta.k   = matrix(dphii.k.deta.k.tmp[[i]],ncol=K1),
#                                  dphii.k.deta.k1 = matrix(dphii.k.deta.k1.tmp[[i]],ncol=K1),
#                                  dpii.k.deta.k = matrix(dpii.k.deta.k.tmp[[i]],ncol=K1))
#     }
#     names(subjectData) = uid
#     subjectData
# }
#
# logLikeCalc5 = function(params,
#                          yval,
#                          x,
#                          x.ppo=NULL,
#                          ppo.k=NULL,
#                          id,
#                          ProfileCol=NA,
#                          ref.muc=NA,
#                          TransIndMtx=NA, ## used in case there are 0 probability transition.  It should be a K-1 by K matrix with 1s in all cells except for those that correspond to 0 probability transtions.
#                          UseGrad=TRUE,
#                          CheeseCalc=FALSE){
#
#     States    = sort(unique(yval))
#     K         = length(States)
#     K1        = K-1
#     uid       = unique(id)
#     N         = length(uid)
#     Ntot      = length(yval)
#     npar      = length(params)
#     alpha.e   = params[1:K1]
#     beta.e    = params[(K1+1):(npar-(K1^2))]
#     gamma.e   = matrix(params[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
#     gamma.mtx = cbind(gamma.e, 0)
#     alpha.len = K1
#     beta.len  = length(beta.e)
#     gamma.len = K1^2
#
#     ## matrices and indices for partial PO model
#     if (!is.null(ppo.k)){
#         x.ppo0 = x.ppo*0
#         n.ppo = ncol(x.ppo)
#         tmp3a = matrix(0, K1, n.ppo)
#     }
#
#     if (!is.matrix(TransIndMtx)){ TransIndMtx = matrix(1, K1, K)}
#
#     ## Z and Y matrices
#     ZMat = YMat = matrix(0, Ntot, K)
#     for (k in 1:K){ YMat[,k] = ifelse(yval==k, 1, 0)
#                     ZMat[,k] = ifelse(yval<=k, 1, 0)
#     }
#
#     ## Marginal cumulative probs and state probs
#     ## The else part regards partial PO
#     cprob.m = prob.m = matrix(NA, Ntot, K)
#     for (k in 1:K1){
#         if (is.null(ppo.k)){cprob.m[,k]   = expit(alpha.e[k] + x %*% beta.e)
#         }else{  x.ppo.k       = x.ppo0
#                 aaa           = which(ppo.k==k)
#                 x.ppo.k[,aaa] = x.ppo[,aaa]
#                 cprob.m[,k]   = expit(alpha.e[k] + cbind(x,x.ppo.k) %*% beta.e)
#         }
#     }
#     cprob.m[,K] = 1
#     prob.m[,1]  = cprob.m[,1]
#     for (k in 2:K){
#         prob.m[,k] = cprob.m[,(k)]- cprob.m[,(k-1)]
#     }
#     one.min.cprob.m = 1-cprob.m
#     one.min.prob.m  = 1-prob.m
#
#     ## Phi and g(Phi) (we only use the first observation for each subject)
#     phi.k = matrix(NA, Ntot, K1)
#     for (k in 1:K1){
#         phi.k[,k] = log(cprob.m[,k]/(prob.m[,(k+1)]))
#     }
#     g.phi.k     = log(1+exp(phi.k))
#     exp.g.phi.k = exp(g.phi.k)
#
#     ## dPhi/deta (we only use the first observation for each subject)
#     ## dpi/deta
#     dphi.k.deta.k = dphi.k.deta.k1 = dpi.k.deta.k = matrix(NA, Ntot, K1)
#     for (k in 1:K1){
#         dphi.k.deta.k[,k]  =  exp.g.phi.k[,k]*one.min.cprob.m[,k]
#         dphi.k.deta.k1[,k] =  -exp.g.phi.k[,k]*one.min.cprob.m[,(k+1)]
#         dpi.k.deta.k[,k]   =  cprob.m[,k]*one.min.cprob.m[,k]
#     }
#
#     ## Some useful matrices to be used in gradient calculations
#     tmp1  = diag(1,K1)
#     tmp2  = rbind(tmp1[-1,], rep(0,K1))
#     tmp4  = matrix(0, K1, K1^2)
#
#     ## To be used in the gradient calculation
#     tmp.mat = matrix(0, K1, K1^2)
#     for (ppp in 1:K1){
#         tmp.mat[ppp,(c(1:K1)+(ppp-1)*K1)] = 1
#     }
#
#     li1 = li2 = rep(0,N)
#     Grad1Vec = Grad2Vec = rep(0,npar)
#     Grad1Mat = Grad2Mat = matrix(0,N, npar)
#
#     subjectData = CreateSubjectData(id, yval,x,x.ppo, ZMat,YMat,cprob.m,prob.m,one.min.cprob.m,phi.k, g.phi.k,dphi.k.deta.k,dphi.k.deta.k1,dpi.k.deta.k,K,K1)
#     blah        = lapply(subjectData, LogLikeiCalc5, gamma.mtx=gamma.mtx, tmp1=tmp1, tmp2=tmp2, tmp3a=tmp3a, tmp4=tmp4, tmp.mat=tmp.mat, ref.muc=ref.muc,
#                          K=K, K1=K1, n.ppo=n.ppo, ppo.k=ppo.k, npar=npar, alpha.len=alpha.len, beta.len=beta.len,
#                          gamma.len=gamma.len, TransIndMtx=TransIndMtx)
#
#
#     li1 = lapply(blah, function(x) x[['logLi1']])
#     loglike1 = -1*sum(c(unlist(li1)))
#     li2 = lapply(blah, function(x) x[['logLi2']])
#     loglike2 = -1*sum(c(unlist(li2)))
#     gradi1 = lapply(blah, function(x) x[['Grad1iVec']])
#     gradient1 = -1*Reduce('+', gradi1)
#     gradi2 = lapply(blah, function(x) x[['Grad2iVec']])
#     gradient2 = -1*Reduce('+', gradi2)
#
#     gradient1[ProfileCol] = 0
#     gradient2[ProfileCol] = 0
#
#     #print(loglike1)
#     #print(loglike2)
#     #print(gradient1)
#     #print(gradient2)
#     loglikelihood  = loglike1+loglike2
#     if (UseGrad==TRUE) attr(loglikelihood,"gradient") = c(gradient1+gradient2)
#     if (CheeseCalc==TRUE){
#         gradi = mapply(function(a,b) a+b, gradi1, gradi2, SIMPLIFY=FALSE)
#         gradi.outer = mapply(outer, gradi, gradi, SIMPLIFY=FALSE)
#         attr(loglikelihood,"cheese") = Reduce("+", gradi.outer)
#     }
#     loglikelihood
# }
#
# LogLikeiCalc5 = function(subjectData, gamma.mtx, tmp1, tmp2, tmp3a, tmp4, tmp.mat, ref.muc, K, K1,
#                          n.ppo, ppo.k, npar, alpha.len, beta.len, gamma.len, TransIndMtx){
#
#         yival            = subjectData[["yival"]]
#         xi               = subjectData[["xi"]]
#         xi.ppo           = subjectData[["xi.ppo"]]
#         mi               = subjectData[["mi"]]
#         ZiMat            = subjectData[["ZiMat"]]
#         YiMat            = subjectData[["YiMat"]]
#         cprobi.m         = subjectData[["cprobi.m"]]
#         probi.m          = subjectData[["probi.m"]]
#         one.min.cprobi.m = subjectData[["one.min.cprobi.m"]]
#         phii.k           = subjectData[["phii.k"]]
#         g.phii.k         = subjectData[["g.phii.k"]]
#         dphii.k.deta.k   = subjectData[["dphii.k.deta.k"]]
#         dphii.k.deta.k1  = subjectData[["dphii.k.deta.k1"]]
#         dpii.k.deta.k    = subjectData[["dpii.k.deta.k"]]
#
#         if (mi==1){ xi       = matrix(xi,nrow=1)
#                     xi.ppo   = matrix(xi.ppo,nrow=1)
#                     ZiMat    = matrix(ZiMat,nrow=1)
#                     YiMat    = matrix(YiMat,nrow=1)
#                     cprobi.m = matrix(cprobi.m, nrow=1)
#                     probi.m  = matrix(probi.m, nrow=1)
#                     phii.k          = matrix(phii.k, nrow=1)
#                     g.phii.k        = matrix(g.phii.k, nrow=1)
#                     dphii.k.deta.k  = matrix(dphii.k.deta.k, nrow=1)
#                     dphii.k.deta.k1 = matrix(dphii.k.deta.k1, nrow=1)
#                     dpii.k.deta.k   = matrix(dpii.k.deta.k, nrow=1)
#         }
#
#         logLi1    = logLi2 = 0
#         logLij2   = rep(0, mi)
#         Grad1iVec = Grad2iVec = rep(0, npar)
#
#
#         ## LogLi1 and dLogLi1/d(alpha,beta)
#         ## LogLi1 only comes from the marginal portion of the model (first observation)
#         for (k in 1:K1){ #print(cprobi.m[1,]);print(probi.m[1,])
#             logLi1 = logLi1 +  ZiMat[1,k]*phii.k[1,k] - ZiMat[1,(k+1)]*g.phii.k[1,k]
#         }
#
#         ## LogLi1/dtheta
#         ## the else term involves the partial proportional odds part of the model
#         xi1     = xi[1,]
#         xi1.ppo = xi.ppo[1,]
#         tmp3    = matrix(rep(xi1,K1), nrow=K1, byrow=TRUE)
#
#         if (is.null(ppo.k)){ deta.k.dtheta  = cbind(tmp1,tmp3,tmp4)
#                              deta.k1.dtheta = cbind(tmp2,tmp3,tmp4)
#         }else{
#             tmp3b = tmp3a
#             for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xi1.ppo[rrr]}
#                                   deta.k.dtheta         = cbind(tmp1,tmp3,tmp3b, tmp4)
#                                   deta.k1.dtheta        = cbind(tmp2,tmp3,tmp3b, tmp4)
#         }
#
#         for (k in 1:K1){
#             Grad1iVec = Grad1iVec+(ZiMat[1,k]-(ZiMat[1,(k+1)]*cprobi.m[1,k]/cprobi.m[1,(k+1)]))*(dphii.k.deta.k[1,k]*deta.k.dtheta[k,] + dphii.k.deta.k1[1,k]*deta.k1.dtheta[k,])
#         }
#
#         if (mi>1){
#             ## logLi2 comes from observations 2 to mi
#             YiMat2 = YiMat[,-K]
#
#             ############ Altered for new ref.muc (reference state in conditional model)
#             if (!is.na(ref.muc)){
#                 probi.m   = cbind(probi.m[,-ref.muc], probi.m[,ref.muc])
#                 YiMat.tmp = cbind(YiMat[,-ref.muc], YiMat[,ref.muc])
#                 YiMat2    = YiMat.tmp[,-K]
#                 yival2    = ifelse(yival<ref.muc, yival,
#                                    ifelse(yival==ref.muc, K, yival-1))
#                 yival     = yival2
#             }
#
#             for (j in 2:mi){
#                 mm           = probi.m[j,]
#                 mm.lag       = probi.m[(j-1),]
#                 Yit          = YiMat2[j,]
#                 Yit1         = YiMat2[(j-1),]
#                 xit          = xi[j,]
#                 xit1         = xi[(j-1),]
#                 xit.ppo      = xi.ppo[j,]
#                 xit1.ppo     = xi.ppo[(j-1),]
#                 dpi.deta     = dpii.k.deta.k[j,]
#                 dpi.deta.lag = dpii.k.deta.k[(j-1),]
#
#                 #############################################################################################################
#                 #Deltaij2    = Delta_calc(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx)
#                 Deltaij2    = Delta_calcWith0TxProbs(mm=mm, mm.lag=mm.lag, gamma.mat=gamma.mtx, CalcTxMtx=TransIndMtx)
#                 #############################################################################################################
#
#                 #lpij.c     = Deltaij2 + t(gamma.e) %*% YiMat2[(j-1),]
#                 lpij.c      = Deltaij2 + gamma.mtx[,yival[(j-1)]] #%*% YiMat2[(j-1),]
#                 probij.cK   = 1/(1+sum(exp(lpij.c)))
#
#                 logLij2[j]  = Yit %*% lpij.c + log(probij.cK)
#                 muit.c      = exp(lpij.c)*probij.cK
#
#                 #############################################################################################################
#                 #hmat = hmat_calc(Deltaij2, gamma.mtx)
#                 hmat = hmat_calc(Deltaij2, gamma.mtx)*TransIndMtx
#
#                 #############################################################################################################
#
#                 # df.dDelta   = matrix(0, K1, K1)
#                 # ## Left side of system of equations that solve for dDelta/dtheta
#                 # for (l in 1:K1){df.dDelta[l,l] = sum(hmat[l,]*(1-hmat[l,])*mm.lag)}
#                 # ## upper triangle for df.dDelta
#                 # for (l in 1:(K1-1)){
#                 #     for (m in (l+1):K1){
#                 #         df.dDelta[l,m] = -sum(hmat[l,]*hmat[m,]*mm.lag)
#                 #     }}
#                 # ## lower triangle for
#                 # df.dDelta[lower.tri(df.dDelta)] =  t(df.dDelta)[lower.tri(df.dDelta)]
#                 df.dDelta = dfdDelta_calc(mm.lag=mm.lag, hmat=hmat)
#
#                 ## right side of system of equations that solve for dDelta/dtheta (only alpha and beta)
#                 ## the else part is for partial PO
#                 if (is.null(ppo.k)){
#                     dpidtheta.lag = dpidtheta_calc1(tmp1, dpi.deta.lag, xit1)
#                     dpidtheta     = dpidtheta_calc1(tmp1, dpi.deta, xit)
#                 }else{
#                     tmp3b = tmp3a
#                     tmp3c = tmp3a
#                     for (rrr in 1:n.ppo){ tmp3b[ppo.k[rrr],rrr] = xit1.ppo[rrr]
#                                           tmp3c[ppo.k[rrr],rrr] = xit.ppo[rrr]
#                     }
#                     dpidtheta.lag = dpidtheta_calc2(tmp1, dpi.deta.lag, xit1, tmp3b)
#                     dpidtheta     = dpidtheta_calc2(tmp1, dpi.deta,     xit,  tmp3c)
#                 }
#                 ########################################################################
#                 dmumdtheta.lag = dmumdtheta_calc(dpidtheta.lag)
#                 ## for new ref.muc #############################
#                 if (!is.na(ref.muc)) dmumdtheta.lag = rbind(dmumdtheta.lag[-ref.muc,], dmumdtheta.lag[ref.muc,])
#                 h.dmumdtheta.lag     = hmat %*% dmumdtheta.lag
#
#                 ########################################################################
#                 dmumdtheta = dmumdtheta_calc(dpidtheta)
#                 ## for new ref.muc #############################
#                 if (is.na(ref.muc)){dmumdtheta  = dmumdtheta[-K,]
#                 }else{              dmumdtheta  = dmumdtheta[-ref.muc,]}
#
#                 ########################################################################
#                 rhs.alpha.beta  = dmumdtheta - h.dmumdtheta.lag
#                 dhdgamma.mm.lag = dhdgamma_mmlag_calc(hmat, mm.lag, tmp.mat)
#                 rhs             = cbind(rhs.alpha.beta, dhdgamma.mm.lag)
#                 ddelta.dtheta   = solve(df.dDelta, rhs)
#
#                 GradiAdd = c((Yit-muit.c) %*% ddelta.dtheta + c(rep(0, alpha.len+beta.len), c(rep(Yit-muit.c, each=K1) * rep(Yit1, K1))))
#                 Grad2iVec = Grad2iVec + GradiAdd
#             }
#         }
#
#         out = list(logLi1=logLi1,
#                    logLi2=sum(logLij2),
#                    Grad1iVec=Grad1iVec,
#                    Grad2iVec=Grad2iVec)
#         out
# }
#
# CalcVarCov = function(MOD, epsilon, yval, x, id, ProfileCol=NA, ref.muc=NA, TransIndMtx=NA, x.ppo=NULL, ppo.k=NULL){
#     npar        = length(MOD$estimate)
#     eps.mtx     = diag(rep(epsilon, npar))
#     grad.at.max = MOD$gradient
#     ObsInfo.tmp = ObsInfo = matrix(NA, npar, npar)
#
#     ## Observed Information
#     for (j in 1:npar){
#         temp            = logLikeCalc5(MOD$estimate+eps.mtx[j,], yval=yval, x=x, x.ppo=x.ppo, ppo.k=ppo.k, id=id, ProfileCol = ProfileCol, ref.muc=ref.muc, TransIndMtx=TransIndMtx, UseGrad=TRUE)
#         ObsInfo.tmp[j,] = (attr(temp,"gradient")-grad.at.max)/(epsilon)
#     }
#     for (m in 1:npar){ for (n in 1:npar){ ObsInfo[m,n] =  (ObsInfo.tmp[m,n]+ObsInfo.tmp[n,m])/2}}
#
#     ## remove rows and columns of information matrix with no effective information.  This is usually caused by absorbing states
#     rmv.noInfo = which(diag(ObsInfo)<.01)
#
#     ## Model based covariance
#     if (length(rmv.noInfo)>0){ print(rmv.noInfo)
#                                mod.cov = solve(ObsInfo[-rmv.noInfo,-rmv.noInfo])
#     }else{                     mod.cov = solve(ObsInfo)}
#
#     ## Robust, sandwich covariance
#     cheese  = attr(logLikeCalc5(params=MOD$estimate, yval=yval, x=x, x.ppo=x.ppo, ppo.k=ppo.k, id=id, ProfileCol=ProfileCol, ref.muc=ref.muc, TransIndMtx=TransIndMtx, CheeseCalc=TRUE), "cheese")
#     if (length(rmv.noInfo)>0){
#         cheese  = cheese[-rmv.noInfo,-rmv.noInfo]
#     }
#     rob.cov = mod.cov %*% cheese %*% mod.cov
#
#     out = list(mod.cov=mod.cov, rob.cov=rob.cov, cheese=cheese)
#     out
#
# }
