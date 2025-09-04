####################################################################################################
##### Create a Function to summarize the OMTM results ##############################################
####################################################################################################
omtm1  <- function(params, yval, XMat, id, ppo.k=NULL, x.ppo=NULL, u=NULL,
                   weight = NULL,
                   ref.muc = NA, TransIndMtx=NA,
                   stepmax = 0.5, UseGrad=TRUE, iterlim=250, check.analyticals = FALSE,
                   print.level=0, ProfileCol = NA){

    # if weights is null assign a weight of 1
    if(is.null(weight)){weight <- rep(1, length(yval))}

    # fit the model
    mod <- nlm(logLikeCalc, p = params, yval = yval,  x = XMat, wt = weight,
               id = id, UseGrad=UseGrad,
               ppo.k = ppo.k, x.ppo = x.ppo, u=u, ref.muc = ref.muc, TransIndMtx=TransIndMtx,
               stepmax=stepmax, iterlim=iterlim, check.analyticals = check.analyticals,
               print.level=print.level, ProfileCol = ProfileCol)

    # variance and covariance matrix
    tmp           <- CalcVarCov(MOD = mod, epsilon = 1e-7, yval = yval, x = XMat, wt = weight,
                                id = id, ProfileCol = NA, ref.muc = ref.muc,
                                ppo.k = ppo.k, x.ppo = x.ppo, u=u, TransIndMtx=TransIndMtx)

    # get the estimated parameter from the nlm call
    est.par     <- mod$estimate
    npar        <- length(est.par)
    K1          <- length(unique(yval)) - 1
    n.alpha     <- K1
    n.beta.ppo  <- ifelse(is.null(x.ppo), 0, ncol(x.ppo))
    n.beta      <- ncol(XMat)+n.beta.ppo
    n.gamma     <- npar - n.beta - n.alpha

    alpha       <- est.par[1:n.alpha]
    beta        <- est.par[(n.alpha+1):(n.alpha+n.beta)]

    ##
    gamma.vec    <- est.par[(n.alpha+n.beta+1):npar]
    gamma.vec.sp <- split(gamma.vec, ceiling(seq_along(gamma.vec)/(K1^2)))
    gamma        <- lapply(gamma.vec.sp, function(x,...) matrix(x, K1, K1, byrow=TRUE))

    # output
    out          <- NULL
    out$loglik   <- -mod$minimum
    out$alpha    <- alpha
    out$beta     <- beta
    out$gamma    <- gamma
    out$vcov     <- tmp$mod.cov
    out$rob.vcov <- tmp$rob.cov
    out$no.info  <- tmp$low.info


    out$control  <- c(mod$code, mod$iterations, length(table(id)), length(id),
                      ref.c = ifelse(is.na(ref.muc), K1 + 1, ref.muc))
    names(out$control) <- c("convergence_code", "n_iter", "n_subjects", "n_obs",
                            "cond_ref")

    # create class omtm for marginalized transition models
    class(out)    <- "omtm"

    # output
    return(out)
}


# For now the function uses as inputs the ouputs from nlm and CalcVarCov (this can change if nlm and
# CalcVarCov are merged in a unique function)
#
# omtm  <- function(log.lik, params, yval, XMat, id, ppo.k, x.ppo, ref.muc = NA, stepmax = 0.25){
#
#   # fit the model
#   mod <- nlm(log.lik, p = params, yval = yval,  x = XMat,
#              id = id, UseGrad=TRUE, ppo.k = ppo.k, x.ppo = x.ppo, ref.muc = ref.muc,
#              stepmax=stepmax, iterlim=250, check.analyticals = FALSE, print.level=0, hessian=FALSE)
#
#   # variance and covariance matrix
#   tmp           <- CalcVarCov(MOD = mod, epsilon = 1e-7, yval = yval, x = XMat,
#                               id = id, ProfileCol = NA, ref.muc = ref.muc, TransIndMtx = NA,
#                               ppo.k = ppo.k, x.ppo = x.ppo)
#
#   # get the estimated parameter from the nlm call
#   est.par     <- mod$estimate
#   npar        <- length(est.par)
#   K1          <- length(unique(yval)) - 1
#   alpha       <- est.par[1:K1]
#   beta        <- est.par[(K1+1):(npar-(K1^2))]
#   gamma       <- matrix(est.par[(npar+1-(K1^2)):npar], K1, K1, byrow=TRUE)
#
#   # output
#   out          <- NULL
#   out$loglik   <- -mod$minimum
#   out$alpha    <- alpha
#   out$beta     <- beta
#   out$gamma    <- gamma
#   out$vcov     <- tmp$mod.cov
#   out$rob.vcov <- tmp$rob.cov
#
#
#   out$control  <- c(mod$code, mod$iterations, length(table(id)), length(id),
#                     ref.c = ifelse(is.na(ref.muc), K1 + 1, ref.muc))
#   names(out$control) <- c("convergence_code", "n_iter", "n_subjects", "n_obs",
#                           "cond_ref")
#
#   # create class omtm for marginalized transition models
#   class(out)    <- "omtm"
#
#   # output
#   return(out)
#   }

summary.omtm <-
    function(object, robust = FALSE,...) {
        mod             <- object
        K1              <- length(mod$alpha)
        ref.c           <- mod$control["cond_ref"]
        ests            = c(mod$alpha, mod$beta, c(unlist(lapply(mod$gamma, function(x) t(x)))))
        se.mod          = sqrt(diag(mod$vcov))
        se.rob          = sqrt(diag(mod$rob.vcov))
        ind             = mod$no.info
        l.gamma         = length(mod$gamma)

        if (length(ests) != length(se.mod)){
            for (i in 1:length(ind)){ se.mod = append(se.mod, 1e8, after=ind[i]-1)
                                      se.rob = append(se.rob, 1e8, after=ind[i]-1)
            }
        }
        z.stat.mod = ests/se.mod
        z.stat.rob = ests/se.rob

        p.mod         = 2*pnorm(abs(z.stat.mod), lower.tail = FALSE)
        p.rob         = 2*pnorm(abs(z.stat.rob), lower.tail = FALSE)

        tab.mod = data.frame(Estimate=ests, StdErr=se.mod, Z = z.stat.mod, "p-value" = p.mod)
        tab.rob = data.frame(Estimate=ests, StdErr=se.rob, Z = z.stat.rob, "p-value" = p.rob)


        row.names.1 <- c(paste0("Intercept:Y<=", 1:K1), paste0("b_",1:length(mod$beta)))

        row.names.2.tmp <- paste0("gamma:", rep(1:(K1+1), each = (K1+1)), 1:(K1+1))
        if (l.gamma>1){
            for (LL in 2:l.gamma){
                row.names.2.tmp <- c(row.names.2.tmp,
                                     paste0("gamma:", rep(1:(K1+1), each = (K1+1)), 1:(K1+1), paste0(rep("'",LL-1), collapse="")))
            }
        }
        row.names.2     <- row.names.2.tmp[-grep(ref.c, row.names.2.tmp)]

        row.names(tab.mod) = row.names(tab.rob) = c(row.names.1, row.names.2)

        L = length(c(mod$alpha, mod$beta))
        cs.tab.mod  = tab.mod[1:L,]
        cs.tab.rob  = tab.rob[1:L,]
        dep.tab.mod = tab.mod[-c(1:L),]
        dep.tab.rob = tab.rob[-c(1:L),]

        if(robust == FALSE){
            out = list(class = class(mod),  control=mod$control,
                       cs.table = tab.mod[1:L,], dependence.table = tab.mod[-c(1:L),])
        }
        if(robust == TRUE){
            out = list(class = class(mod),  control=mod$control,
                       cs.table = tab.rob[1:L,], dependence.table = tab.rob[-c(1:L),])
        }
        class(out) <- "summary.omtm"
        return(out)
    }


#
#
# summary.omtm <-
#   function(object, robust = FALSE,...) {
#     mod             <- object
#     K1              <- length(mod$alpha)
#     ref.c           <- mod$control["cond_ref"]
#     if(robust == FALSE){
#       # Mean Model
#       m.mod              <- c(mod$alpha, mod$beta)
#       mean.res           <- data.frame(estimate = m.mod,
#                                        mod.se = sqrt(diag(mod$vcov))[seq_along(m.mod)])
#       rownames(mean.res) <- c(paste0("Intercept:Y<=", 1:K1), paste0("b_",1:length(mod$beta)))
#       mean.res$z         <- mean.res$estimate/mean.res$mod.se
#       mean.res$p         <- 2*pnorm(abs(mean.res$z), lower.tail = FALSE)
#
#       names(mean.res)    <- c('Estimate','Model.SE')
#       names(mean.res) = c('Estimate','Model.SE', 'z','p-value')
#       # Conditional Model
#       #c.mod              <- c(unlist(t(mod$gamma)))
#       c.mod              <- c(unlist(lapply(mod$gamma, function(x) t(x))))
#
#       if(nrow(mod$vcov) == length(c(m.mod, c.mod))){
#         cond.res           <- data.frame(estimate = c.mod,
#                                          mod.se = c(sqrt(diag(mod$vcov))[(length(m.mod) + 1):
#                                                                            nrow(mod$vcov)]))
#       }
#       if(nrow(mod$vcov) != length(c(m.mod, c.mod))){
#       cond.res           <- data.frame(estimate = c.mod,
#                                     mod.se = c(rep(NA, K1),
#                                                sqrt(diag(mod$vcov))[(length(m.mod) + 1):
#                                                                       nrow(mod$vcov)]))
#       }
#       cond.res$z         <- cond.res$estimate/cond.res$mod.se
#       cond.res$p         <- 2*pnorm(abs(cond.res$z), lower.tail = FALSE)
#
#       all.comb           <- paste0("gamma:", rep(1:(K1+1), each = (K1+1)), 1:(K1+1))
#       #all.comb           <- paste0("gamma:", 1:length(c.mod))
#       ref.comb.remove    <- all.comb[-grep(ref.c, all.comb)]
#       rownames(cond.res) <- ref.comb.remove
#       names(cond.res) = c('Estimate','Model.SE', 'z','p-value')
#       out = list(class = class(object),  control=object$control,
#                  mean.table = mean.res, cond.table = cond.res)
#     }
#     if(robust == TRUE){
#       # Mean Model
#       m.mod              <- c(mod$alpha, mod$beta)
#       c.mod              <- c(mod$gamma)
#       mean.res           <- data.frame(estimate = m.mod,
#                                        mod.se = sqrt(diag(mod$rob.vcov))[seq_along(m.mod)])
#       mean.res$z         <- mean.res$estimate/mean.res$mod.se
#       mean.res$p         <- 2*pnorm(abs(mean.res$z), lower.tail = FALSE)
#
#       names(mean.res)    <- c('Estimate','Model.SE')
#       names(mean.res) = c('Estimate','Rob.Model.SE', 'z','p-value')
#       rownames(mean.res) <- c(paste0("Intercept:Y<=", 1:K1), paste0("b_",1:length(mod$beta)))
#       # Conditional Model
#       #c.mod              <- c(mod$gamma)
#       c.mod              <- c(unlist(lapply(mod$gamma, function(x) t(x))))
#
#       if(nrow(mod$vcov) == length(c(m.mod, c.mod))){
#         cond.res           <- data.frame(estimate = c.mod,
#                                          mod.se = sqrt(diag(mod$rob.vcov))[(length(m.mod) + 1):
#                                                                              nrow(mod$rob.vcov)])
#       }
#       if(nrow(mod$vcov) != length(c(m.mod, c.mod))){
#         cond.res           <- data.frame(estimate = c.mod,
#                                          mod.se = c(rep(NA, K1),
#                                                     sqrt(diag(mod$rob.vcov))[(length(m.mod) + 1):
#                                                                                nrow(mod$rob.vcov)]))
#       }
#       cond.res$z         <- cond.res$estimate/cond.res$mod.se
#       cond.res$p         <- 2*pnorm(abs(cond.res$z), lower.tail = FALSE)
#
#       all.comb           <- paste0("gamma:", rep(1:(K1+1), each = (K1+1)), 1:(K1+1))
#       #all.comb           <- paste0("gamma:", 1:length(c.mod))
#       ref.comb.remove    <- all.comb[-grep(ref.c, all.comb)]
#       rownames(cond.res) <- ref.comb.remove
#       names(cond.res) = c('Estimate','Rob.Model.SE', 'z','p-value')
#
#       out = list(class = class(object),  control=object$control,
#                  mean.table = mean.res, cond.table = cond.res)
#     }
#     class(out) <- "summary.omtm"
#     return(out)
#     }


print.summary.omtm <-
  function(x,...) {
    cat("\nCross-sectional Model Parameters:\n")
    printCoefmat(x$cs.table)
    cat("\nDependence Model Parameters:\n")
    printCoefmat(x$dependence.table)
    cat('\n')
    cat('Number of Observations:            ',x$control[4],'\n')
    cat('Number of Groups:                  ',x$control[3],'\n')
    cat('Convergence Status (nlm code):     ',x$control[1],'\n')
    cat('Number of Iterations:              ',x$control[2],'\n')
    cat('Conditional Model Reference Level: ',x$control[5])
    cat('\n')
  }


