#' @name survAWKMT2-package
#' @aliases  survAWKMT2-package
#' @docType  package
#' @title Two-Sample Tests Based on Weighted Differences of Kaplan-Meier Curves
#' @description
#' Tests for equality of two survival functions based on integrated weighted differences of two Kaplan-Meier curves.
#' @author Miki Horiguchi, Hajime Uno
#' @references
#' Uno H, Tian L, Claggett B, Wei LJ. A versatile test for equality of two survival functions based on weighted differences of Kaplan-Meier curves.
#' Statistics in Medicine 2015, 34, 3680-3695.
#' @seealso
#' survival
#' @examples
#' D        = survival::pbc[1:312, c(2,3,4)] #The pbc data from 'survival' package
#' D$status = as.numeric(D$status==2)
#' D$trt    = as.numeric(D$trt==2)
#' names(D) = c("time", "status", "arm")
#' tau      = max(D[D[,2]==1,1])
#' nmethod  = 10 #Recommended to specify at least 10000 (default) or larger.
#'
#' a = AWKMT2(indata=D, tau=tau, c_first=0, c_last=4, c_by=0.1, method="permutation",
#'            nmethod=nmethod, seed=1, v1=TRUE, v2=TRUE, test="1_side")
#' print(a)
NULL

#' Adaptively Weighted Kaplan-Meier Tests
#'
#' Performs the two-sample tests based on adaptively weighted differences between two Kaplan-Meier curves proposed by Uno, Tian, Claggett and Wei (2015).
#'
#'
#' @usage
#' AWKMT2(indata, tau, c_first=0, c_last=4, c_by=0.1, method="permutation",
#'        nmethod=10000, seed=1, v1=TRUE, v2=TRUE, test="1_side")
#'
#' @param indata          A data matrix (data frame). The 1st column is time-to-event variable, the 2nd column is event indicator (1=event, 0=censor), and the 3rd column is the treatment indicator (1=treat, 0=control). No missing values are allowed in this data matrix.
#' @param tau             A numeric value to specify the time interval of interest. The end of study time will be a general choice.
#' @param c_first         A first number in range to specify the search area of "c" for the versatile tests by Uno et al. (2015). Default is \code{0}.
#' @param c_last          A last number in range to specify the search area of "c" for the versatile tests by Uno et al. (2015). Default is \code{4}.
#' @param c_by            A number to specify the search area of "c" for the versatile tests by Uno et al. (2015). Default is \code{0.1}.
#' @param method          A name of the resampling method. It supports \code{"permutation"} (default) and \code{"perturbation"}.
#' @param nmethod         A number of iterations for the resampling. Recommended to specify at least \code{10000} (default) or larger.
#' @param seed            An integer value, used for the random number generation in the resampling procedures. Default is \code{1}.
#' @param v1              Choice of the test statistic. When \code{TRUE} (default), v1 proposed by Uno et al. (2015) is used as a test statistic.
#' @param v2              Choice of the test statistic. When \code{TRUE} (default), v2 proposed by Uno et al. (2015) is used as a test statistic.
#' @param test            Specify \code{"1_side"} for the one-sided test where the alternative hypothesis is that treatment group is superior to control group with respect to survival.
#'                        Specify \code{"2_side"} for the two-sided test where the alternative hypothesis is that treatment group is not equal to control group with respect to survival.
#'                        Default is \code{"1_side"}.
#' @return A list with components:
#' @return \item{resampling_method}{The resampling method.}
#' @return \item{crude_pvalue_T1_1_side}{The one-sided crude p-value of the test based on v1 in Uno et al. (2015).}
#' @return \item{crude_pvalue_T2_1_side}{The one-sided crude p-value of the test based on v2 in Uno et al. (2015).}
#' @return \item{crude_pvalue_T1_2_side}{The two-sided crude p-value of the test based on v1 in Uno et al. (2015).}
#' @return \item{crude_pvalue_T2_2_side}{The two-sided crude p-value of the test based on v2 in Uno et al. (2015).}
#' @return \item{bona_fide_pvalue_T1_1_side}{The one-sided bona-fide p-value of the test based on v1 in Uno et al. (2015).}
#' @return \item{bona_fide_pvalue_T2_1_side}{The one-sided bona-fide p-value of the test based on v2 in Uno et al. (2015).}
#' @return \item{bona_fide_pvalue_T1_2_side}{The two-sided bona-fide p-value of the test based on v1 in Uno et al. (2015).}
#' @return \item{bona_fide_pvalue_T2_2_side}{The two-sided bona-fide p-value of the test based on v2 in Uno et al. (2015).}
#' @seealso
#'  survival
#' @references
#'  Uno H, Tian L, Claggett B, Wei LJ. A versatile test for equality of two survival functions based on weighted differences of Kaplan-Meier curves.
#'  Statistics in Medicine 2015, 34, 3680-3695.
#' @name   AWKMT2
#' @aliases AWKMT2
#' @import survival
#' @importFrom stats rexp rmultinom
#' @examples
#'  D        = survival::pbc[1:312, c(2,3,4)] #The pbc data from 'survival' package
#'  D$status = as.numeric(D$status==2)
#'  D$trt    = as.numeric(D$trt==2)
#'  names(D) = c("time", "status", "arm")
#'  tau      = max(D[D[,2]==1,1])
#'  nmethod  = 10 #Recommended to specify at least 10000 (default) or larger.
#'
#'  a = AWKMT2(indata=D, tau=tau, c_first=0, c_last=4, c_by=0.1, method="permutation",
#'             nmethod=nmethod, seed=1, v1=TRUE, v2=TRUE, test="1_side")
#'  print(a)
#'
NULL

###############################################################################################
# survAWKMT2 package (version 1.0.0)
#
#
# funcKM2 (hidden)
# funcAWKMT2_c1 (hidden)
# AWKMT2 (main function)
###############################################################################################

###############################################################################################
# funcKM2 (hidden)
###############################################################################################
funcKM2 <- function(data, t_idx, status, x, weight=NULL){
  if(is.null(weight)){
    weight = rep(1, nrow(data))
  }

  data     = data[order(data$time), ]
  data_wt  = cbind(data, weight)

  Y = rep(0,length(t_idx)) #n.risk
  N = rep(0,length(t_idx)) #n.event
  C = rep(0,length(t_idx)) #n.censor
  S = rep(0,length(t_idx)) #surv
  H = rep(0,length(t_idx)) #forSE
  D = rep(0,length(t_idx)) #forSE
  E = rep(0,length(t_idx)) #SE

  #i=1
  Y[1] = sum(data_wt[ ,4])
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0

  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]

    if(sum(x==t_idx[i] & status==1)>0){
      N[i] = sum(x==t_idx[i] & status==1)*data_wt[which(x==t_idx[i] & status==1)[1],4]
    }else{
      N[i] = 0
    }

    if(sum(x==t_idx[i] & status==0)>0){
      C[i] = sum(x==t_idx[i] & status==0)*data_wt[which(x==t_idx[i] & status==0)[1],4]
    }else{
      C[i] = 0
    }

    if(Y[i]<0) Y[i] = 0

    if(Y[i]==0){
      S[i] = S[i-1]
    }else{
      S[i] = S[i-1]*(1-(N[i]/Y[i]))
    }

    if(Y[i]*(Y[i]-N[i])==0){
      H[i] = 0
    }else{
      H[i] = N[i]/(Y[i]*(Y[i]-N[i]))
    }

    if(S[i]<0) S[i] = 0

    D[i] = sum(H[2:i])
    E[i] = sqrt((S[i]**2)*D[i])

    if(is.na(S[i])) S[i] = 0
    if(is.na(E[i])) E[i] = 0
  }

  result = list()
  result$t_idx    = t_idx
  result$n_risk   = Y
  result$n_event  = N
  result$n_censor = C
  result$surv     = S
  result$SE       = E

  out           = cbind(result$t_idx, result$n_risk, result$n_event, result$n_censor, result$surv, result$SE)
  colnames(out) = c("t_idx", "n_risk", "n_event", "n_censor", "surv", "SE")
  return(out)
}
NULL

####################################################################################################################
# funcAWKMT2_c1 (hidden)
####################################################################################################################
funcAWKMT2_c1<- function(indata, t_idx, tau1, tau2, crange, test, type, obs_survdiff=NULL){
  #-- Get status1 and x1 (arm=1) --
  g1      = indata[indata$arm==1, ]
  g1      = g1[order(g1$time), ]
  status1 = g1$status
  x1      = g1$time

  #-- Get status0 and x0 (arm=0) --
  g0      = indata[indata$arm==0, ]
  g0      = g0[order(g0$time), ]
  status0 = g0$status
  x0      = g0$time

  #-- Get stats1 using funcKM2 --
  if(type=="observe" | type=="permutation"){
    wk1 = funcKM2(data=g1, t_idx=t_idx, status=status1, x=x1, weight=NULL)
    wk0 = funcKM2(data=g0, t_idx=t_idx, status=status0, x=x0, weight=NULL)
  }

  if(type=="perturbation"){
    n1  = nrow(g1)
    n0  = nrow(g0)
    wt1 = rexp(n1)
    wt0 = rexp(n0)
    wk1 = funcKM2(g1, t_idx, status1, x1, weight=wt1)
    wk0 = funcKM2(g0, t_idx, status0, x0, weight=wt0)
  }

  wk           =  cbind(wk1, wk0[,-1])
  colnames(wk) = c("t_idx", "n_risk1", "n_event1", "n_censor1", "surv1", "SE1", "n_risk0", "n_event0", "n_censor0", "surv0", "SE0")

  #-- Get stats2 --
  times       = wk[ ,"t_idx"]
  dt_diff     = c(0, diff(times))
  dt_jump     = wk[ ,"n_event1"] + wk[ ,"n_event0"]
  tau_idx     = times > tau1 & times <= tau2
  wk          = cbind(wk, dt_diff, dt_jump, tau_idx)
  wk          = wk[wk[ ,"tau_idx"]==1,]

  survdiff    = wk[ ,"surv1"] - wk[ ,"surv0"]
  survdiff_se = sqrt(wk[ ,"SE1"]^2 + wk[ ,"SE0"]^2)

  if(type=="perturbation"){
    survdiff = survdiff - obs_survdiff
  }

  #-- Get Z1(one-sided) and Z2(two-sided) --
  z1            = survdiff/survdiff_se
  z1[is.na(z1)] = 0
  z2            = abs(z1)

  wk = cbind(wk, survdiff, survdiff_se, z1, z2)
  n_crange = length(crange)

  VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
  }

  #-- Get V1 and V2 (one-sided) --
  if(test=="1_side"){
    tmp.Z   = VTM(z1, n_crange) #Z1
    tmp.cc  = t(VTM(crange, length(z1)))
    tmp.wt  = pmax(tmp.Z, tmp.cc)
    tmp.wtz = tmp.wt * tmp.Z
    #--- integration with respect to time ---
    T1_1_side = tmp.wtz %*% wk[ ,"dt_diff"] #V1
    #--- integration with respect to jump ---
    T2_1_side = tmp.wtz %*% wk[ ,"dt_jump"] #V2

    if(type=="observe"){
      OUT                 = list()
      OUT$stats           = cbind(T1_1_side, T2_1_side)
      OUT$obs_survdiff    = survdiff
    }else{
      OUT = cbind(T1_1_side, T2_1_side)
    }
  }

  #-- Get V1 and V2 (two-sided) --
  if(test=="2_side"){
    tmp.Z2   = VTM(z2, n_crange) #Z1
    tmp.cc2  = t(VTM(crange, length(z2)))
    tmp.wt2  = pmax(tmp.Z2, tmp.cc2)
    tmp.wtz2 = tmp.wt2 * tmp.Z2
    #--- integration with respect to time ---
    T1_2_side = tmp.wtz2 %*% wk[ ,"dt_diff"] #V1
    #--- integration with respect to jump ---
    T2_2_side = tmp.wtz2 %*% wk[ ,"dt_jump"] #V2

    if(type=="observe"){
      OUT                 = list()
      OUT$stats           = cbind(T1_2_side, T2_2_side)
      OUT$obs_survdiff    =  survdiff
    }else{
      OUT = cbind(T1_2_side, T2_2_side)
    }
  }
  return(OUT)
}
NULL

#' @export
################################################################################################################################################
# AWKMT2 (main function)
################################################################################################################################################
AWKMT2 <- function(indata, tau, c_first=0, c_last=4, c_by=0.1, method="permutation", nmethod=10000, seed=1, v1=TRUE, v2=TRUE, test="1_side"){
  rank2 <- function(x){
    rank(x, ties.method="min")
  }

  #-- Get tau1 --
  tau1 = 0

  #-- Get tau2 --
  tau2 = tau

  #-- Get unique_time, n_times --
  t_idx   = unique(sort(c(indata$time, 0, tau1, tau2, tau)))
  n_times = length(t_idx)

  #-- Get crange --
  crange   = seq(c_first, c_last, by=c_by)

  #-- Get observed V1(c) and V2(c) using funcAWKMT2_c1 --
  if(test=="1_side"){
    obs_ours = funcAWKMT2_c1(indata, t_idx, tau1, tau2, crange, test="1_side", type="observe")
  }

  if(test=="2_side"){
    obs_ours = funcAWKMT2_c1(indata, t_idx, tau1, tau2, crange, test="2_side", type="observe")
  }

  obs_survdiff    = obs_ours$obs_survdiff

  if(v1==TRUE & v2==TRUE){
    #-- Get the null distribution of V1(c) and V2(c) using resampling method --
    aa_T1 = obs_ours$stats[ ,1]
    aa_T2 = obs_ours$stats[ ,2]

    bb_T1 = matrix(NA, nrow=nmethod, ncol=length(crange))
    bb_T2 = matrix(NA, nrow=nmethod, ncol=length(crange))

    if(method=="permutation"){
      set.seed(seed)
      for (i in 1:nmethod){
        perm_arm            = indata[sample(1:nrow(indata), size=nrow(indata), replace=FALSE), 3]
        perm_data           = data.frame(indata[,1:2], perm_arm)
        colnames(perm_data) = c("time", "status", "arm")

        if(test=="1_side"){
          bb = funcAWKMT2_c1(perm_data, t_idx, tau1, tau2, crange, "1_side", type="permutation")
        }

        if(test=="2_side"){
          bb = funcAWKMT2_c1(perm_data, t_idx, tau1, tau2, crange, "2_side", type="permutation")
        }

        bb_T1[i,] = bb[ ,1]
        bb_T2[i,] = bb[ ,2]
      }
    }

    if(method=="perturbation"){
      n = length(indata$time)
      set.seed(seed)
      for (i in 1:nmethod){
        if(test=="1_side"){
          bb = funcAWKMT2_c1(indata, t_idx, tau1, tau2, crange, test="1_side", type="perturbation", obs_survdiff)
        }

        if(test=="2_side"){
          bb = funcAWKMT2_c1(indata, t_idx, tau1, tau2, crange, test="2_side", type="perturbation", obs_survdiff)
        }

        bb_T1[i,] = bb[ ,1]
        bb_T2[i,] = bb[ ,2]
      }
    }

    #-- Get crude p-value --
    cc_T1 = rbind(aa_T1, bb_T1)
    cc_T2 = rbind(aa_T2, bb_T2)

    #-- T1 --
    dd_T1    = apply(-cc_T1, 2, rank2)-1
    d2_T1    = dd_T1[1, ]
    crude_T1 = min(d2_T1/nmethod)

    #-- T2 --
    dd_T2    = apply(-cc_T2, 2, rank2)-1
    d2_T2    = dd_T2[1, ]
    crude_T2 = min(d2_T2/nmethod)

    #-- Get the null distribution of pb to choose the threshold value for claiming a statistical significance based on pb --
    #-- T1 --
    ee_T1        = apply(-bb_T1, 2, rank2)-1
    ours_bona_T1 = apply(ee_T1, 1, min)/(nmethod-1)

    aa_crude_pvalue_T1  = rep(crude_T1, nmethod)
    ours_bona_pvalue_T1 = mean(ours_bona_T1 < aa_crude_pvalue_T1)

    #-- T2 --
    ee_T2 = apply(-bb_T2, 2, rank2)-1
    ours_bona_T2 = apply(ee_T2, 1, min)/(nmethod-1)

    aa_crude_pvalue_T2  = rep(crude_T2, nmethod)
    ours_bona_pvalue_T2 = mean(ours_bona_T2 < aa_crude_pvalue_T2)


    ours_bona_pvalue = c(ours_bona_pvalue_T1, ours_bona_pvalue_T2)
  }else{
    #-- Get the null distribution of V1(c) or V2(c) using resampling method --
    if(v1==TRUE & v2==FALSE){
      aa_T = obs_ours$stats[ ,1]
    }

    if(v1==FALSE & v2==TRUE){
      aa_T = obs_ours$stats[ ,2]
    }

    bb_T = matrix(NA, nrow = nmethod, ncol = length(crange))

    if(method=="permutation"){
      set.seed(seed)
      for (i in 1:nmethod){
        perm_arm            = indata[sample(1:nrow(indata), size = nrow(indata), replace = FALSE), 3]
        perm_data           = data.frame(indata[ ,1:2], perm_arm)
        colnames(perm_data) = c("time", "status", "arm")

        if(test=="1_side"){
          bb = funcAWKMT2_c1(perm_data, t_idx, tau1, tau2, crange, "1_side", type="permutation")
        }

        if(test=="2_side"){
          bb = funcAWKMT2_c1(perm_data, t_idx, tau1, tau2, crange, "2_side", type="permutation")
        }

        if(v1==TRUE & v2==FALSE){
          bb_T[i, ] = bb[ ,1]
        }

        if(v1==FALSE & v2==TRUE){
          bb_T[i, ] = bb[ ,2]
        }
      }
    }

    if(method=="perturbation"){
      n = length(indata$time)
      set.seed(seed)
      for (i in 1:nmethod){
        if(test=="1_side"){
          bb = funcAWKMT2_c1(indata, t_idx, tau1, tau2, crange, test="1_side", type="perturbation", obs_survdiff)
        }

        if(test=="2_side"){
          bb = funcAWKMT2_c1(indata, t_idx, tau1, tau2, crange, test="2_side", type="perturbation", obs_survdiff)
        }

        if(v1==TRUE & v2==FALSE){
          bb_T[i, ] = bb[ ,1]
        }

        if(v1==FALSE & v2==TRUE){
          bb_T[i, ] = bb[ ,2]
        }
      }
    }

    #-- Get crude p-value --
    cc_T = rbind(aa_T, bb_T)

    #-- T --
    dd_T    = apply(-cc_T, 2, rank2)-1
    d2_T    = dd_T[1, ]
    crude_T = min(d2_T/nmethod)

    #-- Get the null distribution of pb to choose the threshold value for claiming a statistical significance based on pb --
    #-- T --
    ee_T        = apply(-bb_T, 2, rank2)-1
    ours_bona_T = apply(ee_T, 1, min)/(nmethod-1)

    aa_crude_pvalue_T  = rep(crude_T, nmethod)
    ours_bona_pvalue_T = mean(ours_bona_T < aa_crude_pvalue_T)

    ours_bona_pvalue   = ours_bona_pvalue_T
  }

  #-- Output --
  output = list()
  if(test=="1_side"){
    if(v1==TRUE & v2==TRUE){
      output$resampling_method          = method
      output$crude_pvalue_T1_1_side     = crude_T1
      output$crude_pvalue_T2_1_side     = crude_T2
      output$bona_fide_pvalue_T1_1_side = ours_bona_pvalue_T1
      output$bona_fide_pvalue_T2_1_side = ours_bona_pvalue_T2
    }
    if(v1==TRUE & v2==FALSE){
      output$resampling_method          = method
      output$crude_pvalue_T1_1_side     = crude_T
      output$bona_fide_pvalue_T1_1_side = ours_bona_pvalue_T
    }
    if(v1==FALSE & v2==TRUE){
      output$resampling_method          = method
      output$crude_pvalue_T2_1_side     = crude_T
      output$bona_fide_pvalue_T2_1_side = ours_bona_pvalue_T
    }
  }

  if(test=="2_side"){
    if(v1==TRUE & v2==TRUE){
      output$resampling_method          = method
      output$crude_pvalue_T1_2_side     = crude_T1
      output$crude_pvalue_T2_2_side     = crude_T2
      output$bona_fide_pvalue_T1_2_side = ours_bona_pvalue_T1
      output$bona_fide_pvalue_T2_2_side = ours_bona_pvalue_T2
    }
    if(v1==TRUE & v2==FALSE){
      output$resampling_method          = method
      output$crude_pvalue_T1_2_side     = crude_T
      output$bona_fide_pvalue_T1_2_side = ours_bona_pvalue_T
    }
    if(v1==FALSE & v2==TRUE){
      output$resampling_method          = method
      output$crude_pvalue_T2_2_side     = crude_T
      output$bona_fide_pvalue_T2_2_side = ours_bona_pvalue_T
    }
  }

  return(output)
}
