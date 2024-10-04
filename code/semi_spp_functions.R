library(spatstat)
library(mgcv)
library(stats)

#### Description
## statistical inference of the target parameter of a semiparametric spatial point pattern

#### Values
## a lists containing the point estimate and the standard error of the target parameter, `theta_hat` is the point estimate, `se_hat` is the standard error

#### Arguments
## pp: a `ppp` object in spatstat; the observed point pattern
## target: string or vector of strings; name of the target covariates
## nuisance: string or vector of strings; name of the nuisance covariates
## data： list of piexel image in spatstat; list of covariates
## method: 'mpl' or 'logi'; 'mpl' indicates using quadrature approximation, 'logi' indicates using logistic approximation
## clusters: 'Poisson', 'LGCP', 'Thomas' etc; indicates the observed point pattern is assumed to be Poisson or LGCP
## modelpcf: function; when it is NULL, the pcf will be estimated, when it is not NULL, it will be used as the pcf of the spatial point process 

semi_ppm_gam = function(pp,target="cov_target",nuisance="cov_nuisance",data,method="mpl",clusters="Poisson",modelpcf=NULL){
  # define the semiparametric trend for gam
  target_trend = paste(target, collapse = "+")
  nuisance_trend = paste(nuisance,collapse = ",")
  nuisance_trend = paste("s(",nuisance_trend,")")
  trend = paste("~",target_trend,"+",nuisance_trend)
  
  if(method=="mpl"){
    # ppmfit is the linear model ppm fitted with either "logi" method or "mpl" method
    ppmfit =  ppm(Q = pp,trend = as.formula(trend), covariates = data, use.gam = TRUE, method="mpl")
    # gamfit is the semiparametric gam fitted model
    gamfit = ppmfit$internal$glmfit
    theta_hat = coef(ppmfit)[target]
    se_hat = sqrt(vcov.semi_ppm_gam(gamfit = gamfit,ppmfit = ppmfit,target=target,nuisance=nuisance,clusters = clusters,method=method,modelpcf=modelpcf))
    return(list("theta_hat"=theta_hat,"se_hat"=se_hat))
  }
  if(method=="logi"){
    trend_linear = as.formula(paste("~",paste(c(target,nuisance),collapse="+")))
    # ppmfit is the linear model ppm fitted with either "logi" method or "mpl" method
    ppmfit = ppm(Q = pp,trend = trend_linear ,covariates = data, method="logi")
    glmdata = ppmfit$internal$glmdata
    # gamfit is the semiparametric gam fitted model
    gamfit = gam(as.formula(paste('.logi.Y',trend)),data = glmdata, family = binomial(), subset = .logi.ok, weights = .logi.w)
    theta_hat = coef(gamfit)[target]
    se_hat = sqrt(vcov.semi_ppm_gam(gamfit = gamfit,ppmfit = ppmfit,target=target,nuisance=nuisance,clusters = clusters,method=method,modelpcf=modelpcf))
    return(list("theta_hat"=theta_hat,"se_hat"=se_hat))
  }
}


#### Description
## returns the estimated asymptotic variance of the point estimates of the target parameter

#### Values
## a real-valued number 

#### Arguments
## ppmfit: a linear fitted model returned by `ppm` in spatstat;
## target: string or vector of strings; name of the target covariates
## nuisance: string or vector of strings; name of the nuisance covariates
## method: 'mpl' or 'logi'; 'mpl' indicates using quadrature approximation, 'logi' indicates using logistic approximation
## what: string; objective function used to estimate the pcf. Default is `K`, indicates using K function to estimate pcf
## clusters: 'Poisson', 'LGCP', 'Thomas' etc; indicates the observed point pattern is assumed to be Poisson or LGCP
## modelpcf: function; when it is NULL, the pcf will be estimated, when it is not NULL, it will be used as the pcf of the spatial point process, and the function will not estimate the pcf 

vcov.semi_ppm_gam = function(gamfit,ppmfit,target="cov_target",nuisance="cov_nuisance",method="mpl",clusters="Poisson",what = "K",modelpcf=NULL){
  if(clusters=="Poisson"){
    if(method=="mpl"){
      return(vcov.mpl(gamfit,ppmfit))
    }
    if(method=="logi"){
      return(vcov.logi(gamfit,ppmfit))
    }
  }
  if(clusters!="Poisson"){
    return(vcov.lgcp(gamfit,ppmfit,target=target,nuisance=nuisance,clusters=clusters,what=what,modelpcf=modelpcf))
  }
}

#### Description
## Though named as lgcp, it deal with any non-Poisson processes
## Returns the estimated asymptotic variance of the point estimates of the target parameter when the observed point pattern is not Poisson

#### Values
## a real-valued number 

#### Arguments
## ppmfit: a linear fitted model returned by `ppm` in spatstat;
## target: string or vector of strings; name of the target covariates
## nuisance: string or vector of strings; name of the nuisance covariates
## clusters: 'Poisson', 'LGCP', 'Thomas' etc; indicates the observed point pattern is assumed to be Poisson or LGCP
## what: string; objective function used to estimate the pcf. Default is `K`, indicates using K function to estimate pcf
## modelpcf: function; when it is NULL, the pcf will be estimated, when it is not NULL, it will be used as the pcf of the spatial point process, and the function will not estimate the pcf 


vcov.lgcp = function(gamfit,ppmfit,target,nuisance,clusters = "LGCP",what = "K",modelpcf=NULL){
  cov.df = gamfit$model
  lambda = GAMpredict(gamfit = gamfit,ppmfit = ppmfit)
  wt = ppmfit$Q$w

  weighted_loess <- loess(as.formula(paste(target,"~",nuisance)),weights = lambda*wt, data = cov.df)
  lfd = predict(weighted_loess,newdata = cov.df)
  lfd_projection = cov.df[,target]-lfd
  
  sensitivity= t(lfd_projection) %*% (lfd_projection * lambda * wt)
  
  
  if(is.null(modelpcf)){
    # estimate pcf by minicontrast
    modelpcf = Kpcf.semi_ppm_gam(gamfit,ppmfit,what=what,clusters = clusters)
  }
  
  # weight matrix
  U = union.quad(ppmfit$Q)
  weight_mat = matrix(modelpcf(c(pairdist(U))) - 1, npoints(U), npoints(U))
  pseudo_info  = t(lfd_projection * lambda * wt) %*% weight_mat %*% (lfd_projection * lambda * wt)
  return ((sensitivity+pseudo_info)/sensitivity^2)
}

#### Description
## returns the estimated asymptotic variance of target estimation when the observed point pattern is LGCP, the method is `mpl` 

#### Values
## a real-valued number 

#### Arguments
## ppmfit: a linear fitted model returned by `ppm` in spatstat;
## target: string or vector of strings; name of the target covariates

vcov.mpl = function(gamfit,ppmfit){
  cov.df = gamfit$model
  lambda = GAMpredict(gamfit = gamfit,ppmfit = ppmfit)
  wt = ppmfit$Q$w
  
  weighted_loess <- loess(cov_target~cov_nuisance,weights = lambda*wt, data = cov.df)
  lfd = predict(weighted_loess,newdata = cov.df)
  lfd_projection = cov.df$cov_target-lfd
  
  sensitivity= t(lfd_projection) %*% (lfd_projection * lambda * wt)
  return(1/sensitivity)
}

#### Description
## returns the estimated asymptotic variance of target estimation when the observed point pattern is LGCP, the method is `mpl` 

#### Values
## a real-valued number 

#### Arguments
## ppmfit: a linear fitted model returned by `ppm` in spatstat;
## target: string or vector of strings; name of the target covariates
## how: "stratrand"; the type of the dummy process generated to do logistic approximation

vcov.logi = function(gamfit,ppmfit,how="stratrand"){
  # get projection
  cov.df = gamfit$model
  lambda = GAMpredict(gamfit = gamfit,ppmfit = ppmfit)
  wt = ppmfit$Q$w
  weighted_loess <- loess(cov_target~cov_nuisance,weights = lambda*wt, data = cov.df)
  lfd = predict(weighted_loess,newdata = cov.df)
  lfd_projection = cov.df$cov_target-lfd
  
  # get A1 and Slog
  rho = ppmfit$internal$logistic$rho
  Slog = t(lfd_projection) %*% (lfd_projection * lambda*rho/(lambda+rho)^2)
  A1log = t(lfd_projection) %*% (lfd_projection * lambda*rho*rho/(lambda+rho)^3)
  Sigma1log = A1log
  
  # fitted again and get another fitted model result
  Q2 <- quadscheme.logi(data=ppmfit$Q$data, dummytype = "stratrand",nd = ppmfit$internal$logistic$nd)
  ppmfit2 = ppm(Q = Q2,trend = ppmfit$trend ,covariates = ppmfit$covariates, method="logi")
  glmdata2 = ppmfit2$internal$glmdata
  gamfit2 = gam(gamfit$formula,data = glmdata2, family = binomial(), subset = .logi.ok, weights = .logi.w)
  
  # new fitted intensity
  lambda2 = GAMpredict(gamfit = gamfit2,ppmfit = ppmfit2)
  # new model projection
  cov.df2 = gamfit2$model
  wt2 = ppmfit2$Q$w
  weighted_loess2 <- loess(cov_target~cov_nuisance,weights = lambda2*wt2, data = cov.df2)
  lfd2 = predict(weighted_loess2,newdata = cov.df2)
  lfd_projection2 = cov.df2$cov_target-lfd2
  
  lamdum <- lambda[!is.data(ppmfit$Q)]
  lamdum2 <- lambda2[!is.data(Q2)]
  
  mdum <- lfd_projection[!is.data(ppmfit$Q)]
  mdum2 <- lfd_projection2[!is.data(Q2)]
  
  wlam <- mdum * rho*lamdum/(lamdum+rho)
  wlam2 <- mdum2 * rho*lamdum2/(lamdum2+rho)
  Sigma2log <- t(wlam-wlam2)%*%(wlam-wlam2)/(2*rho*rho)
  
  sensitivity = Slog*(1/(Sigma1log+Sigma2log))*Slog
  
  return(1/sensitivity)
}

# the gam fitted model either by "logi" or "mpl"
# the returned value is a vector of fitted intensity for each point in quadrature point
GAMpredict = function(gamfit,ppmfit){
  mm = model.matrix(gamfit)
  coefs = gamfit$coefficients
  eta = as.vector(mm %*% coefs)
  linkinv = family(gamfit)$linkinv
  answer <- linkinv(eta)
  
  if(family(gamfit)$family=="binomial"){
    answer <- ppmfit$internal$glmdata$.logi.B[1] * answer/(1-answer)
  }
  return(answer)  
}


# return the estimated pcf function
Kpcf.semi_ppm_gam <- function(gamfit,ppmfit, what=c("K", "pcf", "kernel"),clusters = c("LGCP","Thomas")) {
  what <- match.arg(what)
  clusters <- match.arg(clusters)
  
  w <- as.owin(ppmfit, from="covariates")
  
  
  lambda = GAMpredict(gamfit,ppmfit) # the fitted intensity is not at the points of quadrature, but at the 
  lambda = lambda[1:ppmfit$Q$data$n] # the first pp$n is the fitted intensity at ppp
  model = clusterfit(X=ppmfit$Q$data,clusters = clusters,lambda = lambda) # mini contrast
  
  # Extract function definition from internal table
  tableentry <- spatstatClusterModelInfo(clusters)
  fun <- tableentry[["pcf"]]
  
  # Extract model parameters
  par <- model$par
  # Extract covariance model (if applicable)
  cm <- model$covmodel
  model <- cm$model
  margs <- cm$margs
  #
  f <- function(r) as.numeric(fun(par=par, rvals=r,
                                  model=model, margs=margs))
  return(f)
}


print_result = function(bias_lst,se_lst){
  cat('Bias: ')
  cat(mean(bias_lst))
  cat('\n')
  cat('rMSE: ')
  cat(sqrt(mean(bias_lst^2)))
  cat('\n')
  # cat('mSE: ')
  cat('medianSE: ')
  # cat(mean(se_lst))
  cat(median(se_lst))
  cat('\n')
  cat('CP90: ')
  cat(sum(abs(bias_lst)<=qnorm(0.95)*se_lst)/length(se_lst))
  cat('\n')
  cat('CP95: ')
  cat(sum(abs(bias_lst)<=qnorm(0.975)*se_lst)/length(se_lst))
  cat('\n')
}

