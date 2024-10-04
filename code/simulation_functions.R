source('semi_spp_functions.R')



#### Description: read in the generated gaussian process as target and nuisance covariates
generate_simulate_covariate = function(field_1 = NULL,field_2 = NULL,is.small=TRUE,is.independent=TRUE){
  if(is.null(field_1)& is.null(field_2)){
    simulated_field_1_large = read.csv('simulated_field_1_large.csv')
    simulated_field_1_large$X = NULL
    simulated_field_1_large = as.numeric(simulated_field_1_large)
    simulated_field_1_large = matrix(simulated_field_1_large,nrow = 100,ncol = 100)
    
    simulated_field_2_large = read.csv('simulated_field_2_large.csv')
    simulated_field_2_large$X = NULL
    simulated_field_2_large = as.numeric(simulated_field_2_large)
    simulated_field_2_large = matrix(simulated_field_2_large,nrow = 100,ncol = 100)
  }
  else{
    simulated_field_1_large = field_1
    simulated_field_2_large = field_2
  }
  
  
  simulated_field_1_small = simulated_field_1_large[1:50,1:50]
  simulated_field_2_small = simulated_field_2_large[1:50,1:50]
  
  if(is.small &is.independent){
    cov_target = as.im(simulated_field_1_small,unit.square())
    cov_nuisance = as.im(simulated_field_2_small,unit.square())
  }
  
  if(is.small & (!is.independent)){
    cov_target = as.im(simulated_field_1_small,unit.square())
    cov_nuisance = as.im(simulated_field_2_small*simulated_field_1_small,unit.square())
  }
  
  if((!is.small) &is.independent){
    cov_target = as.im(simulated_field_1_large[1:100,1:100],W=owin(c(0,2),c(0,2)))
    cov_nuisance = as.im(simulated_field_2_large[1:100,1:100],W=owin(c(0,2),c(0,2)))
  }
  
  if((!is.small) &(!is.independent)){
    cov_target = as.im(simulated_field_1_large,W=owin(c(0,2),c(0,2)))
    cov_nuisance = as.im(simulated_field_2_large*simulated_field_1_large,W=owin(c(0,2),c(0,2)))
  }
  
  return(as.solist(list("cov_target"=cov_target,"cov_nuisance"=cov_nuisance)))
}


#### Description
## run the simulation

#### Values
## a list containing the bias (bias_lst) and standard error (se_lst) of each every simulation

#### Arguments
## seed: integer; random seed
## iter: integer; number of repetitions of simulations
## lambda: `pixel image` in spatstat; intensity function of the simulated spatial point process
## theta: real-valued number; true value of the target parameter
## print_t: integer; print the current simulation summary every `print_t` times of simulations
## clusters: `LGCP` or `Poisson`; type of spatial point process we want to generate
## ppmodel: `semi` or `oracle` or `naive`; the model we fit for the simulated process
## var: real-valued number; variance parameter of the PCF 
## scale: real-valued number; scale parameter of the PCF 
## formula: formula; the oracle model of the intensity function
## pcfknown: TRUE or FALSE; indicator of whether the PCF is known in the simulation setting

simulation = function(seed=1,iter=1000,lambda,theta,
                      print_t=100,method,clusters,ppmodel="semi",
                      var=NULL,scale=NULL,formula=NULL,pcfknown=FALSE){
  set.seed(seed = seed)
  theta_hat_lst = c()
  se_lst = c()
  for(i in 1:iter){
    if(clusters=="Poisson"){
      pp = rpoispp(lambda)
      modelpcf = NULL
    }
    else{
      pp_LGCP = rLGCP("exponential",0,var=var, scale=scale,win=owin(xrange = lambda$xrange,yrange = lambda$yrange))
      lambda_LGCP = as.im(attr(pp_LGCP, "Lambda"),dimyx=lambda$dim[1])*lambda
      pp = rpoispp(lambda_LGCP)
      if(pcfknown==TRUE){
        modelpcf = function(r){exp(var*exp(-r/scale))}
      }
      else{
        modelpcf = NULL
      }
    }
    if(ppmodel=='semi'){
      semi_ppm_gam_fit = semi_ppm_gam(pp=pp,target="cov_target",nuisance="cov_nuisance",data=data,method=method,clusters=clusters,modelpcf=modelpcf)
      theta_hat_lst = append(theta_hat_lst,semi_ppm_gam_fit$theta_hat)
      se_lst = append(se_lst,semi_ppm_gam_fit$se_hat[1])
    }
    if(ppmodel=="naive"){
      
      if(clusters=="Poisson"){
        fit = ppm(pp~cov_target+cov_nuisance,method=method)
        s = summary(fit)
        theta_hat = s$coefs.SE.CI['cov_target','Estimate']
        se = s$coefs.SE.CI['cov_target','S.E.']
        
        theta_hat_lst = append(theta_hat_lst,theta_hat)
        se_lst = append(se_lst,se)
      }
      if(clusters=="LGCP"){
        fit = kppm(pp~cov_target+cov_nuisance,clusters = clusters)
        s = summary(fit)
        theta_hat = s$coefs.SE.CI['cov_target','Estimate']
        se = s$coefs.SE.CI['cov_target','S.E.']
        
        theta_hat_lst = append(theta_hat_lst,theta_hat)
        se_lst = append(se_lst,se)
      }
    }
    if(ppmodel=="oracle"){
      if(clusters=="Poisson"){
        fit = ppm(formula,method=method)
        s = summary(fit)
        theta_hat = s$coefs.SE.CI['cov_target','Estimate']
        se = s$coefs.SE.CI['cov_target','S.E.']
        
        theta_hat_lst = append(theta_hat_lst,theta_hat)
        se_lst = append(se_lst,se)
      }
      if(clusters=="LGCP"){
        fit = kppm(formula,clusters = clusters)
        s = summary(fit)
        theta_hat = s$coefs.SE.CI['cov_target','Estimate']
        se = s$coefs.SE.CI['cov_target','S.E.']
        
        theta_hat_lst = append(theta_hat_lst,theta_hat)
        se_lst = append(se_lst,se)
      }
    }
    if(i%%print_t==0){
      cat(i)
      cat(' ')
      cat('\n')
      print_result(bias_lst = theta_hat_lst-theta,se_lst)
    }
  }
  return(list('bias_lst' = theta_hat_lst-theta,'se_lst'=se_lst))
  
}



