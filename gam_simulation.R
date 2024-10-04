source('semi_spp_functions.R')

# the column name of dataframe is clusters.nuisance_type.window_size.ind_or_dep.method.statistics.ppmodel
# for example poisson.linear.small.ind.logi.bias.semi

seed=1
iter=1000
print_t = 100

# pcfknown=TRUE
pcfknown=FALSE

df_result = data.frame(row.names = 1:iter)

for(window_size in c("small","large")){
  for(ind_or_dep in c("dep")){
  # for(ind_or_dep in c("ind","dep")){
    
    # generating covaiate
    is.small = (window_size=="small")
    is.independent = (ind_or_dep=="ind")
    
    data = generate_simulate_covariate(is.small = is.small, is.independent = is.independent)
    cov_target = data$cov_target 
    cov_nuisance = data$cov_nuisance
    
    # for(clusters in c("Poisson","LGCP")){
    # for(clusters in c("Poisson")){
    for(clusters in c("LGCP")){
      # for(nuisance_type in c("linear","poly")){
        for(nuisance_type in c("poly")){
        
        # generate the intensity for the point process
        alpha = 400
        theta = 0.3
        var=NULL
        scale=NULL
        formula=NULL
        
        if(nuisance_type=="linear"){
          eta = 0.3
          lambda = alpha*exp(theta*cov_target + eta*cov_nuisance)
          formula = as.formula('pp~cov_target+cov_nuisance')
        }
        if(nuisance_type=="poly"){
          eta = 0.09
          lambda = alpha*exp(theta*cov_target - eta*cov_nuisance^2)
          formula = as.formula('pp~cov_target+I(cov_nuisance^2)')
        }
        if(clusters=="LGCP"){
          var = 0.2
          scale = 0.2
          lambda = lambda*exp(-var/2)
        }
        
        # for(method in c("mpl","logi")){
        for(method in c("logi")){
          # for(ppmodel in c("semi","naive","oracle")){
          # for(ppmodel in c("semi","oracle")){
          # for(ppmodel in c("oracle")){
          # for(ppmodel in c("semi")){
          for(ppmodel in c("naive")){
            result = simulation(seed=seed,iter=iter,lambda=lambda,theta=theta,print_t=print_t,method=method,clusters=clusters,ppmodel=ppmodel,var=var,scale=scale,formula=formula,pcfknown = pcfknown)
            for(statistics in c("bias","se")){
              colname = paste(c(window_size,ind_or_dep,clusters,nuisance_type,method,ppmodel,statistics),collapse = ".")
              print(colname)
              if(statistics=="bias"){
                df_result[colname] = result$bias_lst
              }
              if(statistics=="se"){
                df_result[colname] = result$se_lst
              }
            }
          }
        }
      }
    }
  }
}





