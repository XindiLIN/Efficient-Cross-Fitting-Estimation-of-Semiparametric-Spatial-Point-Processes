#### generate Table 3 in the manuscript: the simulation result of log-Gaussian Cox processes
#### system times: user: 28392.074 seconds, system: 4920.887, elapsed: 35749.550 
#### takes about 9.93 hours to run 
#### implemented on macbook pro 13.3 inch with M1 chip and 8G memory, macOS Ventura Versioin, 13.4.1(c)
#### There are 13 warnings like: In newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS,  ... :Fitting terminated with step failure - check results carefully
source('code/simulation_functions.R')

seed=1
iter=1000
print_t = 100
clusters = "LGCP"
method = "logi"
ppmodel = "semi"
alpha = 400
theta = 0.3
var = 0.2
scale = 0.2
formula=NULL

df_result = data.frame(row.names = 1:iter)


#### run `iter`-times simulations
for(window_size in c("small","large")){
  for(ind_or_dep in c("ind","dep")){
    # generating covaiate
    is.small = (window_size=="small")
    is.independent = (ind_or_dep=="ind")
    data = generate_simulate_covariate(is.small = is.small, is.independent = is.independent)
    cov_target = data$cov_target 
    cov_nuisance = data$cov_nuisance
    
    for(nuisance_type in c("linear","poly")){
      # generate the intensity for the point process
      if(nuisance_type=="linear"){
        eta = 0.3
        lambda = alpha*exp(theta*cov_target + eta*cov_nuisance)
      }
      if(nuisance_type=="poly"){
        eta = 0.09
        lambda = alpha*exp(theta*cov_target - eta*cov_nuisance^2)
      }
      lambda = lambda*exp(-var/2)
      for(pcfknown in c(TRUE,FALSE)){
        # run the simulation
        result = simulation(seed=seed,iter=iter,lambda=lambda,
                            theta=theta,print_t=print_t,method=method,
                            clusters=clusters,ppmodel=ppmodel,var=var,scale=scale,
                            formula=formula,pcfknown = pcfknown)
        # store the the bias and SE of each simulation into a Table
        for(statistics in c("bias","se")){
          # the tag of the simulaion
          colname = paste(c(window_size,ind_or_dep,nuisance_type,pcfknown,statistics),collapse = ".")
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

#### aggregate the simulation results into a Table

window = c()
covar = c()
nuisance = c()
bias = c()
rMSE = c()
meanSE = c()
meanSE_knownpcf = c()
CP90 = c()
CP90_knownpcf = c()
CP95 = c()
CP95_knownpcf = c()

num_setting = length(colnames(df_result))/2

for(i in 1:num_setting){
  bias_lst = df_result[,2*i-1]
  se_lst = df_result[,2*i]
  setting_name = colnames(df_result)[2*i]
  setting_name = strsplit(setting_name, "\\.")
  
  PCFknown = setting_name[[1]][4]
  
  if(PCFknown == 'FALSE'){
    # get setting name
    window = append(window,setting_name[[1]][1])
    covar = append(covar,setting_name[[1]][2])
    nuisance = append(nuisance,setting_name[[1]][3])
    # get summary statistics for the simulation setting
    bias = append(bias,mean(bias_lst))
    rMSE = append(rMSE,sqrt(mean(bias_lst^2)))
    meanSE = append(meanSE,mean(se_lst))
    CP90 = append(CP90,sum(abs(bias_lst)<=qnorm(0.95)*se_lst)/length(se_lst))
    CP95 = append(CP95,sum(abs(bias_lst)<=qnorm(0.975)*se_lst)/length(se_lst))
  }
  if(PCFknown == 'TRUE'){
    meanSE_knownpcf = append(meanSE_knownpcf,mean(se_lst))
    CP90_knownpcf = append(CP90_knownpcf,sum(abs(bias_lst)<=qnorm(0.95)*se_lst)/length(se_lst))
    CP95_knownpcf= append(CP95_knownpcf,sum(abs(bias_lst)<=qnorm(0.975)*se_lst)/length(se_lst))
    
  }
  
}

# adjust the scale of the summary
CP90 = CP90 * 100
CP95 = CP95 * 100
CP90_knownpcf = CP90_knownpcf * 100
CP95_knownpcf = CP95_knownpcf * 100
bias = bias * 100


# output and save the table
Table_3 = data.frame(Window = window,Covar = covar,Nuisance = nuisance,Bias = bias, rMSE = rMSE,
                     meanSE = meanSE,meanSE_star = meanSE_knownpcf,
                     CP90=CP90,CP90_star = CP90_knownpcf,
                     CP95=CP95,CP95_star = CP95_knownpcf)
write.csv(Table_3, "output/Table_3.csv", row.names=FALSE)

