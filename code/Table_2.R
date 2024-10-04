#### generate Table 2 in the manuscript: the simulation result of Poisson processes
#### system.time: user: 12062.444  system: 182.741   elapsed: 12258.322 
#### takes about 3.35 hours 
#### implemented on macbook pro 13.3 inch with M1 chip and 8G memory, macOS Ventura Versioin, 13.4.1(c)
#### There are 10 warnings like: In newton(lsp = lsp, X = G$X, y = G$y, Eb = G$Eb, UrS = G$UrS,  ... :Fitting terminated with step failure - check results carefully   
source('code/simulation_functions.R')

seed=1
iter=1000
print_t = 100 # report the summary every `print_t` simulations

pcfknown = FALSE
clusters = "Poisson"
method = "logi"
ppmodel = "semi"
alpha = 400
theta = 0.3
var=NULL
scale=NULL
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
      # run the simulation
      result = simulation(seed=seed,iter=iter,lambda=lambda,
                          theta=theta,print_t=print_t,method=method,
                          clusters=clusters,ppmodel=ppmodel,var=var,scale=scale,
                          formula=formula,pcfknown = pcfknown)
      # store the the bias and SE of each simulation into a Table
      for(statistics in c("bias","se")){
        # the tag of the simulaion
        colname = paste(c(window_size,ind_or_dep,nuisance_type,statistics),collapse = ".")
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

#### aggregate the simulation results into a Table

window = c()
covar = c()
nuisance = c()
bias = c()
rMSE = c()
meanSE = c()
CP90 = c()
CP95 = c()

num_setting = length(colnames(df_result))/2

for(i in 1:num_setting){
  setting_name = colnames(df_result)[2*i]
  setting_name = strsplit(setting_name, "\\.")
  
  # get setting name
  window = append(window,setting_name[[1]][1])
  covar = append(covar,setting_name[[1]][2])
  nuisance = append(nuisance,setting_name[[1]][3])
  
  bias_lst = df_result[,2*i-1]
  se_lst = df_result[,2*i]
  
  # get summary statistics for the simulation setting
  bias = append(bias,mean(bias_lst))
  rMSE = append(rMSE,sqrt(mean(bias_lst^2)))
  meanSE = append(meanSE,mean(se_lst))
  CP90 = append(CP90,sum(abs(bias_lst)<=qnorm(0.95)*se_lst)/length(se_lst))
  CP95 = append(CP95,sum(abs(bias_lst)<=qnorm(0.975)*se_lst)/length(se_lst))
}

# adjust the scale of the summary
CP90 = CP90 * 100
CP95 = CP95 * 100
bias = bias * 100

# output and save the table
Table_2 = data.frame(Window = window,Covar = covar,Nuisance = nuisance,Bias = bias, rMSE = rMSE,meanSE = meanSE,CP90=CP90,CP95=CP95)
write.csv(Table_2, "output/Table_2.csv", row.names=FALSE)

