#### generate Table 4 in the manuscript: the simulation result of model comparison
#### system.time() user: 11179.098    system: 1183.815   elapsed: 12647.875 
#### takes about 3.105 hours to run
#### implemented on macbook pro 13.3 inch with M1 chip and 8G memory, macOS Ventura Versioin, 13.4.1(c)
   
source('simulation_functions.R')

seed=1
iter=1000
print_t = 100
ind_or_dep = "dep"
nuisance_type = "poly"
method = "logi"
alpha = 400
theta = 0.3
pcfknown = FALSE
var = 0.2
scale = 0.2
formula = as.formula('pp~cov_target+I(cov_nuisance^2)')

df_result = data.frame(row.names = 1:iter)


#### run `iter`-times simulations
for(window_size in c("small","large")){
  # generating covaiate
  is.small = (window_size=="small")
  is.independent = (ind_or_dep=="ind")
  
  data = generate_simulate_covariate(is.small = is.small, is.independent = is.independent)
  cov_target = data$cov_target 
  cov_nuisance = data$cov_nuisance
  # generate intensity function
  
  for(clusters in c("Poisson","LGCP")){
    eta = 0.09
    lambda = alpha*exp(theta*cov_target - eta*cov_nuisance^2)
    if(clusters=="LGCP"){lambda = lambda*exp(-var/2)}
    for(ppmodel in c("semi","naive","oracle")){
      # run the simulation
      result = simulation(seed=seed,iter=iter,lambda=lambda,
                          theta=theta,print_t=print_t,method=method,
                          clusters=clusters,ppmodel=ppmodel,var=var,scale=scale,
                          formula=formula,pcfknown = pcfknown)
      # store the the bias and SE of each simulation into a Table
      for(statistics in c("bias","se")){
        # the tag of the simulaion
        colname = paste(c(window_size,clusters,ppmodel,statistics),collapse = ".")
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
process = c()
model = c()
bias = c()
rMSE = c()
meanSE = c()
CP90 = c()
CP95 = c()

num_setting = length(colnames(df_result))/2

for(i in 1:num_setting){
  bias_lst = df_result[,2*i-1]
  se_lst = df_result[,2*i]
  setting_name = colnames(df_result)[2*i]
  setting_name = strsplit(setting_name, "\\.")
  # get setting name
  window = append(window,setting_name[[1]][1])
  process = append(process,setting_name[[1]][2])
  model = append(model,setting_name[[1]][3])
  # get summary statistics for the simulation setting
  bias = append(bias,mean(bias_lst))
  rMSE = append(rMSE,sqrt(mean(bias_lst^2)))
  meanSE = append(meanSE,mean(se_lst))
  CP90 = append(CP90,sum(abs(bias_lst)<=qnorm(0.95)*se_lst)/length(se_lst))
  CP95 = append(CP95,sum(abs(bias_lst)<=qnorm(0.975)*se_lst)/length(se_lst))
}

model[model=="naive"] = "para"

# adjust the scale of the summary
CP90 = CP90 * 100
CP95 = CP95 * 100
bias = bias * 100

# output and save the table
Table_4 = data.frame(Window = window,Process = process,Model=model,Bias = bias, rMSE = rMSE,
                     meanSE = meanSE,CP90=CP90,CP95=CP95)
write.csv(Table_4, "Table_4.csv", row.names=FALSE)

