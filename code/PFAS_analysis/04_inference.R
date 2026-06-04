library(mgcv)
library(ggplot2)
library(gratia)
source('code/semi_spp_functions.R')

load("output/PFAS/fitted_models.RData")

target_covariate_names   = covariate_names[1:5]
nuisance_covariate_names = covariate_names[6:8]

cov.df  = gamfit_semi$model
lambda  = GAMpredict(gamfit = gamfit_semi, ppmfit = NULL)
wt      = wt_semi

# Partial-out nuisance covariates via weighted loess
lfd_projection = list()
for (target_cov in target_covariate_names) {
  formula <- as.formula(paste0(target_cov, "~", paste(nuisance_covariate_names, collapse = "+")))
  weighted_loess <- loess(formula, weights = lambda * wt, data = cov.df)
  lfd_projection[[target_cov]] = predict(weighted_loess, newdata = cov.df)
}

lfd_projection = as.data.frame(lfd_projection)
lfd_projection = cov.df[, target_covariate_names] - lfd_projection
lfd_projection = as.matrix(lfd_projection)

sensitivity = t(lfd_projection) %*% (lfd_projection * lambda * wt)
semi_se  <- sqrt(diag(solve(sensitivity)))

##### build inference results table
semi_est <- coef_semi[target_covariate_names]
para_est <- coef_para[target_covariate_names]
para_se  <- sqrt(diag(vcov_para)[target_covariate_names])

results_df <- data.frame(
  Coefficient    = target_covariate_names,
  Semi_Estimate  = round(semi_est, 4),
  Semi_SE        = round(semi_se,  4),
  Para_Estimate  = round(para_est, 4),
  Para_SE        = round(para_se,  4),
  row.names      = NULL
)

print(results_df)
write.csv(results_df, "output/inference_results.csv", row.names = FALSE)

###### plot nuisance nonlinear effect
draw(gamfit_semi, select = "s(Median.Earnings.Inflation.Adj.2020)", rug = FALSE, ci_alpha = 0) +
  geom_line(color = "#0479A8", linewidth = 1.0) +
  coord_cartesian(xlim = c(15000, 66000)) +
  labs(
    x = "Median Earnings Per Year ($)",
    y = "Estimated Nuisance Function"
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank())
