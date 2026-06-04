library(spatstat)
library(mgcv)

load("output/PFAS/prepared_data.RData")
load("output/PFAS/cov_extrapolation_im.RData")

covariate_names <- c("DistanceToAirport",
                     "DistanceToMilitaryBase",
                     "DistanceToLandfill",
                     "DistanceToApparel",
                     "DistanceToFurniture",
                     "PopDen",
                     "Median.Earnings.Inflation.Adj.2020",
                     "Bachelor.Degree.Percent.2020")

parametric_formula = as.formula(paste("pp_PFAS_PA ~", paste(covariate_names, collapse = "+")))

semiparametric_formula = as.formula(paste0(
  "pp_PFAS_PA ~",
  paste(covariate_names[1:5], collapse = "+"), "+",
  paste("s(", covariate_names[-(1:5)], ")", collapse = "+")
))

ppmfit_para = ppm(parametric_formula, data = covariate_images, use.gam = TRUE, method = "mpl")
ppmfit_semi = ppm(semiparametric_formula, data = covariate_images, use.gam = TRUE, method = "mpl")

gamfit_semi <- ppmfit_semi$internal$glmfit
gamfit_para <- ppmfit_para$internal$glmfit
coef_semi   <- ppmfit_semi$coef
coef_para   <- ppmfit_para$coef
vcov_para   <- vcov(ppmfit_para)
wt_semi     <- ppmfit_semi$Q$w
