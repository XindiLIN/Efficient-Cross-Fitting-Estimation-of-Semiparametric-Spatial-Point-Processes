#### system.time:  user: 32.075  system: 47.729   elapsed: 349.628 
#### takes 5.82 minutes
#### implemented on macbook pro 13.3 inch with M1 chip and 8G memory, macOS Ventura Versioin, 13.4.1(c)

source('code/semi_spp_functions.R')
set.seed(1)

#### Beilschmiedia pendula

pp = bei
target = "elev"
nuisance = "grad"
data = bei.extra
method = "logi"
clusters = "Thomas"

## semiparametric model
bei.semi = semi_ppm_gam(pp=bei,target = "elev",nuisance = "grad",data = bei.extra,method = "logi",clusters = "Thomas")

bei.est.semi = bei.semi$theta_hat
bei.se.semi = bei.semi$se_hat

# adjust the scale
bei.est.semi = bei.est.semi * 100
bei.se.semi = bei.se.semi * 100

bei.lower.semi = bei.est.semi - qnorm(0.975)*bei.se.semi
bei.upper.semi = bei.est.semi + qnorm(0.975)*bei.se.semi

## parametric model
bei.para = kppm(bei~elev+grad,data = bei.extra,clusters = "Thomas")

bei.est.para = coef(bei.para)["elev"]
bei.se.para = sqrt(diag(vcov(bei.para))["elev"])

# adjust the scale
bei.est.para = bei.est.para * 100
bei.se.para = bei.se.para * 100

bei.lower.para = bei.est.para - qnorm(0.975)*bei.se.para
bei.upper.para = bei.est.para + qnorm(0.975)*bei.se.para

## output the results
bei.result = data.frame(Estimate = c(bei.est.para,bei.est.semi), 
           SE = c(bei.se.para,bei.se.semi), 
           CI95.lower = c(bei.lower.para,bei.lower.semi),
           CI95.upper = c(bei.upper.para,bei.upper.semi),row.names = c("para","semi"))

write.csv(bei.result, "output/Beilschmiedia.csv")



#### Capparis frondosa

## load the Capparis point pattern

load("data/bci.tree1.rdata")
cappfr = bci.tree1[bci.tree1$sp=="cappfr",]
cappfr = cappfr[!duplicated(cappfr[,c("gx","gy")]),]
pp_cap = ppp(cappfr$gx,cappfr$gy,window =owin(c(0,1000),c(0,500)))

## semiparametric model
cap.semi = semi_ppm_gam(pp=pp_cap,target = "elev",nuisance = "grad",data = bei.extra,method = "logi",clusters = "Thomas")

cap.est.semi = cap.semi$theta_hat
cap.se.semi = cap.semi$se_hat

# adjust the scale

cap.est.semi = cap.est.semi * 100
cap.se.semi = cap.se.semi * 100

cap.lower.semi = cap.est.semi - qnorm(0.975)*cap.se.semi
cap.upper.semi = cap.est.semi + qnorm(0.975)*cap.se.semi

## parametric model
cap.para = kppm(pp_cap~elev+grad,data = bei.extra,clusters = "Thomas")

cap.est.para = coef(cap.para)["elev"]
cap.se.para = sqrt(diag(vcov(cap.para))["elev"])

# adjust the scale

cap.est.para = cap.est.para * 100
cap.se.para = cap.se.para * 100

cap.lower.para = cap.est.para - qnorm(0.975)*cap.se.para
cap.upper.para = cap.est.para + qnorm(0.975)*cap.se.para

## adjust the scale

## output the results
cap.result = data.frame(Estimate = c(cap.est.para,cap.est.semi), 
                        SE = c(cap.se.para,cap.se.semi), 
                        CI95.lower = c(cap.lower.para,cap.lower.semi),
                        CI95.upper = c(cap.upper.para,cap.upper.semi),row.names = c("para","semi"))

write.csv(cap.result, "output/Capparis.csv")

