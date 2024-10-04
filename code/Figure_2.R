##### plot logi nuisance estimation plot
#### system times: user: 6.405 seconds, system: 0.207, elapsed: 6.617 
#### takes about 6.617 seconds to run 
#### implemented on macbook pro 13.3 inch with M1 chip and 8G memory, macOS Ventura Versioin, 13.4.1(c)

library("dplyr")
library("gratia")
library(ggplot2)
library(spatstat)

# load data of capparis
load("bci.tree1.rdata")
cappfr = bci.tree1[bci.tree1$sp=="cappfr",]
cappfr = cappfr[!duplicated(cappfr[,c("gx","gy")]),]
pp_cap = ppp(cappfr$gx,cappfr$gy,window =owin(c(0,1000),c(0,500)))

ppmfit_bei =  ppm(bei~elev+s(grad), data = bei.extra,use.gam = TRUE, method="mpl")
ppmfit_cap =  ppm(pp_cap~elev+s(grad), covariates = bei.extra,use.gam = TRUE, method="mpl")

gamfit_bei = ppmfit_bei$internal$glmfit
gamfit_cap = ppmfit_cap$internal$glmfit

newdata = data.frame(grad=seq(0,0.3,length.out=100),elev=seq(0,0.3,length.out=100))

fitted_bei = data.frame(grad=seq(0,0.3,length.out=100),.fitted=predict(gamfit_bei,newdata = newdata,type="terms",exclude = "elev"))
fitted_cap = data.frame(grad=seq(0,0.3,length.out=100),.fitted=predict(gamfit_cap,newdata = newdata,type="terms",exclude = "elev"))

fitted_bind = fitted_cap %>% bind_rows(fitted_bei) %>% mutate(species=c(rep("Capparis ",dim(fitted_cap)[1]),rep("Beilschmiedia",dim(fitted_bei)[1])))

png("output/Figure_2.png", width = 800, height = 600)

ggplot(data = fitted_bind,aes(x=grad,y=s.grad.,group=species))+geom_line(aes(colour=species,linetype=species)) +labs(y="Value of the Estimated Nuisance Function")

dev.off()
