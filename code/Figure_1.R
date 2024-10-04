#### system.time: user: 5.776  system: 0.171   elapsed: 5.958 
#### takes 4.723 seconds
#### implemented on macbook pro 13.3 inch with M1 chip and 8G memory, macOS Ventura Versioin, 13.4.1(c)
source('code/semi_spp_functions.R')
set.seed(1)

# load data of capparis
load("data/bci.tree1.rdata")
cappfr = bci.tree1[bci.tree1$sp=="cappfr",]
cappfr = cappfr[!duplicated(cappfr[,c("gx","gy")]),]
pp_cap = ppp(cappfr$gx,cappfr$gy,window =owin(c(0,1000),c(0,500)))


ppmfit.bei.semi =  ppm(bei~elev+s(grad), data = bei.extra,use.gam = TRUE, method="mpl")
intensity.bei.semi = predict(ppmfit.bei.semi)

ppmfit.bei.para =  ppm(bei~elev+grad, data = bei.extra, method="mpl")
intensity.bei.para = predict(ppmfit.bei.para)

common.min.bei = min(min(intensity.bei.semi$v),min(intensity.bei.para$v))
common.max.bei = max(max(intensity.bei.semi$v),max(intensity.bei.para$v))


ppmfit.cap.semi =  ppm(pp_cap~elev+s(grad), covariates = bei.extra,use.gam = TRUE, method="mpl")
intensity.cap.semi = predict(ppmfit.cap.semi)

ppmfit.cap.para =  ppm(pp_cap~elev+grad, data = bei.extra, method="mpl")
intensity.cap.para = predict(ppmfit.cap.para)

common.min.cap = min(min(intensity.cap.semi$v),min(intensity.cap.para$v))
common.max.cap = max(max(intensity.cap.semi$v),max(intensity.cap.para$v))

png("output/Figure_1.png", width = 1200, height = 600)

par(mfrow = c(2, 2))

# bei
plot(intensity.bei.para,main="Parametric Estimate of Intensity Function",
     zlim=c(common.min.bei,common.max.bei),ribbon = TRUE)
plot(bei,cols = "white", cex = .2, pch = 16,add=TRUE)
mtext("Beilschmiedia", side = 2, line = 1, cex = 1.2)

plot(intensity.bei.semi,main="Semiparametric Estimate of Intensity Function",
     zlim=c(common.min.bei,common.max.bei),ribbon = TRUE)
plot(bei,cols = "white", cex = .2, pch = 16,add=TRUE)

# cap
plot(intensity.cap.para,zlim=c(common.min.cap,common.max.cap),ribbon = TRUE,main="")
plot(pp_cap,cols = "white", cex = .2, pch = 16,add=TRUE)
mtext("Capparis", side = 2, line = 1, cex = 1.2)

plot(intensity.cap.semi,zlim=c(common.min.cap,common.max.cap),ribbon = TRUE,main="")
plot(pp_cap,cols = "white", cex = .2, pch = 16,add=TRUE)

dev.off()






