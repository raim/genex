
## ADJUST THESE PATHS to your git location

## relative path from <gitrepo>/scripts
source("../R/genex.R")
path <- "../data"

## input data
vanilla.time <- file.path(path,'Time series',
                          '190610 time course vanillate.csv')
vanilla.dose <- file.path(path,'Dose response',
                          '190513 vanillate dose response EVC corrected_header.csv')

### MODEL EVALUATION PARAMETER
## select between methods by Forrey (J Comp Phys 1997)
## or Stephane Laurent (stackexchange)
## NOTE: laurent is more stable for low delta!!

method <- "laurent" # "forrey"  

### GENERAL EXPERIMENT PARAMETERS 

## growth and protein degradation rates
mu <- 0.033  # h-1, 1 division per day - mu = log(2)/24
dy <- 0.064  # h-1, protein degradation rate, ASV-tagged mVenus 
beta <- mu + dy

## INDUCER HALF-LIFE - constant inducer over time
#delta <- 0 # for rhamnose?

### INDUCER-SPECIFIC PARAMETERS
y0 <- 0 # intial protein concentration
n <- 4  # Hill factor, ~ number of binding sites 
K <- 100 # half-maximal concentration 
l <- 0  # uninduced "leaky" transcription
v <- 700 # maximal induced transcription

## INDUCER HALF-LIFE - constant inducer over time

## NOTE: with current parameters, stability problems at
##  delta<.001 and delta>9 (see scan below)
#delta <- 0.05

## experimental conditions
#I0 <- c(0,10^seq(-1,4,.1)) # initial inducer concentration
#I0 <- 500
## USE NLS TO FIT DATA
timecourse= read.csv(vanilla.time,
                     header=TRUE, sep=";",dec=",")
timecourse[,2:10]
timecourse[1,2:ncol(timecourse)] <- 0

van_dose= read.csv(vanilla.dose,
                   header=TRUE,dec=',',sep=';')

van_t <- data.frame(time=rep(timecourse[,1],9), yn=unlist(timecourse[,2:10]))
van_d <- data.frame(I0=rep(van_dose[,1],18), yn=unlist(van_dose[,2:19]))
twfour <- data.frame(I0=rep(van_dose[,1],9), yn=unlist(van_dose[,2:10]),time=24)
feight <- data.frame(I0=rep(van_dose[,1],9), yn=unlist(van_dose[,11:19]),time=48)


time <-van_t[,1]
conc <-van_d[,1]
concfour <- twfour[,1]

## FIRST FIT K and v
## evaluation function for fit
f <- function(I0,K,v,delta) {
    #if ( delta < .005 ) delta <- 0
    y <- fexpr(time=24, delta=delta, I0=I0,  beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method=method)
    
    y
}

## input for fit
start <- list(K=200, v=700, delta=.05) # start value for estimation

## fitd
nlfit <- nls(yn ~ f(I0,K,v,delta), data=twfour, start=start) 

f.K <- coefficients(nlfit)["K"]
f.v <- coefficients(nlfit)["v"]
f.delta <- coefficients(nlfit)["delta"]

plot(concfour,twfour[,2], cex=1.5,pch=4, ylim=c(0,8000),
     xlab="inducer concentration", ylab="fluorescence/OD, au")   
I0 <- c(0,2*10^seq(-1,3,.1))
yi1 <- fexpr(time=24, I0=I0, delta=f.delta, beta=beta, 
             y0=y0, n=n, K=f.K, l=l, v=f.v, method=method)
lines(I0,yi1)
points(concfour,feight[,2],pch=3,col=2)
yi1 <- fexpr(time=48, I0=I0, delta=f.delta, beta=beta, 
             y0=y0, n=n, K=f.K, l=l, v=f.v, method=method)
lines(I0,yi1,col=2)


matplot(time,van_t[,2],type='p',cex=1.5,pch=4,
        xlab="time [h]", ylab=expression("RFU/OD"[750]),bty = "l")
yi1 <- fexpr(time=0:116, I0=500, delta=f.delta/2, beta=beta, 
             y0=y0, n=n, K=f.K, l=l, v=f.v, method=method)
lines(0:116, yi1)



## NLS FIT OF DATA, FOR DELTA AND ONLY

## evaluation function for fit
f2 <- function(time,delta,K) {
  #if ( delta < .005 ) delta <- 0
  y <- fexpr(time=time, delta=delta, I0=500,  beta=beta, 
             y0=y0, n=n, K=K, l=l, v=f.v, method=method)
  
  y
}

#plot(dati)

plot(concfour,twfour[,2], cex=1.5,pch=4, ylim=c(0,8000),
     xlab="inducer concentration", ylab="fluorescence/OD, au")   
plot(concfour,feight[,2],type='p',cex=1.5,pch=4,
     xlab="inducer concentration", ylab="fluorescence/OD, au") 

matplot(time,van_t[,2],type='p',cex=1.5,pch=4,xlab="time [h]", ylab=expression("RFU/OD"[750]),bty = "l")


## input for fit
start <- list(delta=f.delta/5, K=f.K) # start value for estimation
## fitd
nlfit <- nls(yn ~ f2(time, delta, K), data=van_t, start=start) 

f2.delta <- coefficients(nlfit)["delta"]
f2.K <- coefficients(nlfit)["K"]
##f2.l <- coefficients(nlfit)["l"]
##f2.v <- coefficients(nlfit)["v"]

## plot results
par(mfcol=c(1,1))

matplot(time,van_t[,2],type='p',cex=1.5,pch=4,xlab="time [h]",
        ylab=expression("RFU/OD"[750]),bty = "l")
lines(0:116, f2(0:116,f.delta,f.K),col=2)
lines(0:116, f2(0:116,f2.delta,f2.K))

##plot(time,yn,type="p",cex=1.5,pch=4,ylim=c(0,max(yn,yt)))
#lines(time,yt, col="darkgray", lwd=5, lty=2)
#lines(time, predict(nlfit,newdata=list(time=time)), col=2, lwd=2)
#par(xpd=TRUE)
legend('topright',inset=c(-0.1,0),
       legend=c("vanillate [uM]: 500","fitted time course"),
       pch=c(4,NA,NA,NA), pt.cex=c(1.5,NA,NA,NA),
       lty=c(NA, 1, NA, NA), lwd=c(NA,1,NA,NA),
       col=c("black",'black', NA, NA),bty = "n")
legend("right",#bg="#FFFFFF99",box.col=NA,
       legend=bquote("residual standard error"~
                       sigma*":"~.(round(summary(nlfit)$sigma,1))))


I0 <- c(0,2*10^seq(-1,3,.1))
plot(concfour,twfour[,2], cex=1.5,pch=4,
     ylim=c(0,8000),xlim=c(0.1,2000), log='',
     xlab="inducer concentration", ylab="fluorescence/OD, au")     
yi1 <- fexpr(time=24, I0=I0, delta=f2.delta, beta=beta, 
             y0=y0, n=n, K=f2.K, l=l, v=f.v, method=method)
lines(I0,yi1)
## original fit of K/v
yi1 <- fexpr(time=24, I0=I0, delta=f.delta, beta=beta, 
             y0=y0, n=n, K=f.K, l=l, v=f.v, method=method)
lines(I0,yi1,col=2)

I0 <- c(0,2*10^seq(-1,3,.1))
plot(concfour,feight[,2],type='p',cex=1.5,pch=4,log='',xlim=c(0.1,2000),
     xlab="inducer concentration", ylab="fluorescence/OD, au")     
yi2 <- fexpr(time=48, I0=I0, delta=f2.delta, beta=beta, 
             y0=y0, n=n, K=f2.K, l=l, v=f.v, method=method)
#ylim <- c(0, max(yi2,na.rm=TRUE))
lines(I0,yi2)
## original fit
yi1 <- fexpr(time=48, I0=I0, delta=f.delta, beta=beta, 
             y0=y0, n=n, K=f.K, l=l, v=f.v, method=method)
lines(I0,yi1,col=2)
 

plot(concfour,twfour[,2], cex=1.5,pch=4,
     ylim=c(0,8000),xlim=c(50,3000), log='x',
     xlab="inducer concentration", ylab="fluorescence/OD, au")     
points(concfour,feight[,2],pch=3,col=2)
for ( tm in seq(12,120,12) ) {
    yi2 <- fexpr(time=tm, I0=I0, delta=f2.delta, beta=beta, 
                 y0=y0, n=n, K=f2.K, l=l, v=f.v, method=method)
    lines(I0,yi2,col=2)
}

