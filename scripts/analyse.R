source("../R/genex.R")

### MODEL EVALUATION PARAMETER
## select between methods by Forrey (J Comp Phys 1997)
## or Stephane Laurent (stackexchange)
## NOTE: laurent is more stable for low delta!!
method <- "laurent" # "forrey" # 

### GENERAL EXPERIMENT PARAMETERS 

## growth and protein degradation rates
mu <- 0.029 # h-1, 1 division per day - mu = log(2)/24
dy <- 0.2  # h-1, protein degradation rate, ASV-tagged mVenus 
beta <- mu + dy

### INDUCER-SPECIFIC PARAMETERS
y0 <- 50 # intial protein concentration
n <- 2   # Hill factor, ~ number of binding sites 
K <- 100 # half-maximal concentration 
l <- 1   # uninduced "leaky" transcription
v <- 100 # maximal induced transcription

## INDUCER HALF-LIFE - constant inducer over time
tI <- 17.5 # h, half-life aTc
delta <- log(2)/tI # h-1, degradation rate aTc
#delta <- 0 # for rhamnose?

## NOTE: with current parameters, stability problems at
##  delta<.001 and delta>9 (see scan below)
#delta <- 0.001

## experimental conditions
time <- 0:150
I0 <- c(0,10^seq(-1,4,.1)) # initial inducer concentration

## activation term as a function of time
par(mai=c(.5,.5,.1,.1), mgp=c(1.3,.4,0),tcl=-.25)
plot(I0, Gauss2F1(1, beta/(n*delta), beta/(n*delta) +1, -exp(n*delta*0)*K^n/I0^n),
     type="l", xlab="initial inducer concentration", ylab="activation term")
for ( tm in seq(0,150,10) )
    lines(I0, Gauss2F1(1, beta/(n*delta), beta/(n*delta) +1, -exp(n*delta*tm)*K^n/I0^n))
arrows(x0=2000,x1=5000,y0=0.8,y1=0.4)
text(3500,0.6, "increasing time",pos=4)

## increasing delta vs. time
ylim <- c(0, 425)
time <- 0:150
plot(1, col=NA,  ylim=ylim,xlim=range(time),
     xlab="time, h", ylab="fluorescence/OD, au")
for ( delta in seq(0,8,.05)/10 ) {
  ## NOTE: instabilities also at delta>8 at later time-points!
  ## sudden decrease to 0
  yt1 <- fexpr(time=time, I0=1000, delta=delta, beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method=method)
  lines(time,yt1)
}
arrows(x0=120,x1=80,y0=350,y1=250)
text(x=100,y=300,expression("increasing"~delta),pos=4)

## increasing inducer at different sampling times
I0 <- seq(0,1000,10) #c(0,10^seq(-2,4,.1))
plot(1, col=NA, ylim=ylim, xlim=c(0.1,max(I0)),
     xlab="inducer concentration", ylab="fluorescence/OD, au")     
for ( tm in seq(0,72,5)) {
  yt1 <- fexpr(time=tm, I0=I0, delta=log(2)/24, beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method=method)
  lines(I0,yt1,col=1+as.numeric(tm>24))
}


## rhamnose increase to high levels over long time
time <- 0:70
delta.rha <- 0.005
deltaP.low <- 0.01
deltaP.high <- 0.05
beta.low <- deltaP.low + mu
beta.high <- deltaP.high + mu
## -> requires a low protein degradation rate
## NOTE: small delta here works bettery with laurent
yt <- fexpr(time=time, I0=100, delta=delta.rha, beta=beta.low, 
            y0=0, n=1, K=100, l=10, v=600, method=method)
## rhamnose step-down experiment -> to fit degradation rate
## -> should yield a high degradation rate?
## use last yt from build-up as y0
yd <- fexpr(time=time, I0=0, delta=delta.rha, beta=beta.high, 
            y0=yt[length(yt)], n=n, K=K, l=l, v=v, method=method)

## plot
par(mfcol=c(1,1))
plot(time, yt, xlim=c(0,2*time[length(time)]),type="l")
legend("topleft",expression("slow increase implies low"~delta[P]~"..."),
       box.col=NA, bg="#FFFFFF99")
lines(time+time[length(time)], yd,col=2)
legend("right",expression("... but step-down implies high"~delta[P]),
       bty="n", text.col=2)
legend("topright","rhamnose time-series")
text(time[length(time)],yt[length(yt)]*.75, 
     bquote(delta[P]*"="*.(beta.low-mu)),pos=2)
text(time[length(time)],yt[length(yt)]*.75, 
     bquote(delta[P]*"="*.(beta.high-mu)),pos=4,col=2)

## USE NLS TO FIT

## generate data with noise
real.delta <- 0.1 # delta of simulated data
start.delta <- 0.05 # start value for fit
## real data at high res
time <- seq(0,70,1)
yt <- fexpr(time=time, I0=1000, delta=real.delta, beta=beta, 
            y0=y0, n=n, K=K, l=l, v=v, method=method)

## "experimental data" at low res, and multiple samples per time-point
ltime <- rep(seq(0,max(time),length.out=10),3)
yn <- fexpr(time=ltime, I0=1000, delta=real.delta, beta=beta, 
            y0=y0, n=n, K=K, l=l, v=v, method=method)
yn <- c(yn + rnorm(length(ltime),mean=0, sd=30))

## NLS FIT OF DATA, FOR DELTA ONLY

## evaluation function for fit
f <- function(time, delta, l, v, K) {  
    #if ( delta < .005 ) delta <- 0
    y <- fexpr(time=time, delta=delta, I0=1000, beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method=method)
    
    y
}

## input for fit
dat <- data.frame(time=ltime, yn=yn)
start <- list(delta=start.delta, l=l, v=v, K=K) #, beta=beta) # start value for estimation

## fit
nlfit <- nls(yn ~ f(time, delta, l, v, K), data=dat, start=start) 

fitted.delta <- coefficients(nlfit)[1]

## plot results
par(mfcol=c(1,1))
plot(ltime,yn,type="p",cex=1.5,pch=4,ylim=c(0,max(yn,yt)))
lines(time,yt, col="darkgray", lwd=5, lty=2)
lines(time, predict(nlfit,newdata=list(time=time)), col=2, lwd=2)
legend("topright",c("data with noise",
                 paste0(c("real","start","fitted"),": ",
                     c(real.delta, start.delta,round(fitted.delta,3)))),
       title=expression("fit"~delta[P]),
       pch=c(4,NA,NA,NA), pt.cex=c(1.5,NA,NA,NA),
       lty=c(NA, 2, NA, 1), lwd=c(NA,5,NA,2),
       col=c("black","darkgray", NA, "red"))
legend("right",#bg="#FFFFFF99",box.col=NA,
       legend=bquote("residual standard error"~
                       sigma*":"~.(round(summary(nlfit)$sigma,1))))

