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


par(mfcol=c(2,1), mai=c(.5,.5,.1,.1), mgp=c(1.3,.4,0),tcl=-.25)
i0 <- c(10000, 100)
ts <- c(24,48)
yt1 <- fexpr(time=time, I0=i0[1], delta=delta, beta=beta, 
             y0=y0, n=n, K=K, l=l, v=v)
yt2 <- fexpr(time=time, I0=i0[2], delta=delta, beta=beta, 
             y0=y0, n=n, K=K, l=l, v=v)
yi1 <- fexpr(time=ts[2], I0=I0, delta=delta, beta=beta, 
             y0=y0, n=n, K=K, l=l, v=v)
yi2 <- fexpr(time=ts[1], I0=I0, delta=delta, beta=beta, 
             y0=y0, n=n, K=K, l=l, v=v)

ys <- c(yt1,yt2,yi1,yi2)
ys <- ys[is.finite(ys)]
ylim <- c(0, max(ys,na.rm=TRUE))

plot(time, yt1, type="b", pch=19, cex=.5, ylim=ylim,
     ylab="fluorescence/OD, a.u.", xlab="time, h")
lines(time, yt2, col=2, type="b", pch=19, cex=.5)
legend("topright",legend=c(as.expression(bquote(I[0] ~":"~.(i0[1]))),
                           bquote(I[0]~":"~.(i0[2]))),
       title=paste0("half-life: ",tI, ",h"), col=1:2,lty=1,pch=19,pt.cex=.5)
abline(v=ts, col=3:4, lty=1)
plot(I0, yi1, type="b", col=4, lty=1, pch=19, cex=.5,
     log="x",ylim=ylim,xlim=c(1,max(I0)),
     xlab="inducer concentration", ylab="fluorescence/OD, au")
lines(I0, yi2, type="b", col=3, lty=1, pch=19, cex=.5)
abline(v=i0,col=1:2)
legend("bottomright",paste("at",ts,"h"),
       title=paste0("half-life: ",tI, ",h"),
       col=3:4,lty=1,pch=19,pt.cex=.5)

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

I0 <- seq(0,1000,10) #c(0,10^seq(-2,4,.1))
plot(1, col=NA, ylim=ylim, xlim=c(0.1,max(I0)),
     xlab="inducer concentration", ylab="fluorescence/OD, au")     
for ( tm in seq(0,72,5)) {
  yt1 <- fexpr(time=tm, I0=I0, delta=log(2)/24, beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method=method)
  lines(I0,yt1,col=1+as.numeric(tm>24))
}

## rhamnose increase to high levels over long time
beta.low <- 0.01
beta.high <- 0.1
## -> requires a low protein degradation rate
## NOTE: small delta here works bettery with laurent
yt <- fexpr(time=time, I0=100, delta=0.001, beta=beta.low, 
            y0=0, n=1, K=100, l=10, v=600, method=method)
## rhamnose step-down experiment -> to fit degradation rate
## -> should yield a high degradation rate?
## use last yt from build-up as y0
yd <- fexpr(time=time, I0=0, delta=0, beta=beta.high, 
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

## use nlm to fit

## generate data with noise
real.delta <- 0.1 # delta of simulated data
start.delta <- 0 # start value for fit
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
f <- function(time, delta) 
  fexpr(time=time, delta=delta, I0=1000, beta=beta, 
        y0=y0, n=n, K=K, l=l, v=v, method=method)
dat <- data.frame(time=ltime, yn=yn)
start <- list(delta=0.5) # start value for estimation
nlfit <- nls(yn ~ f(time, delta),data=dat,start=start) 
fitted.delta <- coefficients(nlfit)

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

