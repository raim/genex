source("../R/genex.R")

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

time <- 0:100
plot(1, col=NA,  ylim=ylim,xlim=range(time),
     xlab="time, h", ylab="fluorescence/OD, au")
for ( delta in seq(0,8,.05)/10 ) {
  ## NOTE: instabilities also at delta>8 at later time-points!
  ## sudden decrease to 0
  yt1 <- fexpr(time=time, I0=1000, delta=delta, beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method="laurent")
  lines(time,yt1)
}

I0 <- seq(0,1000,10) #c(0,10^seq(-2,4,.1))
plot(1, col=NA, ylim=ylim, xlim=c(0.1,max(I0)),
     xlab="inducer concentration", ylab="fluorescence/OD, au")     
for ( tm in seq(0,72,5)) {
  yt1 <- fexpr(time=tm, I0=I0, delta=log(2)/24, beta=beta, 
               y0=y0, n=n, K=K, l=l, v=v, method="laurent")
  lines(I0,yt1,col=1+as.numeric(tm>24))
}

## rhamnose step-down experiment -> to fit degradation rate
yt <- fexpr(time=time, I0=0, delta=0, beta=.2, 
            y0=6000, n=n, K=K, l=l, v=v, method="laurent")
plot(time, yt)
plot(time, log(yt))

## rhamnose increase to high levels over long time
## -> requires a high degradation rate
yt <- fexpr(time=time, I0=100, delta=0, beta=0.01, 
            y0=0, n=1, K=10, l=100, v=500, method="laurent")
plot(time, yt)


## use nlm to fit

## generate data with noise
yt <- fexpr(time=time, I0=1000, delta=0.1, beta=beta, 
            y0=y0, n=n, K=K, l=l, v=v, method="laurent")
yn <- yt + rnorm(length(time),mean=0, sd=50)
plot(time,yt,type="l")
points(time,yn)
