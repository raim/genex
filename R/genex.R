Gauss2F1 <- function(a,b,c,x, method=c("forrey","laurent") ){

  ## TODO: avoid this and allow vectors of equal length
  ## the loop below, checking for x values, must be re-formulated
  if ( sum(as.numeric(c(length(a)>1, length(b)>1, 
                        length(c)>1, length(x)>1)))>1 )
    stop("Gauss2F1: only one parameter can be a vector")
  
  ## NOTE: for |z|>1 we use the transformation recommended in
  ## R.C. Forrey: Computing the Hypergeometric Function
  ## JOURNAL OF COMPUTATIONAL PHYSICS 137,79–100 (1997)
  ## Stephane Laurent's version from stackexchange seems both more
  ## efficient and also much more stable for some small delta
  ## https://stats.stackexchange.com/questions/33451/computation-of-hypergeometric-function-in-r
  
  ## loop through x to check values
  ## TODO: use apply?
  ## TODO: use different x[i] selection for laurent
  y <- rep(list(NA),length(x))
  for ( i in 1:length(x) )
    if ( method[1]=="forrey" ) {
      if( abs(x[i]) < 1 ) {
        y[[i]] <- gsl::hyperg_2F1(a,b,c,x[i])
      } else {
        w <- 1/(1-x[i])
        A <- w^a*gamma(c)*gamma(b-a)/(gamma(b)*gamma(c-a))
        A <- A*gsl::hyperg_2F1(a,c-b,a-b+1,w)
        B <- w^b*gamma(c)*gamma(a-b)/(gamma(a)*gamma(c-b))
        B <- B*gsl::hyperg_2F1(b,c-a,b-a+1,w)
        y[[i]] <- A+B
      } 
    } else if (method[1]=="laurent" ) { # Stephane Laurent's version
      if(x[i]>=0 & x[i]<1){
        y[[i]] <- gsl::hyperg_2F1(a,b,c,x[i])
      }else{
        y[[i]] <- gsl::hyperg_2F1(c-a,b,c,1-1/(1-x[i]))/(1-x[i])^b 
      }
    }
  unlist(y)
}
## induced protein expression in growing cultures, with half-life of inducer
## analytic solution, protein level y(t)
fexpr <- function(delta, time=seq(1,10,.1), I0=10, y0=50,
                  n=2, K=100, l=1, v=100, beta=.1,
                  method=c("forrey","laurent")) {
  
  ## TODO: fix this - they could be vectors but 2F1 interface
  ## does currently not allow it
  if ( any(c(length(n)>1, length(beta)>1,length(delta)>1)))
    stop("fexpr: currently n, beta and delta can not be vectors")
  if ( sum(as.numeric(c(length(time)>1, length(I0)>1, 
                        length(y0)>1, length(K)>1, 
                        length(l)>1, length(v)>1)))>1 ) 
    stop("fexpr: only one parameter can be a vector")
  
  y <- rep(list(NA),length(I0))
  ## TODO: use apply instead of loop
  for ( i in 1:length(I0) ) {
    i0 <- I0[i]  
    if ( delta == 0 | i0 == 0) {
      A <- i0^n/(i0^n + K^n) # induced transcription 
    } else {
      A <- Gauss2F1(1, 
                    beta/(n*delta),
                    beta/(n*delta)+1, 
                    -exp(n*delta*time)*K^n/i0^n, method=method[1])
    }
    if ( any(A<0, na.rm=TRUE) ) warning("non-physical value for A")
    y[[i]] <- y0*exp(-beta*time) + (l + v*A)/beta * (1-exp(-beta*time)) 
  }
  y <- unlist(y)
  ## REMOVE y<0 - not physical, numerical artefact
  if ( any(!is.na(y)) )
    if ( any(y<0, na.rm=TRUE) ) {
      warning(sum(y<0,na.rm=TRUE), " y<0")
      #y[which(y<0)] <- NA
    }
  y
}

