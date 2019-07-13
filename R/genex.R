Gauss2F1 <- function(a,b,c,x, method=c("forrey","laurent") ){
  
  if ( sum(c(length(a)>1, length(b)>1, length(c)>1, length(x)>1))>1 )
    stop("Gauss2F1: only one parameter can be a vector")
  ## NOTE: for |z|>1 we use the transformation recommended in
  ## R.C. Forrey: Computing the Hypergeometric Function
  ## JOURNAL OF COMPUTATIONAL PHYSICS137,79–100 (1997)
  ## Stephane Laurent's version from stackoverflow seems both more
  ## efficient and also much more stable for some small delta
  y <- rep(list(NA),length(x))
  for ( i in 1:length(x) )
    if( abs(x[i]) < 1 ) {
      y[[i]] <- gsl::hyperg_2F1(a,b,c,x[i])
    } else {
        if ( method[1]=="forrey" ) {
            w <- 1/(1-x[i])
            A <- w^a*gamma(c)*gamma(b-a)/(gamma(b)*gamma(c-a))*gsl::hyperg_2F1(a,c-b,a-b+1,w)
            B <- w^b*gamma(c)*gamma(a-b)/(gamma(a)*gamma(c-b))*gsl::hyperg_2F1(b,c-a,b-a+1,w)
            ##if ( is.na(A) ) A <- 0
            ##if ( is.na(B) ) B <- 0
            y[[i]] <- A+B
        } else if (method[1]=="laurent" ) # Stephane Laurent's version
            y[[i]] <- gsl::hyperg_2F1(c-a,b,c,1-1/(1-x[i]))/(1-x[i])^b 
    }
  unlist(y)
}
## induced protein expression in growing cultures, with half-life of inducer
## analytic solution, protein level y(t)
fexpr <- function(delta, time=seq(1,10,.1), I0=10, y0=50,
                  n=2, K=100, l=1, v=100, beta=.1,
                  method=c("forrey","laurent")) {
  
  y <- rep(list(NA),length(I0))
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
    y[[i]] <- y0*exp(-beta*time) + (l + v*A)/beta * (1-exp(-beta*time)) 
  }
  y <- unlist(y)
  ## REMOVE y<0 - not physical, numerical artefact
  if ( any(!is.na(y)) )
    if ( any(y<0, na.rm=TRUE) ) {
      warning(sum(y<0,na.rm=TRUE), " y<0 replaced by NA")
      y[which(y<0)] <- NA
    }
  y
}

