#
# Back to brass tacks.
# Author: Kyle Taylor
#

# fake some data
X <- matrix(c(rep(1,100),rnorm(100,1:100,7)),ncol=2) # design matrix
y <- rnorm(100,1:100,12)                             # response data

solve.regression <- function(b,X){ 
  sum <- rep(0,nrow(X))
  for(i in 2:ncol(X)){
    sum <- sum + (b[i]*X[,i])
  }
  return(sum+b[1])
}

solve.betas <- function(X,y){
  as.vector(base::solve(t(X)%*%X)%*%(t(X)%*%y))
}

solve.residuals <- function(b,X,y){
  y-solve.regression(b,X)
}

np.bootstrap <- function(X,y,n=100,reportSE=T){
  betas <- NA
  for(i in 1:n){
    x <- X
    s <- sample(1:nrow(X),replace=T)
    y.bs <- y[s]
    
    for(j in 1:ncol(X)) {
      x[,j] <- X[s,j]
    }
    if(is.na(betas[1])){
      betas <- solve.betas(x,y.bs)
    } else {
      betas <- rbind(betas, solve.betas(x,y.bs))
    }
  }
  if(reportSE){
    b <- vector()
      for(i in 1:ncol(betas)) b <- append(b,sd(betas[,i])/sqrt(nrow(X)-1))
        return(b)
  } 
  return(betas)
}

solve.ser <- function(X,y){
  M <- diag(nrow(X))-(X%*%base::solve((t(X)%*%X))%*%t(X))
  L <- diag(nrow(X))-((rep(1,nrow(X))%*%t(rep(1,nrow(X))))/nrow(X))
  t(y)%*%(M%*%y) / (length(y)-(ncol(X)-1))
}

solve.rss <- function(X,y){
  nrow(X)-(ncol(X-1))/(nrow(X)) * solve.ser(X,y)
}

r.squared <- function(X,y){
  M <- diag(nrow(X))-(X%*%base::solve((t(X)%*%X))%*%t(X))
  L <- diag(nrow(X))-((rep(1,nrow(X))%*%t(rep(1,nrow(X))))/nrow(X))
  as.vector(1-((t(y)%*%(M%*%y))/(t(y)%*%(L%*%y))))
}

b<-solve.betas(X,y)
bs <- np.bootstrap(X,y,n=5000)
plot(y~X[,2])
lines(solve.regression(b,matrix(c(rep(1,119),0:118),ncol=2)),col="red")
