require(stats)
require(gtools)

source("Auxiliary_functions.R")

################################################################################

#Simulations

par(mar = c(5.1, 4.1, 2.1, 2.1))

#Gaussian example 1

set.seed(211926)

#We generate a true function s from m points
{
  m <- 4
  
  data <- Gen_Data(m, c(0,1))
  
  true_s <- approxfun(data[,1], data[,2])
  
}

#We generate simulated data from s
{
  n <- 100
  
  X <- runif(n, min = 0, max = 1)
  epsilon <- rnorm(n, 0, 0.05)
  
  Y <- true_s(X)
  
  amostra <- as.data.frame(cbind(X, Y + epsilon))
  amostra <- amostra[order(amostra$X),]
  colnames(amostra) <- c("X","Y")

  plot(amostra$X, amostra$Y, pch = 16, col = 'grey', xlab = "X", ylab = "Y")
  curve(true_s, col = 'grey', lwd = 2, add = TRUE)
}

#Plot with CVSS and w
{
  
  w <- w_search(amostra,100)
  
  plot(seq(from = 4, to = 2*(n-1), by = 2)/(n-1), w[[2]], pch = 16, cex = 0.8,
       xlab = "w", ylab = "CVSS(w)")
}

#Smoothers
{
  
  #optimal w
  Dados2 <- smoother(amostra, w[[1]], 'runlin', TRUE)
  
  #w lesser than the optimal value
  Dados3 <- smoother(amostra, ((w[[1]]*(n-1) + 4)/2)/(n-1), 'runlin', TRUE)
  
  #w greater than the optimal value
  Dados4 <- smoother(amostra, ((w[[1]]*(n-1) + 2*(n-1))/2)/(n-1), 'runlin', TRUE)
  
  plot(amostra$X, amostra$Y, pch = 16, col = 'grey', xlab = "X", ylab = "Y")
  curve(true_s, col = 'grey', lwd = 2, add = TRUE)
  
  lines(Dados3$X, Dados3$Ysmoo, col = 'blue', lwd = 2)
  lines(Dados4$X, Dados4$Ysmoo, col = 'red', lwd = 2)
  lines(Dados2$X, Dados2$Ysmoo, col = 'green', lwd = 2)
  
  legend("topright", legend = c("w_opt","w < w_opt", "w > w_opt"),
         col = c("green", "blue", "red"), pch = c(16,16,16),
         pt.cex = 1, 
         cex = 0.8, 
         text.col = "black", 
         horiz = F)
}




#Gaussian example 2

set.seed(211927)

#We generate a true function s from m points
{
  m <- 4
  
  data <- Gen_Data(m, c(0,1))
  
  true_s <- approxfun(data[,1], data[,2])
}

#We generate simulated data from s
{
  n <- 200
  
  X <- runif(n, min = 0, max = 1)
  epsilon <- rnorm(n, 0, 0.05)
  
  Y <- true_s(X)
  
  amostra <- as.data.frame(cbind(X, Y + epsilon))
  amostra <- amostra[order(amostra$X),]
  colnames(amostra) <- c("X","Y")
  
  plot(amostra$X, amostra$Y, pch = 16, col = 'grey', xlab = "X", ylab = "Y")
  curve(true_s, col = 'grey', lwd = 2, add = TRUE)
  
  
  
  w <- w_search(amostra,20)[[1]]
  
  
  #optimal w
  Dados2 <- smoother(amostra, w, 'runlin', FALSE)
  
  #w lesser than the optimal value
  Dados3 <- smoother(amostra, ((w*(n-1) + 4)/2)/(n-1), 'runlin', FALSE)
  
  #w greater than the optimal value
  Dados4 <- smoother(amostra, ((w*(n-1) + 2*n)/2)/(n-1), 'runlin', FALSE)
  
  
  
  lines(Dados3$X, Dados3$Ysmoo, col = 'blue', lwd = 2)
  lines(Dados4$X, Dados4$Ysmoo, col = 'red', lwd = 2)
  lines(Dados2$X, Dados2$Ysmoo, col = 'green', lwd = 2)
  
  legend("topright", legend = c("w_opt","w < w_opt", "w > w_opt"),
         col = c("green", "blue", "red"), pch = c(16,16,16),
         pt.cex = 1, 
         cex = 0.8, 
         text.col = "black", 
         horiz = F)
}

################################################################################

#Logistic case


#Logistic example 1

set.seed(211926)
{
#m is the amount of points used to generate the true function s
m <- 4

#N is the size parameter of the binomial distribution of Y given X=x (we assume N constant in x)
N <- 3

#we generate the set of points which will determine s
#This set indicates the true probabilities associated with some x values.
Dados <- Gen_Data(m, c(0,1))

#Now we calculate the true eta values and the corresponding piecewise linear function
Dados$V3 <- logit(Dados$V2)
true_s <- approxfun(Dados$V1, Dados$V3)

#Now we calculate the function associated with the true probabilities (while not a piecewise linear
#function, it is useful for visualization purposes)
true_prob <- function(x){inv.logit(true_s(x))}

#plot of x vs predictor
curve(true_s, from=0, to=1, xlab="x", ylab="eta(x)", col = 'grey')

#plot of x vs probabilities
curve(true_prob, from=0, to=1, xlab="X", ylab="p(X)", ylim = c(0,1), col = 'grey')
}
{
  n <- 100
  
  X <- runif(n, min = 0, max = 1)
  
  pX <- true_prob(X)
  
  Y <- rbinom(n,size = N,prob = pX)/N
  
  amostra <- as.data.frame(cbind(X,Y))
  amostra <- amostra[order(amostra$X),]
  
  #initial guesses for eta
  init_guess <- as.data.frame(rep(0,n))
  
  eta <- Fitting(init_guess,amostra, 0.01)
  
  est_prob <- approxfun(amostra$X, inv.logit(eta))
  
  curve(true_prob, from=0, to=1, xlab="X", ylab="p(X)", ylim = c(0,1), col = 'grey', lwd = 2)
  points(amostra$X, amostra$Y, pch = 16, col = 'grey') 
  curve(est_prob, from=0, to=1, ylim = c(0,1), col = 'green', add = TRUE, lwd = 2)
}


#Logistic example 2

set.seed(211926)
{
  #m is the amount of points used to generate the true function s
  m <- 8
  
  #N is the size parameter of the binomial distribution of Y given X=x
  N <- function(x){return(floor(50*x+1))}
  
  #we generate the set of points which will determine s
  #This set indicates the true probabilities associated with some x values.
  Dados <- Gen_Data(m, c(0,1))
  
  #Now we calculate the true eta values and the corresponding piecewise linear function
  Dados$V3 <- logit(Dados$V2)
  true_s <- approxfun(Dados$V1, Dados$V3)
  
  #Now we calculate the function associated with the true probabilities
  true_prob <- function(x){inv.logit(true_s(x))}
  
  #plot of x vs predictor
  curve(true_s, from=0, to=1, xlab="X", ylab="eta(X)", col = 'grey')
  
  #plot of x vs probabilities
  curve(true_prob, from=0, to=1, xlab="X", ylab="p(X)", ylim = c(0,1), col = 'grey')
}

{
  n <- 200
  
  X <- runif(n, min = 0, max = 1)
  
  pX <- true_prob(X)
  
  rbinom_sizes <- function(args){return(rbinom(1,args[1],prob = args[2])/args[1])}
  
  Y <- apply(cbind(N(X),pX),FUN= rbinom_sizes, MARGIN = 1)
  
  amostra <- as.data.frame(cbind(X,Y))
  amostra <- amostra[order(amostra$X),]
  
  #initial guesses for eta
  init_guess <- as.data.frame(rep(0,n))
  
  eta <- Fitting(init_guess,amostra, 0.01)
  
  est_prob <- approxfun(amostra$X, inv.logit(eta))
  
  curve(true_prob, from=0, to=1, xlab="X", ylab="p(X)", ylim = c(0,1), col = 'grey', lwd = 2)
  points(amostra$X, amostra$Y, pch = 16, col = 'grey') 
  curve(est_prob, from=0, to=1, ylim = c(0,1), col = 'green', add = TRUE, lwd = 2)
}




