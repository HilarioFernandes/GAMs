#This function returns the limits of neighbourhoods with m and w given
neighbourhoods <- function(m, w){
  
  indexes <- c(1:m)
  
  indexes <- cbind(indexes- floor((floor(w*m)-1)/2),rep(1,m), indexes+ floor((floor(w*m)-1)/2), rep(m,m) )
  
  return(cbind(pmax(indexes[,1], indexes[,2]), pmin(indexes[,3], indexes[,4])))
}

#This function calculates a smoother
#type can assume the values 'movavg' (for moving averages) or 'runlin' (for running lines)
smoother <- function(data, w, type, print_or_not){
  
  data2 <- as.data.frame(data)
  colnames(data2) <- c("X","Y")
  
  m <- nrow(data)
  
  indexes <- neighbourhoods(m, w)
  
  if (print_or_not == TRUE){
    print(indexes)
  }
  
  
  
  
  if (type == 'movavg'){
    
    data2$Ysmoo <- 0
    
    for(i in 1:m){
      data2$Ysmoo[i] <- mean(data[indexes[i,1]:indexes[i,2],2]) 
    }
    
  }
  
  if (type == 'runlin'){
    
    data2$Xbar <- 0
    data2$Ybar <- 0
    data2$beta_1 <- 0
    data2$beta_0 <- 0
    
    for(i in 1:m){
      data2$Xbar[i] <- mean(data[indexes[i,1]:indexes[i,2],1])
      data2$Ybar[i] <- mean(data[indexes[i,1]:indexes[i,2],2]) 
    }
    
    for(i in 1:m){
      data2$beta_1[i] <- sum( ((data2$Y-data2$Ybar)*(data2$X-data2$Xbar))[indexes[i,1]:indexes[i,2]] )/sum( ((data2$X-data2$Xbar)^2)[indexes[i,1]:indexes[i,2]] )
    }
    
    #print(data2[,2])
    #print(data2[,1])
    #print(data2$Xbar)
    
    data2$beta_0 <- data2$Ybar-data2$beta_1*data2$Xbar
    
    data2$Ysmoo <- data2$beta_0 + data2$beta_1*data2$X
    
  }
  if (print_or_not == TRUE){
    print(data2)
  }
  
  return(data2)
  
}


#This function calculates the cross-validation error of a smoother with running lines
#given the data and w
#The data.frame data must possess columns with names "X" and "Y"
CVSS <- function(data,w){
  
  m <- nrow(data)
  
  sum <- 0
  
  for(i in 2:(m-1)){
    temp_data <- smoother(data[-i,],w,'runlin', FALSE)
    
    temp_function <- approxfun(temp_data$X, temp_data$Ysmoo)
    
    sum <- sum + (data$Y[i] - temp_function(data$X[i]))^2
    
    #print(c(data$Y[i], data$X[i], temp_function(data$X[i])))
    
  }
  
  temp_data <- smoother(data[-1,],w,'runlin', FALSE)
  
  sum <- sum + (data$Y[1] - temp_data$X[2])^2
  
  temp_data <- smoother(data[-m,],w,'runlin', FALSE)
  
  sum <- sum + (data$Y[m] - temp_data$X[m-1])^2
  
  return(sum/m)
}

#This function searches for the value of w that minimizes the CVSS.
#Although the cross-validation error can be calculated with w in a real interval,
#it is sufficient to check the CVSS for finitely many values of w.
#These are of the form j/(m-1) for j even between 4 and 2(m-1), where m is the number of 
#rows in data (number of observations).

#tolerance is an useful parameter. Given that tolerance iterations have been executed
#after a local minimum has been found, the algorithm stops. If truncation isn't desirable,
#it is enough to set tolerance <- nrow(data). In some cases tolerance <- 20 is useful.
w_search <- function(data, tolerance){
  
  m <- nrow(data)
  
  temp_CVSS <- Inf
  CVSS_arm <- temp_CVSS
  
  sequence <- seq(from = 4, to = 2*(m-1), by = 2)/(m-1)
  
  count <- 0
  
  for (temp_w in sequence){
    
    count <- count + 1
    
    CVSS_arm <- c(CVSS_arm, CVSS(data, temp_w))
    
    if(CVSS(data, temp_w) < temp_CVSS){
      temp_CVSS <- CVSS(data, temp_w)
      w <- temp_w
      print("new w found:")
      print(paste0("w = ",w,", ","CVSS(w) = ", temp_CVSS))
      
      count <- 0
    }
    
    if (count >= tolerance){
      return(list(w, CVSS_arm[-1]))
    }
  }
  
  return(list(w, CVSS_arm[-1]))
}

#This function generates m points on the set [0,1] x [a,b] such that one of them 
#has a 0 as its first coordinate, and another one has a 1 as its first coordinate.
#The output of this function will be used to generate a true function s.
Gen_Data <- function(m, xlim){
  
  Data <- matrix(runif(2*m, min = xlim[1], max = xlim[2]), ncol = 2)
  Data[1,1] <- 0
  Data[2,1] <- 1
  
  Data <- as.data.frame(Data)
  Data <- Data[order(Data$V1),]
  
}

#This function adjusts the smoother using the iterative process described in the paper
#-init_guess is a vector with starting guesses for eta values
#-sample is a data.frame with the observations X and Y
#-epsilon is the desired precision (sup norm of the difference between two consecutive iterations)
Fitting <- function(init_guess, sample, epsilon){
  
  eta <- init_guess
  eta_hist <- eta
  delta <- Inf
  
  while (delta > epsilon){
    
    probs <- t(inv.logit(t(eta)))
    adj_Y <- eta + (sample$Y - probs)/(probs*(1-probs))
    
    
    #escolhemos o w por validação cruzada
    
    data <- as.data.frame(cbind(sample$X, adj_Y))
    colnames(data) <- c("X","Y")
    
    w <- w_search(data,20)[[1]]
    
    eta_hist <- cbind(eta_hist,smoother(data, w, 'runlin', FALSE)$Ysmoo) 
    eta <- smoother(data, w, 'runlin', FALSE)$Ysmoo
    
    delta <- max(abs(eta_hist[,ncol(eta_hist)]-eta_hist[,ncol(eta_hist)-1]))
    print(paste0("delta = ", delta))
  }
  
  return(eta)
  
}