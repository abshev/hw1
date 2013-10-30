
##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(mvtnorm)
library(coda)

########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################


"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv,
                           niter=10000,burnin=1000,
                           print.every=1000,retune=100,
                           verbose=TRUE){
  #Intitiate variables used in the algorithm:
  #-p: number of parameters
  #-draws: (ninter+burnin)xp matrix to store the chain
  #-propvar: length p vector of proposal variances
  #-accept.count: length p vector that records the number of successful jumps for
  #               each parameter
  p = length(beta.0)
  draws = matrix(nrow = burnin + niter, ncol = p)
  draws[1,] = beta.0
  propvar = 1 / diag(Sigma.0.inv) 
  count = numeric(p)
  accept.count = numeric(p)
  #The chains runs in two loops, i loops over draws, and j loops over parameters
  for(i in 2:burnin){
    draws[i,] = draws[i - 1,]
    for(j in 1:p){
      #For each parameter a metropolis-hastings draw is made
      mh = marginalMH(beta = draws[i,], ind = j, m = m, y = y, X = X,
                      beta.0 = beta.0, Sigma.0.inv = Sigma.0.inv, propvar = propvar)
      #If the draw is successful the accept count is incremented.
      if(mh$jump){
        accept.count[j] = accept.count[j] + 1
      }
      #The new draw is added to the draws matrix and the number of reps is incremented
      draws[i,j] = mh$draw
      #Every retune iterations, the proposal variance is adjusted
      if(i == retune){
        propvar[j] = retuneVar(propvar[j], accept.count[j]/retune)
        accept.count[j] = 0
      }
    }
    #Prints the iteration number every print.every interations
    if(i %% print.every == 0 & verbose){
      print(i)
      flush.console() 
    }
  }
  #The loops are repeated for the post-burnin iterations
  accept.count = c(0,0)
  for(i in 1:niter){
    draws[burnin + i,] = draws[burnin + i - 1,]
    for(j in 1:p){
      mh = marginalMH(beta = draws[burnin + i,], ind = j, m = m, y = y, X = X,
                      beta.0 = beta.0, Sigma.0.inv = Sigma.0.inv, propvar = propvar)
      if(mh$jump){
        accept.count[j] = accept.count[j] + 1
      }
      draws[i + burnin,j] = mh$draw
    }
    if(i %% print.every == 0){
      print(i)
      flush.console() 
    }
  }
  #Calculates the quantiles of b0 and b1.
  post.dists = cbind(beta0 = quantile(ecdf(draws[-(1:(burnin + 1)),1]), 
                                      seq(from = .01, to = .99, by = .01)),
                     beta1 = quantile(ecdf(draws[-(1:(burnin + 1)),2]),
                                      seq(from = .01, to = .99, by = .01)))
  return(post.dists)
}







#THis function runs the MH step
marginalMH = function(beta, ind, m, y, X, beta.0, Sigma.0.inv, propvar){
  #The draw is initiated as the previous value
  draw = beta[ind]
  #This flag let's the main function know if the jump was successful for
  #the purposes of tracking the acceptance rate
  jump = FALSE
  beta.new = beta
  #The new beta is drawn from the proposal density
  beta.new[ind] = rnorm(1, beta[ind], sqrt(propvar[ind]))
  #a is the ratio of likelihoods
  a = exp(logLik(beta.new,  m , y, X, beta.0, Sigma.0.inv, ind) - 
            logLik(beta, m, y, X, beta.0, Sigma.0.inv, ind))
  #accept new jump with probability min(a,1)
  if(runif(1) <= a){
    draw = beta.new[ind]
    jump = TRUE
  }
  return(list(draw = draw, jump = jump))
}

#This function retunes the variance if the acceptance rate is outside the
#acceptable range.
retuneVar = function(var, acceptance){
  if(.25 <= acceptance & acceptance <= .55){
    return(var)
  }
  else{
    return(var * exp(20*(acceptance - .35)))
  }
}

#Computes the value of the log-likelihood given data and parameters
logLik = function(beta, m, y, X, beta.0, Sigma.0.inv, ind){
  p = inv.logit(X%*%beta)
  #fix for machine rounding error
  p[p==1] = .999999
  sum(y * log(p / (1 - p)) + m * log(1 - p) - 
        .5 * ((beta[ind] - beta.0[ind])^2) * Sigma.0.inv[ind, ind])
}

#Inverse logit function
inv.logit = function(x){
  exp(x)/(1+exp(x))
}


#################################################
# Set up the specifications:
beta.0 <- matrix(c(0,0))
Sigma.0.inv <- diag(2)
niter <- 10000
burnin <- 20000
print.every <- 1000
retune <- 50
# etc... (more needed here)
#################################################

# Read data corresponding to appropriate sim_num:
blr = read.csv(paste("~/HW1/BayesLogit/data/blr_data_", sim_num, ".csv", sep = ""))
# Extract X and y:

# Fit the Bayesian model:
fit = bayes.logreg(blr$n, blr$y, cbind(blr$X1, blr$X2), beta.0, Sigma.0.inv, 
                   niter,burnin,print.every,retune,verbose=TRUE) 
# Extract posterior quantiles...

# Write results to a (99 x p) csv file...
write.table(fit, file = paste("~/HW1/BayesLogit/results/blr_res_", 
                              sim_num, ".csv", sep = ""), sep = ",",
            row.names = FALSE, col.names = FALSE)
# Go celebrate.
 
cat("done. :)\n")







