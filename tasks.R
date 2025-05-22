## project bayesian inference

# install.packages(c("readr","coda","runjags","MCMCvis","ggmcmc","basicMCMCplots"))
# install.packages(c("igraph", "pracma", "numDeriv", "R6"))
library(tidyverse)
library(R2OpenBUGS)
library(rjags)
library(nimble)
# nimble



## rebuild data
dat <- matrix(c(380, 419,
                632, 673,
                28, 34,
                363, 396,
                527, 576,
                36, 50,
                282, 332,
                514, 548,
                16, 34,
                446, 490,
                588, 628,
                28, 39,
                400, 441,
                27, 32), nrow = 14, byrow = TRUE)

dat <- as.data.frame(dat)
dat[,2] <- dat[,2] - dat[,1]
# dimnames(dat) <- list(c("North Carolina", "Georgia", "Wisconsin", "Florida", "Mississippi"))
geography <- rep(c("North Carolina", "Georgia", "Wisconsin", "Florida", "Mississippi"), each = 3)
geography <- geography[-length(geography)]
insurance <- rep(c("Any Medicaid", "Private Insurance Only", "Uninsured"), 5)
# remove row: Missisippi has no Any Mediaid data
insurance <- insurance[-13]
dat <- cbind(geography, insurance, dat)
names(dat) <- c("geography", "insurance", "yes", "no")
dat$Sum <- dat$yes + dat$no
dat$geography <- as.factor(dat$geography)
dat$insurance <- as.factor(dat$insurance)


## for more details: modify dataset, create collapsed datasets
dat_margins <- cbind(rbind(ungroup(dat)[,1:2], NA), addmargins(as.matrix(dat[-(1:2)])))
# dataframes collapsed by geography and insurance
dat_geography <- dat %>% group_by(geography) %>%
  summarize(yes = sum(yes),
            no = sum(no),
            sum = sum(Sum))
totN <- 4692
dat_geo_prob <- dat_geography %>% mutate(across(where(is.numeric), ~ .x / totN))
dat_insurance <- dat %>% group_by(insurance) %>%
  summarize(yes = sum(yes),
            no = sum(no),
            sum = sum(Sum))
dat_ins_prob <- dat_insurance %>% mutate(across(where(is.numeric), ~ .x / totN))



####### task 1
get_post_pars <- function(dat, accrate=0) {
  # inform user about argument behaviour of accrate
  message("Acceptance rate should be set to zero if an uninformative prior is chosen")
  
  # initiate new columns
  dat$posterior_mode <- NA
  dat$posterior_mean <- NA
  dat$posterior_median <- NA
  dat$posterior_variance <- NA
  dat$lower_ci <- NA
  dat$upper_ci <- NA
  
  # loop over table rows and calculate posterior summary statistics
  for (i in 1:nrow(dat)) {
    
    y <- dat$yes[i]
    n <- dat$Sum[i]
    
    

    if (accrate == 0) {
      # for uninformative prior: 
      y_prior <- y
    }
    else
      # for informative prior: > xx% acceptance rate
      y_prior <- (dat$Sum[i] * accrate)
      
    # calculate prior values
    alpha_prior <- y_prior + 1
    beta_prior <- n - y_prior + 1 
    
    # calculate scale parameters
    alpha_bar <- alpha_prior + y
    beta_bar <- beta_prior + n - y
    
    # calculate summary statistics
    dat$posterior_mean[i] <- alpha_bar/(alpha_bar+beta_bar)
    dat$posterior_mode[i] <- (alpha_bar-1)/(alpha_bar+beta_bar-2)
    dat$posterior_median[i] <-  qbeta(p = 0.5, shape1 = alpha_bar, shape2 = beta_bar)
    dat$posterior_variance[i] <- (alpha_bar*beta_bar)/((alpha_bar+beta_bar)^2*(alpha_bar+beta_bar+1))
    dat$lower_ci[i] <- qbeta(0.025, alpha_bar, beta_bar)
    dat$upper_ci[i] <- qbeta(0.975, alpha_bar, beta_bar)

  }
  # return extended dataframe
  return(dat)
}

tab_summary_inf_prior <- get_post_pars(dat, 0.9)
tab_summary_uninf_prior <- get_post_pars(dat, 0)


####### task 2
N <- nrow(dat)

dummy_AnyMedicaid <- ifelse(dat$insurance == "Any Medicaid", 1, 0)
dummy_uninsured <- ifelse(dat$insurance == "Uninsured", 1, 0)

model.data <- list('vaccinated' = dat$yes) # 
model.constant <- list('N' = N, 'anymed' = dummy_AnyMedicaid, 'unins' = dummy_uninsured, 'total' = dat$Sum)

model.inits <- list(alpha0=1, alpha1=1, alpha2=1)

model_set <- nimbleCode(
  {
    for (i in 1:N) {
      logit(p[i]) <- alpha0 + alpha1*anymed[i] + alpha2*unins[i]
      #y_ij[i] <- alpha0 + alpha1*anymed[i] + alpha2*unins[i]
      vaccinated[i] ~ dbinom(p[i], total[i])
      #vaccinated[i] ~ dbinom()
    }
    # prior info (only when uniformative)
    alpha0 ~ dbeta(shape1=1, shape2=1)
    alpha1 ~ dbeta(shape1=1, shape2=1)
    alpha2 ~ dbeta(shape1=1, shape2=1)
    # -> how to specify prior knowledge that y is 90% or higher?
    
    
    # prior info: dnorm instead of dbeta?
    #prob ~ dbeta(shape1=1, shape2=1)
    #size ~ dbeta(shape1=1, shape2=1)  
    
    # is it correct to model the size between 0-1?
    # prior should be uninformative
    # but total sample sizes are integers and between 1 and 500
  }
)

model <- nimbleModel(model_set, 
                     constants=model.constant, 
                     data=model.data, 
                     inits=model.inits)

model.compiled <- compileNimble(model)

out <- nimbleMCMC(model,
                  niter=10000,
                  nburnin=5000,
                  thin=1,
                  summary=TRUE)

print(out)
plot(out$samples)


dataNodes <- model$getNodeNames(dataOnly = TRUE)
parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
simNodes <- model$getDependencies(parentNodes, self = FALSE)


####### task 3
# how to check convergence?
samples <- as.data.frame(out$samples)
# why is the length of samples only n=5000?
vacc.mcmc <- as.mcmc.list(out)
gelman.diag(out)
 

####### task 4
# trace plot
plot(samples$prob, type = "l")
library(coda)
par(mfrow=c(2, 1))
traceplot(as.mcmc(out$samples))
densplot(as.mcmc(out$samples))
# summary measures
out$summary
# Mean    Median     St.Dev.   95%CI_low 95%CI_upp
# alpha0 0.9983570 0.9988560 0.001594295 0.994153064 0.9999734
# alpha1 0.9739971 0.9807099 0.024914312 0.906153562 0.9993116
# alpha2 0.1041172 0.0829857 0.083061948 0.004197816 0.3168535


plot_posterior <- function(y, n) {
  # function to plot posterior for every row of dat
  alpha0 <- 1
  beta0 <- 1
  theta<-seq(0,1,0.01)
  likelihood<-dbeta(theta,y+1,n-y+1)
  prior<-dbeta(theta,alpha0,beta0)
  posterior<-dbeta(theta,alpha0+y,beta0+n-y)
  plot(theta,likelihood,type="l",ylim=c(0,35),col="blue",ylab="")
  lines(theta,prior,col="red", lwd=2)
  lines(theta,posterior, col="purple", lwd=2)
}

par(mfrow = c(4, 4))
for (i in 1:nrow(dat)) {
  y <- dat$yes[i]
  n <- dat$Sum[i]
  plot_posterior(y, n)
}


####### task 5
nsamples <- nrow(samples)
n <- length(dat$yes)
predprob <- matrix(0, nsamples, n)
set.seed(123)
for (i in 1:nsamples) {
  model.compiled[["prob"]] <- samples[i, "prob"]
  model.compiled[["size"]] <- samples[i, "size"]
  model.compiled$simulate(simNodes, includeData = TRUE)
  predprob[i, ] <- model.compiled[["vaccinated"]]
}  # does not produce correct results yet -> matrix has fixed values per column 
# at value of y (vaccinated / yes)


####### task 6
# posterior probability for coverage between children of different insurance status
# to get overview of empirical probabilities by insurance status
dat$percent <- dat$yes/dat$Sum

# get mean coefficient estimates from out
# alpha0
out$summary[1]
# alpha1
out$summary[2]
# alpha2
out$summary[3]
# odds for vaccination for anymediaid people vs. private insured
exp(out$summary[2])  # 2.64
# odds for vaccination for uninsured people vs. private insured
exp(out$summary[3])  # 1.11


###### task 7
####### task 2
N <- nrow(dat)

dummy_AnyMedicaid <- ifelse(dat$insurance == "Any Medicaid", 1, 0)
dummy_uninsured <- ifelse(dat$insurance == "Uninsured", 1, 0)

# 5 geography levels
levels_geography <- length(levels(dat$geography))
# vector with repeated levels of length nrow(dat)
n_vector_geography <- as.numeric(dat$geography)

model.data <- list('vaccinated' = dat$yes) # 
model.constant <- list('N' = N, 'anymed' = dummy_AnyMedicaid, 
                       'unins' = dummy_uninsured, 'total' = dat$Sum,
                       'ngeography' = levels_geography, 'vector_geo' = n_vector_geography)

model.inits <- list(alpha0=0, alpha1=1, alpha2=1)

model_set2 <- nimbleCode(
  {
    for (i in 1:N) {
      logit(p[i]) <- alpha0[vector_geo[i]] + alpha1*anymed[i] + alpha2*unins[i]
      #y_ij[i] <- alpha0 + alpha1*anymed[i] + alpha2*unins[i]
      vaccinated[i] ~ dbinom(p[i], total[i])
      #vaccinated[i] ~ dbinom()
    }
    
    for (j in 1:ngeography) {
      alpha0[j] ~ dnorm(mu, sigma^2)
    }


    # prior info (only when uniformative)
    mu ~ dnorm(0, 1)
    sigma ~ dnorm(0, 1)
    
    alpha1 ~ dbeta(shape1=1, shape2=1)
    alpha2 ~ dbeta(shape1=1, shape2=1)
    # -> how to specify prior knowledge that y is 90% or higher?

  }
)

model <- nimbleModel(model_set2, 
                     constants=model.constant, 
                     data=model.data, 
                     inits=model.inits)

model.compiled <- compileNimble(model)

out2 <- nimbleMCMC(model,
                  niter=10000,
                  nburnin=5000,
                  thin=1,
                  summary=TRUE)

out2$summary 
# summary measures for alpha1 and alpha2 unexpectedly low, check model for correctness


###### task 8
# how to obtain posterior probabilities per row from mean coefficient estimates


###### task 9
# caterpillar plot

###### task 10
# estimate vaccination coverage for anymed children in mississippi
  
  
  
  
## old:
## try with openbugs
model_set <- function()
{
  for (i in 1:N) {
    y_ij[i] <- alpha0 + alpha1*insurance[i]
    vaccinated[i] ~ dbinom(prob, size)
    
  }
  # prior info
  prob ~ dbeta(shape1=1, shape2=1)
  size ~ dbeta(shape1=1, shape2=1)
}

# causes error
write.model(model_set, "model_set.txt")
file.show("model_set.txt")