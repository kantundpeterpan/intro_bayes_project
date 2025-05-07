library(tidyr)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(rjags)
library(R2jags)

source("../functions.R")

# Load data
vaccination_data <- data.frame(
  Geography = c("North Carolina", "North Carolina", "North Carolina", 
                "Georgia", "Georgia", "Georgia",
                "Wisconsin", "Wisconsin", "Wisconsin",
                "Florida", "Florida", "Florida",
                "Mississippi", "Mississippi"),
  Insurance = c("Any Medicaid", "Private Insurance Only", "Uninsured",
                "Any Medicaid", "Private Insurance Only", "Uninsured",
                "Any Medicaid", "Private Insurance Only", "Uninsured",
                "Any Medicaid", "Private Insurance Only", "Uninsured",
                "Private Insurance Only", "Uninsured"),
  Vaccinated = c(380, 632, 28, 363, 527, 36, 282, 514, 16, 446, 588, 28, 400, 27),
  SampleSize = c(419, 673, 34, 396, 576, 50, 332, 548, 34, 490, 628, 39, 441, 32)
)

beta_binomial_posterior <- function(
    prior_a, prior_b, k, n
    ){
        post_a = prior_a + k
        post_b = prior_b + n - k
        post_mean = post_a / (post_a + post_b)
        post_mode = (post_a - 1) / (post_a + post_b - 2)
        post_variance = (post_a*post_b)/(post_a + post_b)**2*(post_a+post_b+1)
        # HPD interval
        hpd_ll = qbeta(0.025, post_a, post_b)
        hpd_ul = qbeta(0.975, post_a, post_b)
        results <- data.frame(
            post_a = post_a,
            post_b = post_b,
            post_mean = post_mean,
            post_mode = post_mode,
            post_variance = post_variance,
            hpd_ll = hpd_ll,
            hpd_ul = hpd_ul
        )

        return (results)
    }

plot_two_posteriors <- function(
    a1, b1, a2, b2
){
    x <- seq(0,1,0.001)
    p <- ggplot() + 
         geom_line(aes(x = x, y = dbeta(x, a1,b1))) + 
         geom_line(aes(x = x, y = -dbeta(x, a2, b2)))

    return(p)
}

mcmc.as.data.frame <- function(mcmc){
    df = as.data.frame(as.matrix(mcmc))
    df$chain <- rep(1:nchain(mcmc), each = nrow(mcmc[[1]]))

    return (df)
}

alpha_above90 = 150
beta_above90 = 7
x = seq(0, 1, 0.001)


df <- data.frame(
    x = c(x,x),
    knowledge = c(rep("none", length(x)), rep("above90", length(x))),
    beta = c(dbeta(x, 1, 1), dbeta(x, alpha_above90, beta_above90))
)


plot_priors <- ggplot(df, aes(x = x, y = beta, color = knowledge)) + geom_line()

uni_posterior <- beta_binomial_posterior(
    1,1, vaccination_data$Vaccinated, vaccination_data$SampleSize
)

# uni_posterior <- cbind(
#     vaccination_data[,c('Geography', 'Insurance')],
#     uni_posterior
# )

confident_posterior <- beta_binomial_posterior(
    150,7, vaccination_data$Vaccinated, vaccination_data$SampleSize
)

post_summary <- cbind(
    cbind(
    vaccination_data[,c("Geography", "Insurance")],
    uni_posterior),
    confident_posterior
)

posteriors <- data.frame(
    a1 = uni_posterior$post_a,
    b1 = uni_posterior$post_b,
    a2 = confident_posterior$post_a,
    b2 = confident_posterior$post_b
)

plots <- apply(posteriors, MARGIN=1, FUN=function(row){
    plot_two_posteriors(
        row[1], row[2], row[3], row[4]
    )
    # sum(row)
})


## Question 2 - Bayesian Logistic regression model
#set private insurance as reference level
vaccination_data$Insurance %>% 
    as.factor() %>%
    relevel(ref = "Private Insurance Only") -> vaccination_data$Insurance

vaccination_data$Geography %>% 
    as.factor() %>%
    relevel(ref = "Georgia") -> vaccination_data$Geography


X <- model.matrix(
    ~ Insurance,
     data = vaccination_data
)

jags_str <- "
model
{
  # Likelihood
  for (t in 1:T) {
    y[t] ~ dbin(p[t], K[t])
    logit(p[t]) <- alpha_0 + ins_1 * x_1[t] + ins_2 * x_2[t]
 }

  # Priors
  alpha_0 ~ dnorm(0.0,0.01)
  ins_1 ~ dnorm(0.0,0.01)
  ins_2 ~ dnorm(0.0,0.01)


  # Parameters of interest
  
  ## Question 2
  pi_private <- exp(alpha_0)/(1+exp(alpha_0))
  pi_medicaid <- exp(alpha_0 + ins_1)/(1+exp(alpha_0 + ins_1))
  pi_uninsured <- exp(alpha_0 + ins_2)/(1+exp(alpha_0 + ins_2))

  ## Question 6
  diff_priv_medicaid  <- pi_private - pi_medicaid
  diff_priv_uninsured <- pi_private - pi_uninsured
}
"

# Set up the data
model_data <- list(
    T = dim(vaccination_data)[1],
    y = vaccination_data$Vaccinated,
    x_1 = X[,2],
    x_2 = X[,3],
    # x_3 = X[,4],
    # x_4 = X[,5],
    K = vaccination_data$SampleSize)

# Choose the parameters to watch
model_parameters <- c(
    # "alpha_0",
    # # "geom_1", "geom_2",
    # "ins_1", "ins_2"
    "pi_private", "pi_medicaid", "pi_uninsured",
    "diff_priv_medicaid", "diff_priv_uninsured"
    )

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(jags_str),
  n.chains = 4,
  n.iter = 40000,
  n.burnin = 2000,
  n.thin = 2
)

pi_params <- c(
    "pi_private", "pi_medicaid",
    "pi_uninsured"
)

mcmc <- as.mcmc(model_run)

mcmcdf <- mcmc.as.data.frame(mcmc)

dens_plot_df <- mcmcdf %>% 
    select(pi_params, "chain") %>%
    pivot_longer(
        cols = all_of(pi_params),
        names_to = "parameter",
        values_to = "value"
    )
