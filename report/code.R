library(tidyr)
library(dplyr)
library(ggplot2)
library(kableExtra)
library(gridExtra)
library(rjags)
library(R2jags)
library(bayesplot)
library(glue)
library(latex2exp)
library(stringr)
library(ggridges)

source("../functions.R")

# Load data
vaccination_data <- data.frame(
  # Geography = c("North Carolina", "North Carolina", "North Carolina", 
  #               "Georgia", "Georgia", "Georgia",
  #               "Wisconsin", "Wisconsin", "Wisconsin",
  #               "Florida", "Florida", "Florida",
  #               "Mississippi", "Mississippi"),
    Geography = c("NC", "NC", "NC", 
                "GA", "GA", "GA",
                "WI", "WI", "WI",
                "FL", "FL", "FL",
                "MS", "MS"),
  # Insurance = c("Any Medicaid", "Private Insurance Only", "Uninsured",
  #               "Any Medicaid", "Private Insurance Only", "Uninsured",
  #               "Any Medicaid", "Private Insurance Only", "Uninsured",
  #               "Any Medicaid", "Private Insurance Only", "Uninsured",
  #               "Private Insurance Only", "Uninsured"),
  Insurance = c("Medicaid", "Private", "Uninsured",
                "Medicaid", "Private", "Uninsured",
                "Medicaid", "Private", "Uninsured",
                "Medicaid", "Private", "Uninsured",
                "Private", "Uninsured"),
  Vaccinated = c(380, 632, 28, 363, 527, 36, 282, 514, 16, 446, 588, 28, 400, 27),
  SampleSize = c(419, 673, 34, 396, 576, 50, 332, 548, 34, 490, 628, 39, 441, 32)
)

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
    alpha_above90, beta_above90, vaccination_data$Vaccinated, vaccination_data$SampleSize
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
    relevel(ref = "Private") -> vaccination_data$Insurance

vaccination_data$Geography %>% 
    as.factor() %>%
    relevel(ref = "GA") -> vaccination_data$Geography


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
    logit(p[t]) <- alpha_0 + alpha_1 * x_1[t] + alpha_2 * x_2[t]
 }

  # Priors
  alpha_0 ~ dnorm(0.0,0.01)
  alpha_1 ~ dnorm(0.0,0.01)
  alpha_2 ~ dnorm(0.0,0.01)


  # Parameters of interest
  
  ## Question 2
  pi_private <- exp(alpha_0)/(1+exp(alpha_0))
  pi_medicaid <- exp(alpha_0 + alpha_1)/(1+exp(alpha_0 + alpha_1))
  pi_uninsured <- exp(alpha_0 + alpha_2)/(1+exp(alpha_0 + alpha_2))

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
    "alpha_0", "alpha_1", "alpha_2",
    "pi_private", "pi_medicaid", "pi_uninsured",
    "diff_priv_medicaid", "diff_priv_uninsured"
    )

# Run the model
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(jags_str),
  n.chains = 2,
  n.iter = 40000,
  n.burnin = 2000,
  n.thin = 2
)

alpha_params <- c(
    "alpha_0", "alpha_1", "alpha_2"
)

pi_params <- c(
    "pi_private", "pi_medicaid",
    "pi_uninsured"
)

mcmc <- as.mcmc(model_run)

mcmcdf <- mcmc.as.data.frame(mcmc)

alpha_plot_df <- mcmcdf %>% 
    select(alpha_params, "chain") %>%
    pivot_longer(
        cols = all_of(alpha_params),
        names_to = "parameter",
        values_to = "value"
    )

dens_plot_df <- mcmcdf %>% 
    select(pi_params, "chain") %>%
    pivot_longer(
        cols = all_of(pi_params),
        names_to = "parameter",
        values_to = "value"
    )


# Question 4

alpha_plot_df %>%
    group_by(parameter) %>%
    summarize(
      mean = mean(value),
      median = median(value),
      # mode = ?,
      sd = sd(value)
    ) -> alpha_params_summary

alpha_hpddf <- HPDinterval(
  as.mcmc(
    do.call(
      rbind, mcmc[,alpha_params]
    )
  )
) %>% as.data.frame()
alpha_hpddf$parameter <- rownames(alpha_hpddf)
rownames(alpha_hpddf) <- NULL

alpha_summary <- merge(
  alpha_params_summary, alpha_hpddf,
  by = 'parameter'
  )

rownames(alpha_summary) <- alpha_summary$parameter
# rownames(alpha_summary) <- c("$\alpha_0$", "$\alpha_1$", "$\\alpha_2$")
alpha_summary$parameter <- NULL
colnames(alpha_summary) <- c('Mean', 'Median', 'SD', 'LL', 'UL')

alpha_summary %>% round(3) -> alpha_summary

# Question 5

dens_plot_df %>%
    group_by(parameter) %>%
    summarize(
      mean = mean(value),
      median = median(value),
      # mode = ?,
      sd = sd(value)
    ) -> pi_params_summary

hpddf <- HPDinterval(
  as.mcmc(
    do.call(
      rbind, mcmc[,pi_params]
    )
  )
) %>% as.data.frame()
hpddf$parameter <- rownames(hpddf)
rownames(hpddf) <- NULL

final_summary <- merge(
  pi_params_summary, hpddf,
  by = 'parameter'
  )

rownames(final_summary) <- final_summary$parameter
final_summary$parameter <- NULL
colnames(final_summary) <- c('Mean', 'Median', 'SD', 'LL', 'UL')

final_summary %>% round(3) -> final_summary

vaccination_data$link <- c(rep(c('pi_medicaid', 'pi_private', 'pi_uninsured'), 4),
                           'pi_private', 'pi_uninsured')
final_summary$link <- rownames(final_summary)


q5 <- merge(vaccination_data, final_summary, by='link') %>% 
  merge(.,post_summary, by=c('Geography', 'Insurance')) %>%
  select(colnames(vaccination_data),-link, Mean, post_mean, post_mean.1) 

q5$diff_logreg_beta11 <- q5$Mean - q5$post_mean
q5$diff_logreg_betaconf <- q5$Mean - q5$post_mean.1

q5 <- q5 %>%
  add_row(Geography="MS", Insurance="Medicaid",
   diff_logreg_beta11 = NULL, diff_logreg_betaconf = NULL,
   Mean=final_summary['pi_medicaid', 'Mean']) %>%
   arrange(Geography)

q5_compbeta11_plot <- q5 %>% 
  select(Geography, Insurance, diff_logreg_beta11) %>% 
  add_row(Geography="MS", Insurance="Medicaid", diff_logreg_beta11 = NULL)  %>% 
  ggplot(aes(x = Geography, y = diff_logreg_beta11, color = Insurance, fill = Insurance)) + 
  geom_bar(stat = 'identity', position='dodge') +
  theme_minimal() + theme(legend.position = "none")

q5_compbetaconf_plot <- q5 %>% 
  select(Geography, Insurance, diff_logreg_betaconf) %>% 
  add_row(Geography="MS", Insurance="Medicaid", diff_logreg_betaconf = NULL)  %>% 
  ggplot(aes(x = Geography, y = diff_logreg_betaconf, color = Insurance, fill = Insurance)) + 
  geom_bar(stat = 'identity', position='dodge') +
  theme_minimal() + theme(legend.position = "top")

cols <- c(
  "Geo.", "Ins.", "$k$", "$n$", 
  "$\\bar\\pi_{logreg}$",
  "$\\bar\\pi_{Beta(1,1)}$",
  sprintf("$\\bar\\pi_{Beta(%s,%s)}$", alpha_above90, beta_above90),
  "$\\Delta_{Beta(1,1)}$",
  sprintf("$\\Delta_{Beta(%s,%s)}$", alpha_above90, beta_above90)
)
## Get greek letter labels for 
## facet labels
## https://ggplot2.tidyverse.org/reference/labellers.html
colnames(q5) <- cols
tex_labels <- list()

for(l in cols[8:9]){
  tex_labels[l] = TeX(l)
}

q5_long <- q5 %>% select(c(1,2,8,9)) %>% 
  pivot_longer(cols = c(3,4))

q5_long$name <- tex_labels[q5_long$name] %>% unlist() %>% as.character() 

q5_compplot <-  q5_long %>% 
  ggplot(aes(x = Geo., y = value, fill = Ins.)) + 
  geom_bar(stat='identity', position='dodge') + 
  facet_wrap(~name, labeller = label_parsed) + 
  theme_minimal() + theme(legend.position = "top")

fig_q5_comp_cap = paste("Differences between posterior means of vaccine coverage per region and insurance status obtained from logistic regression and conjugate pair modeling.
Differences are most pronounced for the Uninsured group in MS, NC, WC  when using a $Beta(1,1) prior$ and all states when using a ", sprintf("$Beta(%s, %s)$", alpha_above90, beta_above90), "prior.")

# Question 6
## Density plots and eCDF of differences

diff_plot_df <- mcmcdf %>% 
  select(diff_priv_medicaid, diff_priv_uninsured) %>%
  pivot_longer(
    cols=c(diff_priv_medicaid, diff_priv_uninsured),
    names_to='parameter',
    values_to='value'
  )

diff_dens_plot <- ggplot(diff_plot_df, aes(x = value)) +
  geom_density(aes(y = after_stat(density)), color = "blue", fill = "blue", alpha = 0.3) +
  geom_line(stat = "ecdf", aes(y = after_stat(y) * max(density(x)$y)), color = "red") +
  scale_y_continuous(
    name = "Density",
    sec.axis = sec_axis(~./max(density(diff_plot_df$value)$y), name = "ECDF") # This assumes max density is similar across facets, might need adjustment
  ) +
  facet_wrap(~ parameter, scales = "free") +
  theme_classic() +
  labs(x = expression("Delta"), title = "Posterior densities and ECDF of differences in vaccine coverages between insurance groups")

## Probabilites
## P(pi_priv > pi_medicaid) = P(pi_priv - pi_medicaid > 0)
p_priv_greater_medicaid <- sum(mcmcdf$diff_priv_medicaid>0) / dim(mcmcdf)[1]

## P(pi_priv > pi_uninsured) = P(pi_priv - pi_uninsured > 0)
p_priv_greater_uninsured <- sum(mcmcdf$diff_priv_uninsured>0) / dim(mcmcdf)[1]


# Question 7
## Include region in the logistic regression model

q7_jags_str <- "
model
{
  # Likelihood
  for (t in 1:T) {
    y[t] ~ dbin(p[t], K[t])
    logit(p[t]) <- inprod(beta, X[t,])
 }

  # Priors
  for (i in 1:7){
    beta[i] ~ dnorm(0.0,0.01)
  }

  # Vaccine coverage
  for (t in 1:T){
    #logodds
    lo[t] <- inprod(beta, X[t,])
    #convert to probability / vaccine coverage
    coverage[t] <- exp(lo[t]) / (1 + exp(lo[t]))
  }
}
"

q7_X <- model.matrix(
    ~ Geography + Insurance,
     data = vaccination_data
)

# Set up the data
q7_model_data <- list(
    T = dim(vaccination_data)[1],
    y = vaccination_data$Vaccinated,
    X = q7_X,
    K = vaccination_data$SampleSize)


# Choose the parameters to watch
q7_model_parameters <- c(
    "beta", "lo", "coverage"
    )

# Run the model
q7_model_run <- jags(
  data = q7_model_data,
  parameters.to.save = q7_model_parameters,
  model.file = textConnection(q7_jags_str),
  n.chains = 2,
  n.iter = 40000,
  n.burnin = 2000,
  n.thin = 2
)

# convert output to data.frame
q7_mcmc <- as.mcmc(q7_model_run)
q7_mcmcdf <- mcmc.as.data.frame(q7_mcmc)

# extract indices corresponding to coverage parameters
pi_idx <- grepl("coverage", colnames(q7_mcmcdf))

# use index to extract columns
q7_mcmcdf_coverage <- q7_mcmcdf[,pi_idx]
q7_mcmcdf_coverage$chain <- q7_mcmcdf$chain

# columns are string ordered ("coverage[13]") etc ...
map_mcmc_vacc_data <- str_match(colnames(q7_mcmcdf)[pi_idx], '\\d+') %>% as.numeric()

q7_pi_long_df <- q7_mcmcdf_coverage %>%                  
    select(map_mcmc_vacc_data, "chain") %>%
    pivot_longer(
        cols = all_of(map_mcmc_vacc_data),
        names_to = "parameter",
        values_to = "value"
    )

# compute posterior summary measures
q7_pi_summary <- q7_pi_long_df %>% 
  select(-chain) %>%
  group_by(parameter) %>%
  summarize(
      mean = mean(value),
      median = median(value),
      # mode = ?,
      sd = sd(value)
    )

# compute HPDs
q7_pi_hpddf <- HPDinterval(
  as.mcmc(
    do.call(
      rbind, q7_mcmc[,colnames(q7_mcmcdf)[pi_idx]]
    )
  )
) %>% as.data.frame()
q7_pi_hpddf$parameter <- rownames(q7_pi_hpddf)

# create final table as in Q5
q7_post_summary <- merge(
  q7_pi_summary, q7_pi_hpddf, by = "parameter"
  )
q7_post_summary$parameter <- NULL
colnames(q7_post_summary) <- c('Mean', 'Median', 'SD', 'LL', 'UL')
q7_post_summary$Geography <- vaccination_data$Geography[map_mcmc_vacc_data]
q7_post_summary$Insurance <- vaccination_data$Insurance[map_mcmc_vacc_data]

q7 <- merge(vaccination_data, q7_post_summary, by = c('Geography', 'Insurance')) %>% 
merge(., post_summary, by=c('Geography', 'Insurance')) %>%
select(colnames(vaccination_data),-link, Mean, post_mean, post_mean.1)

q7$diff_logreg_beta11 <- q7$Mean - q7$post_mean
q7$diff_logreg_betaconf <- q7$Mean - q7$post_mean.1

# compute posterior for coverage MS/Medicaid using additive model
ms_medi_contrast <- c(1,0,1,0,0,1,0)
post_coverage_ms_medi <- (ms_medi_contrast %*% t(q7_mcmcdf[,1:7])) %>% {exp(.) / (1 + exp(.)) }

q7 <- q7 %>%
  add_row(Geography="MS", Insurance="Medicaid",
   diff_logreg_beta11 = NULL, diff_logreg_betaconf = NULL,
   Mean=mean(post_coverage_ms_medi)) %>%
   arrange(Geography)

colnames(q7) <- cols
tex_labels <- list()

for(l in cols[8:9]){
  tex_labels[l] = TeX(l)
}

q7_long <- q7 %>% select(c(1,2,8,9)) %>% 
  pivot_longer(cols = c(3,4))

q7_long$name <- tex_labels[q7_long$name] %>% unlist() %>% as.character() 

q7_compplot <-  q7_long %>% 
  ggplot(aes(x = Geo., y = value, fill = Ins.)) + 
  geom_bar(stat='identity', position='dodge') + 
  facet_wrap(~name, labeller = label_parsed) + 
  theme_minimal() + theme(legend.position = "top")

# Q8
#replace coverage[x] by Geo_insurance
colnames(q7_mcmcdf_coverage) <- c(vaccination_data[map_mcmc_vacc_data,] %>% mutate(key = paste(Geography, "_", Insurance, sep = '')) %>% select(key) %>% as.vector() %>% unlist(), 'chain')

# reference df
nc <- q7_mcmcdf_coverage[,grepl("NC", colnames(q7_mcmcdf_coverage))]
colnames(nc) <- c("Medicaid", "Private", "Uninsured")

# compute ratios 
ratio_list <- list()

for (col in colnames(q7_mcmcdf_coverage[,!grepl("NC", colnames(q7_mcmcdf_coverage))])){
  print(col)
  if (col != 'chain'){
    if (grepl("Medicaid", col)){
      ratio_list[[col]] <- (q7_mcmcdf_coverage[,col] / nc["Medicaid"]) %>% unlist() %>% as.vector()
    } else if (grepl("Private", col)){
      ratio_list[[col]] <- (q7_mcmcdf_coverage[,col] / nc["Private"]) %>% unlist() %>% as.vector()
    } else if (grepl("Uninsured", col)){
      ratio_list[[col]] <- (q7_mcmcdf_coverage[,col] / nc["Uninsured"]) %>% unlist() %>% as.vector()
    }
  }
}
ratio_df <- as.data.frame(ratio_list)

q8_hpddf <- ratio_df %>% as.mcmc() %>% HPDinterval()

ratio_df_long <- ratio_df %>% 
  pivot_longer(
    cols = all_of(colnames(ratio_df)),
    names_to = 'parameter',
    values_to = 'value'
  )

ridge_plot <- ratio_df_long %>%
  ggplot(aes(x = value, y = parameter, fill = parameter)) + 
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
