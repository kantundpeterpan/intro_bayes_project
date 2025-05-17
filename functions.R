beta_binomial_posterior <- function(
    prior_a, prior_b, k, n
    ){
        post_a = prior_a + k
        post_b = prior_b + n - k
        post_mean = post_a / (post_a + post_b)
        post_mode = (post_a - 1) / (post_a + post_b - 2)
        post_sd = sqrt((post_a*post_b)/((post_a + post_b)**2*(post_a+post_b+1)))
        # 95% equal tail interval
        hpd_ll = qbeta(0.025, post_a, post_b)
        hpd_ul = qbeta(0.975, post_a, post_b)
        results <- data.frame(
            post_a = post_a,
            post_b = post_b,
            post_mean = post_mean,
            post_mode = post_mode,
            post_sd = post_sd,
            hpd_ll = hpd_ll,
            hpd_ul = hpd_ul
        )

        return (results)
    }

plot_two_posteriors <- function(
    a1, b1, a2, b2
){
    x <- seq(0,1,0.001)

    # Calculate y values for both distributions
    y1 <- dbeta(x, a1, b1)
    y2 <- -dbeta(x, a2, b2) # Keeping the negative as in your original code

    # Create a data frame with the x, y values, and an identifier for each distribution
    plot_data <- data.frame(
        x = c(x, x), # Combine x values
        density = c(y1, y2), # Combine y values
        dist = factor(rep(c("Beta(1,1)", sprintf("Beta(%s,%s)", alpha_above90, beta_above90) ), each = length(x))) # Create identifier
    )

    # Use the combined data frame and map color to the 'dist' variable
    p <- ggplot(plot_data, aes(x = x, y = density, color = dist)) +
         geom_line() + theme(legend.position = "none") # Now only one geom_line call is needed
        #  scale_color_brewer(palette = "Set1") +
        #  labs(color = "Distribution") # Add a meaningful legend title

    return(p)
}

mcmc.as.data.frame <- function(mcmc){
    df = as.data.frame(as.matrix(mcmc))
    df$chain <- rep(1:nchain(mcmc), each = nrow(mcmc[[1]]))

    return (df)
}