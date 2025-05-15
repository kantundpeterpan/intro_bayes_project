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