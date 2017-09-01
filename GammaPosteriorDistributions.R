#   Simulates the Bayesian updating process when
#   the prior is Gamma distributed with Gamma(a, b) and 
#   the data are poisson distributed.
#   
#   Input: The parameters in the gamma distribution for the prior and the data points.
#
#   Output: Plot of the posterior distributions for each step in the updating process.

GammaPosteriorDistributions <- function(gamma.a, gamma.b, data) {
  
  simulation.number <- 1000
  n <- length(data)
  theta.values <- seq(min(data), max(data), length.out = simulation.number)
  
  prior <- dgamma(theta.values, gamma.a, gamma.b)
  posterior <- matrix(NA, simulation.number, n)
  
  for (j in 1:n) {
    posterior[, j] <- dgamma(theta.values, gamma.a + cumsum(data)[j], gamma.b + j)
  }
  
  # Prior and last posterior should be black and solid
  colours <- c(1, seq(2, (length(data) - 1)), 1)
  linetypes <- c(1, rep(2, length(data) - 1), 1)
  par(bty = "l")
  par(cex.lab = 1.4142)
  
  matplot(theta.values, cbind(prior, posterior),
          type = "l",
          lty = linetypes,
          col = colours,
          xlab = "",
          ylab = "Posterior distribution",
          main = paste("Prior: ", "Gamma(", gamma.a, ", ", gamma.b, ")", sep = ""))
  
  posterior.mean <- c()
  for (i in seq(length(data))) {
    calculate.mean <- (gamma.a + cumsum(data)[i]) / (gamma.b + i)
    posterior.mean <- c(posterior.mean, calculate.mean)
  }
  
  names(posterior.mean) <- seq(length(data))
  
  posterior.sd <- c()
  for (i in seq(length(data))) {
    calculate.sd <- sqrt((gamma.a + cumsum(data)[i]) / (gamma.b + i)^2)
    posterior.sd <- c(posterior.sd, calculate.sd)
  }
  names(posterior.sd) <- seq(length(data))
  
  prior.mean <- gamma.a / gamma.b
  prior.sd <- sqrt(gamma.a / gamma.b^2)
  
  out <- list(posterior.mean = posterior.mean, posterior.sd = posterior.sd, prior.mean = prior.mean, prior.sd = prior.sd)
  return(out)
}