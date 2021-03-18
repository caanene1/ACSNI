# Plot density first
# All realisation of e is similar
N2 <- N[N$X3 >= 0, ]
plot(density(N2$X3))

#
library(fitdistrplus)
library(logspline)

## Explore the distribution using fine and grey plot
descdist(N2$X3, discrete = FALSE)
# Observed largely beta, with potential for gamma

## Fit distributions to the data
fit.beta <- fitdist(N2$X3, "beta")
fit.gamma <- fitdist(N2$X3, "gamma")
fit.normal <- fitdist(N2$X3, "norm")

## Plot to inspect
plot(fit.beta)
plot(fit.gamma)
plot(fit.normal)

## Get the AIC of the models
fit.beta$aic
fit.gamma$aic
fit.normal$aic

## Use Kolmogorov-Smirnov test simulation to confirm the results
number <- 50000

stats <- replicate(number, {
  # simulates random variates having a specified beta distribution
  r <- rbeta(n = length(N2$X3),
             shape1 = fit.beta$estimate["shape1"],
             shape2 = fit.beta$estimate["shape2"])

  # The estimated parameters
  estfit.beta <- fitdist(r, "beta")

  # Estimate statistics
  as.numeric(ks.test(r
                     , "pbeta"
                     , shape1 = estfit.beta$estimate["shape1"]
                     , shape2 = estfit.beta$estimate["shape2"])$statistic)
})

# Plot
plot(ecdf(stats), las = 1, main = "KS-test statistic simulation (CDF)",
     col = "darkorange", lwd = 1.7)
grid()

## Get the p-value
fit <- logspline(stats)

### Get the p-value
1 - plogspline(ks.test(N2$X3, "pbeta"
                       , shape1 = fit.beta$estimate["shape1"]
                       , shape2 = fit.beta$estimate["shape2"])$statistic
               , fit
)

# library(gamlss)
# library(gamlss.dist)
# library(gamlss.add)
fit <- fitDist(N2$X3, k = 2, type = "realplus", trace = FALSE, try.gamlss = TRUE)
# summary(fit)
