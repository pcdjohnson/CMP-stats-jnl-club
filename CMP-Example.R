#devtools::install_github("jreduardo/cmpreg")
#devtools::install_github("cran/compoisson") # archive version - package no longer on CRAN
rm(list = ls())
library(COMPoissonReg) # glm.cmp
library(degreenet) # simcmp
library(MASS) # glm.nb
library(scales) # alpha
#library(compoisson) 
library(glmmTMB)
#library(cmpreg)
# Other options:
#   try mpcmp: Mean-Parametrized Conway-Maxwell Poisson (COM-Poisson) Regression?
#   There’s an experimental CMP in brms (a Stan frontend) https://github.com/paul-buerkner/brms/issues/607


# Z function
Z <- function(lambda, nu, stop.at = 100) {
  i <- 0:stop.at
  sum(lambda^i / (factorial(i)^nu))
}

Z(lambda = 1, nu = 1, stop.at = 5)
Z(lambda = 1, nu = 1, stop.at = 10)
Z(lambda = 1, nu = 1, stop.at = 100)
Z(1, 1)
exp(1)


# Counts in biology are often overdispersed so a poisson GLM is a poor default 

# Simulate some overdispersed (neg binom) counts
n <- 100
x <- gl(2, n/2)
set.seed(8723645)
y.nb <- rnbinom(n, mu = 5, size = 0.1)

# Fit the neg binom GLM
fit.nb <- glm.nb(y.nb ~ x)
summary(fit.nb)

# Fit the Poisson GLM
fit.pois.nb <- glm(y.nb ~ x, family = poisson)
summary(fit.pois.nb)

# Simulate from the Poisson GLM
y.sim.pois.nb <- simulate(fit.pois.nb)[[1]]

# Show that Poisson greatly underestimates the amount of noise
par(mfrow = c(1, 2))
col <- alpha("blue", 0.3)
stripchart(y.nb ~ x, vertical = TRUE, method = "jitter", 
           main = "What the data look like",
           pch = 16, col = col)
legend("topright", legend = paste0("P=", round(coef(summary(fit.nb))["x2", "Pr(>|z|)"], 3)))
stripchart(y.sim.pois.nb ~ x, vertical = TRUE, method = "jitter", 
           main = "What the data look like\nto a Poisson GLM",
           pch = 16, col = col)
legend("topright", legend = paste0("P=", round(coef(summary(fit.pois.nb))["x2", "Pr(>|z|)"], 3)))


# Simulate underdispersed data where both the mean and dispersion differ between groups
y.cmp <- 
  c(simcmp(n = n/2, v = c(5, 0.5), maxdeg = 10000), 
    simcmp(n = n/2, v = c(5.5, 2), maxdeg = 10000)) - 1

lapply(levels(x), function(i) barplot(table(y.cmp[x == i]), main = paste0("x=", i)))
tapply(y.cmp, x, mean)
tapply(y.cmp, x, var)
tapply(y.cmp, x, sd)


# Fit CMP GLM
fit.cmp <- glm.cmp(y.cmp ~ x, formula.nu = ~ x)
fit.cmp
exp(coef(fit.cmp))
nu(fit.cmp)
predict(fit.cmp) # uses eqn 7 (bad) approximation in Lynch et al.

# What does the Poisson GLM make of the underdispersed counts?
fit.pois.cmp <- glm(y.cmp ~ x, family = poisson)
summary(fit.pois.cmp)

# Simulate from the CMP GLM
y.sim.pois.cmp <- simulate(fit.pois.cmp)[[1]]

# Show how the Poisson GLM overestimates the amount of noise
par(mfrow = c(1, 2))
stripchart(y.cmp ~ x, vertical = TRUE, method = "jitter", col = col, pch = 16,
           main = "What the data look like")
legend("topleft", 
       legend = 
         c(
           paste0("lambda P=", round(as.numeric(summary(fit.cmp)$DF["X:x2", "p.value"]), 3)),
           paste0("nu P=", round(as.numeric(summary(fit.cmp)$DF["S:x2", "p.value"]), 3))))
stripchart(y.sim.pois.cmp ~ x, vertical = TRUE, method = "jitter", col = col, pch = 16,
           main = "What the data look like\nto a Poisson GLM")
legend("topright", legend = paste0("P=", round(coef(summary(fit.pois.cmp))["x2", "Pr(>|z|)"], 3)))


# Try glmmTMB (I didn't explore this option much as for some 
# data sets it took 10+ min)
system.time(
  fit.tmb <- glmmTMB(y.cmp ~ x, dispformula = ~ x, family = "compois")
)
summary(fit.tmb) # no sig. diff. in conditional model - differs from glm.cmp - why?
exp(fixef(fit.tmb)$cond) # gives mean and multiplicative x effect
exp(fixef(fit.tmb)$disp) # I haven't worked out the parameterisation for glmmTMB



