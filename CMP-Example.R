
# good talk by Galit Shmuéli on the CMP:
# https://youtu.be/Nt7gmKmxjxA
library(COMPoissonReg)
library(degreenet)
library(compoisson)
library(glmmTMB)
y <- simcmp(n=100, v=c(1,sqrt(5)), maxdeg=10000) - 1
barplot(table(y))
mean(y)
var(y)


fit <- glmmTMB(y ~ 1, family = "compois")
summary(fit)
exp(fixef(fit)$cond)

sigma(fit)

names(fit)
fit$fit

fit2 <- glm.cmp(y ~ 1)
coef(fit2)
predict(fit2)
nu(fit2)[1]


mean(y^nu(fit2))

simulate(fit)[[1]]



