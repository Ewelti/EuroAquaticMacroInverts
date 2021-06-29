#install.packages("brms")
library(brms)

#### two-stage models ####

#fit trend model to a single time-series
fit1 <- brm(Response ~ Year, family = poisson())

#see what the default priors are
get_prior(Response ~ Year, family = poisson())

#set priors now
prior = c(set_prior("normal(0,0.5)", class = "ar"),
          set_prior("normal(0,5)", class = "b"))

#including autocorrelation (of the residuals)
fit1 <- brm(Response ~ Year, autocor = cor_ar(~Year, p = 1),family = poisson())

#including random effects
fit1 <- brm(Response ~ Year + (1|Year), autocor = cor_ar(~Year, p = 1),family = poisson())


#combine trends in mixed model

prior1 <- prior(normal(0,10), class = b) +
  prior(cauchy(0,2), class = sd)

fit1 <- brm(trend ~ zAge + zBase * Trt + (1|patient),
            data = epilepsy, family = gaussian(), prior = prior1)

fit1 <- brm(count ~ zAge + (1 + zAge|patient),
            data = epilepsy, family = gaussian(), prior = prior1)

### one-stage model #####

get_prior(Response ~ Year, family = poisson())

#set priors now
prior = c(set_prior("normal(0,0.5)", class = "ar"),
          set_prior("normal(0,5)", class = "b"))

#including autocorrelation (of the residuals)
fit1 <- brm(Response ~ Year, autocor = cor_ar(~Year, p = 1),family = poisson())

#### model fits ####

ranef(fit1)

summary(fit1)

plot(fit1,ask = FALSE)

plot(conditional_effects(fit1),ask=FALSE)

#### model checking ####

loo(fit1)

pp_check(fit1)

WAIC(fit1)

#### model comparison ####

fit1 <- brm(rating ~ treat + period + carry,
            data = inhaler)
fit1 <- add_criterion(fit1, "waic")

fit2 <- brm(rating ~ treat + period + carry + (1|subject),
            data = inhaler)
fit2 <- add_criterion(fit2, "waic")

# compare both models
loo_compare(fit1, fit2, criterion = "waic")