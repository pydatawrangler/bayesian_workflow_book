# generate fake data
N <- 100
Y <- rnorm(N, 1.6, 0.2)
hist(Y)

# compile model
library(rstan)

model <- stan_model('first_model.stan')

# pass data to stan and run model
options(mc.cores=4)
fit <- sampling(model, list(N=N, Y=Y), iter=200, chains=4)

# diagnose
print(fit)

# graph
params <- extract(fit)
hist(params$sigma)

library(shinystan)
launch_shinystan(fit)
