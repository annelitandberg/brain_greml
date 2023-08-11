library(gaston)

n_sample = 1000 # amount of participants in sample
n_site = 10000 # amount of vertices

B = matrix(rnorm(n_sample * n_site), n_sample, n_site) # generating data with 10000 brain regions per individual (n = 1000)
B = scale(B) # standardizing each column
BB = tcrossprod(B) / n_site # creating a cortical thickness relatedness matrix
b = B %*% (sqrt(1 / n_site) * rnorm(n_site)) # sampling vertex-wise effects

res = rnorm(n_sample) # residuals
y = b + res # phenotype

rs = list(B = BB) # list of relatedness matrices
X = matrix(1, length(y), 1) # covariates with fixed effects
fit = lmm.aireml(y, X, rs, min_tau = c(0)) # fitting model (variance decomposition)
fit$tau # total vertex-wise effect
fit$sigma2 # residual effect
fit$BLUP_beta # fixed effect
