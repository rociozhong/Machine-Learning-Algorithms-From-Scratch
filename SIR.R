# Write my own code to fit the sliced inverse regression, validate it by
# comparing to the “dr” package. Use use 10 as the number of slices.


rm(list = ls())

# response variable y 
# indep variables: X: n*p
# H slices

fun_sir = function(X, y, H){
   n = length(y)
   # Normalize X
   X = scale(X, center=T, scale=F)
   e_tmp = eigen(cov(X))
   sqrt_sigma = e_tmp$vector %*% diag(1/sqrt(e_tmp$value)) %*% t(e_tmp$vector)
   X = X %*% sqrt_sigma
   
   # sort the data by y
   dta = as.data.frame(cbind(y, X))
   sort_dta = dta[order(dta[, 1]), ]
   # split the data into H slices
   slice_dta = split(sort_dta, factor(sort(rank(row.names(sort_dta)) %% H)))
   # calculate mean for each y slice 
   smean_y = as.vector(unlist(lapply(slice_dta, function(x) mean(x[ ,1]))))
   # calculate mean vector (including all vector means for each x slice)
   smean_x = lapply(slice_dta, function(x) apply(x[, 2: ncol(x)], 2, mean))
   
   
   # split the vector into chunks of size p  = ncol(X)
   # p = ncol(X)
   # split(smean_x, ceiling(seq_along(smean_x) / p))
   # sample mean of x_bar
   x_bar = apply(X, 2, sum) / n 
   
   # Compute the covariance matrix for the slice means of x, weighted by the slice sizes
   temp_cov = lapply(smean_x, function(x)  length(x) * (x - x_bar) %*% t(x - x_bar))
   covm_slice_mean = Reduce('+', temp_cov) / n
   
   
   # compute the sample covariance for x_i's
   # temp_sample_cov = apply(X, 1, function(x) (x - x_bar))
   # mat = list()
   # for(i in 1: ncol(temp_sample_cov)){mat[[i]] = temp_sample_cov[, i] %*% t(temp_sample_cov[, i])}
   # sample_covx = Reduce('+', mat) / n
   # sample_covx = cov(X)
   # A = solve(sample_covx) %*% covm_slice_mean
   A = covm_slice_mean
   ev = eigen(A)
   result = x %*% ev$vectors[, c(1:2)]
   return(list(ev = ev, result = result))
}

set.seed(21)
n = 1000; p = 10
x = matrix(rnorm(n*p), n, p)
b = matrix(c(1, 1, rep(0, p-2)))
y = sin(x %*% b) + rnorm(n)

final_re = fun_sir(x, y, 10)

library(MASS)
library(dr)
fit.sir = dr(y~., data = data.frame(x, y), method = "sir", nslices = 10)


# underlying model that cannot be detected by SIR
# lets do a new function 
set.seed(111)
n = 1000
p = 10

x = matrix(rnorm(n*p), n, p)

b1 = matrix(c(1, 0, 1, rep(0, p-3)))
b2 = matrix(c(0, 1, 0, 1, 1, rep(0, p-5)))

dir1 = x %*% b1
dir2 = x %*% b2

link = 5*cos(dir1)^2 - (dir2)^3
y = link + rnorm(n)

fun_sir(x, y, 10)$ev$vector
fit_sir2 = dr(y ~., data = data.frame(x, y), method = "sir", nslices = 10)










