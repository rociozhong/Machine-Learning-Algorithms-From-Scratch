library(MASS)
library(glmnet) # glmnet package handles the centering and scaling automatically, 
# it will do transformation before tuning the algo, and will transform the obtained results back to the original scale
set.seed(1)
N = 400
P = 20

Beta = c(1:5/5, rep(0, P-5))
Beta0 = 0.5

# genrate X
V = matrix(0.5, P, P)
diag(V) = 1

X = as.matrix(mvrnorm(N, mu = 3*runif(P)-1, Sigma = V))

# create artifical scale of X
X = sweep(X, 2, 1:10/5, "*")

# genrate Y
y = Beta0 + X %*% Beta + rnorm(N)

# check OLS
lm(y ~ X)

# now start the Lasso 
# First we scale and center X, and record them.
# Also center y and record it. don't scale it. 
# now since both y and X are centered at 0, we don't need to worry about the intercept anymore. 
# this is because for any beta, X %*% beta will be centered at 0, so no intercept is needed. 
# However, we still need to recover the real intercept term after we are done estimating the beta. 
# The real intercept term can be recovered by using the x_center, x_scale, y2, and the beta parameter we estimated.

x_center = colMeans(X)
x_scale = apply(X, 2, sd)
X2 = scale(X)

bhat = rep(0, ncol(X2)) # initialize it
ymean = mean(y)
y2 = y - ymean

# now start to write functions 
# prepare the soft thresholding function (should be just one line, or a couple of)

soft_th <- function(b, pen)
{
  return( sign(b) * max(abs(b) - pen, 0) )
}

# initiate lambda. This is one way to do it, the logic is that I set the first lambda as the largetst gradient. 

#lambda = exp(seq(log(max(abs(cov(X2, y2)))), log(0.001), length.out = 100))
lambda = glmnet(X, y)$lambda

LassoFit <- function(myX, myY, mybeta, mylambda, tol = 1e-10, maxitr = 500, N)
{
	# initia a matrix to record the objective function value
	f = rep(0, maxitr)
	
	for (k in 1:maxitr)
	{
		# compute residual
		r = myY - myX %*% mybeta
		
		# I need to record the residual sum of squares
		f[k] = mean(r*r)
		
		for (j in 1:ncol(myX))
		{
			# add the effect of jth variable back to r 
			# so that the residual is now the residual after fitting all other variables
		  tmp_beta = mybeta
		  tmp_beta[j] = 0
			rj = myY - myX %*% tmp_beta
			
			# apply the soft thresholding function to the ols estimate of the jth variable 
			mybeta[j] = soft_th(sum(rj * myX[, j]), mylambda * N) / sum(myX[, j] * myX[, j])
		}
		if (k > 10)
		{
			# this is just my adhoc way of stoping rule, you dont have to use it
			if (sum(abs(f[(k-9):k] - mean(f[(k-9):k]))) < tol) break;
		}
	}
	return (mybeta)
}

# test your function on a large lambda (penalty) level. 
# this should produce a very spase model. 
# these are not the beta in the original scale of X

LassoFit(X2, y2, mybeta = rep(0, ncol(X2)), mylambda = lambda[10], tol = 1e-7, maxitr = 500, N = nrow(X))

# now initiate a matrix that records the fitted beta for each lambda value 

beta_all = matrix(NA, ncol(X), length(lambda))

# this vecter stores the intercept of each lambda value
beta0_all = rep(NA, length(lambda))

# initial a zero vector for bhat, 
# then throw that into the fit function using the largest lambda value. 
# that will return the fitted beta, then use this beta on the next (smaller) lambda value
# iterate until all lambda values are used

bhat = rep(0, ncol(X2)) # initialize it

for (i in 1:length(lambda)) # loop from the largest lambda value
{
	bhat = LassoFit(X2, y2, bhat, lambda[i], N = nrow(X))
	
	# data is scaled, figure out how to scale that back 
	# save the correctly scaled beta into the beta matrix 
	
	beta_all[, i] = diag(1 / x_scale) %*% bhat
	
	# recalculte the intercept term in the original, uncentered and unscaled X

	beta0_all[i] = ymean - as.numeric(x_center %*% beta_all[, i])
}

# now you have the coefficient matrix 
# each column correspond to one lambda value 
rbind("intercept" = beta0_all, beta_all)


# you should include a similar plot like this in your report
# feel free to make it look better
matplot(colSums(abs(beta_all)), t(beta_all), type="l")

# this plot should be identical (close) to your previous plot
plot(glmnet(X, y))

# set your lambda to their lambda value and rerun your algorithm 
lambda = glmnet(X, y)$lambda

# then this distance should be pretty small 
# my code gives distance no more than 0.01
max(abs(beta_all - glmnet(X, y)$beta))
max(abs(beta0_all - glmnet(X, y)$a0))

