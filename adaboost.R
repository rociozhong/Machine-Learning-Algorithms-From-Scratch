
## write my own code for a one-dimensional adaboost using a stump model as the weak learner.
rm(list = ls())

# a stump model, with weights
f = function(x, y, w){
  score = -Inf
  cut = -Inf
  fL = -Inf
  fR = -Inf
  for (c in x) {
    left = (x <= c)
    right = (x > c)
    pLeft = 0
    pRight = 0
    if (sum(left) > 0) {
      pLeft = sum(w[left][y[left] == 1]) / sum(w[left]) 
    }
    if (sum(right) > 0) {
      pRight = sum(w[right][y[right] == 1]) / sum(w[right]) 
    }
    giniLeft = pLeft * (1 - pLeft)
    giniRight = pRight * (1 - pRight)
    tmp = ((-giniLeft * sum(w[left])) + (-giniRight * sum(w[right]))) / sum(w)
    if (tmp > score) {
      score = tmp
      cut = c
      fL = 2 * as.numeric((pLeft > 0.5)) - 1
      fR = 2 * as.numeric((pRight > 0.5)) - 1
    }
  }
  return(c(cut, fL, fR))
}



###### part b ##########
set.seed(3)
n = 300
x = runif(n)
y = (rbinom(n, 1, (sin (4 * pi * x) + 1) / 2) - 0.5) * 2

T = 100
cut = rep(0, T)
classifier = matrix(0, n, T)
epsilon = rep(0, T)
alpha = rep(0, T)
Z = rep(0, T)
w = rep(1/n, n)

for (t in 1: T){
    f_out = f(x, y, w)
    cut[t] = f_out[1]
    classifier[, t] = sapply(x, function(z) ifelse(z <= cut[t], f_out[2], f_out[3]))
    epsilon[t] = sum(((y != classifier[, t]) * 1) * w)
    alpha[t] = 1/2 * log((1- epsilon[t]) / epsilon[t])
    Z[t] = 2 * sqrt(epsilon[t] * (1 - epsilon[t]))
    w =  w / Z[t] * exp( -alpha[t] * y * classifier[, t])
}  

F = matrix(0, n, T)
for (t in 1: T){
  F[, t] = alpha[t] * classifier[, t]
}

output = sign(apply(F, 1, sum))
table(output, y)
sum(output !=y ) / n # 0.1733333


# generate another indep dataset, run the previous code again withou set.seed
rm(list = ls())

f = function(x, y, w){
  score = -Inf
  cut = -Inf
  fL = -Inf
  fR = -Inf
  for (c in x) {
    left = (x <= c)
    right = (x > c)
    pLeft = 0
    pRight = 0
    if (sum(left) > 0) {
      pLeft = sum(w[left][y[left] == 1]) / sum(w[left]) 
    }
    if (sum(right) > 0) {
      pRight = sum(w[right][y[right] == 1]) / sum(w[right]) 
    }
    giniLeft = pLeft * (1 - pLeft)
    giniRight = pRight * (1 - pRight)
    tmp = ((-giniLeft * sum(w[left])) + (-giniRight * sum(w[right]))) / sum(w)
    if (tmp > score) {
      score = tmp
      cut = c
      fL = 2 * as.numeric((pLeft > 0.5)) - 1
      fR = 2 * as.numeric((pRight > 0.5)) - 1
    }
  }
  return(c(cut, fL, fR))
}


n = 300
x = runif(n)
y = (rbinom(n, 1, (sin (4 * pi * x) + 1) / 2) - 0.5) * 2

T = 100
cut = rep(0, T)
classifier = matrix(0, n, T)
epsilon = rep(0, T)
alpha = rep(0, T)
Z = rep(0, T)
w = rep(1/n, n)

for (t in 1: T){
  f_out = f(x, y, w)
  cut[t] = f_out[1]
  classifier[, t] = sapply(x, function(z) ifelse(z <= cut[t], f_out[2], f_out[3]))
  epsilon[t] = sum(((y != classifier[, t]) * 1) * w)
  alpha[t] = 1/2 * log((1- epsilon[t]) / epsilon[t])
  Z[t] = 2 * sqrt(epsilon[t] * (1 - epsilon[t]))
  w =  w / Z[t] * exp( -alpha[t] * y * classifier[, t])
}  

F = matrix(0, n, T)
for (t in 1: T){
  F[, t] = alpha[t] * classifier[, t]
}

output = sign(apply(F, 1, sum))
sum(output !=y ) / n # 0.16,  0.1133333, 0.1666667, 0.1766667


