
#################### question 1 #########################
rm(list = ls())
require(graphics)
summary(faithful)
attach(faithful)
plot(faithful[,1] ~ faithful[,2], pch = 22, bg = "gray", xlab = "Waiting time [min]", ylab = "Eruption duration [min]")
lines(ksmooth(faithful[,2], faithful[,1], kernel = "normal", bandwidth = 9), col = "red", lwd = 2)

#result_ks = ksmooth(faithful[,2], faithful[,1], kernel = "normal", bandwidth = 9)

n = nrow(faithful)
x = faithful[,2]
y = faithful[,1]
# Nadaraya-Watson Kernel regression estimator for a one dimensional problem

mfit = function(z,x,y,h){
  u = (x - z) / h
  ks = dnorm(u)
  w = ks / sum(ks)
  return(m = sum(w * y))
}


temp = unlist(lapply(x, function(a) mfit(a, x = x, y = y, h = 6)))

plot(temp ~ x, pch = 22, bg = "gray", xlab = "Waiting time [min]", ylab = "Eruption duration [min]")
lines(ksmooth(faithful[,2], faithful[,1], kernel = "normal", bandwidth = 6), col = "red", lwd = 2)

# part b
rm(list = ls())
getwd()
setwd("/Users/rociozhong/Library/Mobile Documents/com~apple~CloudDocs/STAT_542")
data = read.csv("Video_Games_Sales_as_at_22_Dec_2016.csv", header = T, sep = ",", na.strings=c("","NA"))
names(data)


mydta = data[, c("Global_Sales", "Critic_Score", "Critic_Count", "User_Score")]
str(mydta)
mydta[ mydta == "tbd"] = NA 
mydta = na.omit(mydta)
mydta$User_Score = as.numeric(as.character(mydta$User_Score))
mydta$y = log( mydta$Global_Sales + 1) 

# cross-validation selection of bandwidth

CV = function(x, y, h, weight = false){#Gaussian kernel
  tmp = order(x)
  x = x[tmp]
  y = y[tmp]
  k = function(u){dnorm(u)} #Gaussian kernel
  mfit = function(z, x, y, k, h){
    u = (x-z) / h
    ks = k(u)
    w = ks / sum(ks)
    return(list(w = w, m = sum(w*y)))
  }
  
  if(!weight)wt = 1
  else wt = (abs(x-0.5) <= 0.4)
  n = length(y)
  mhi = rep(0, n)
  mh = rep(0, n)
  S = matrix(0, n, n)
  for(i in 1:n){
    mhi[i] = mfit(x[i], x[-i], y[-i], k, h)$m 
    
    a = mfit(x[i], x, y, k, h)
    mh[i] = a$m
    S[i,] = a$w
  }
  
  cv = mean((y - mhi)^2 * wt)
  gcv = n*sum((y - mh)^2) / (n - sum(diag(S)))^2
  
  return(list(cv = cv, gcv = gcv))
}


bw = seq(0.1, 10, 0.1)
result_cv = rep(0, length(bw))

video_function = function(x){
for( i in seq_along(bw)){
  result_cv[i] = CV(x = x, y = mydta$y, h = bw[i], weight = FALSE)$cv
}
return(c(bw[which.min(result_cv)], min(result_cv, na.rm = T)))
}


video_function(x = mydta$Critic_Score) # bandwidth = 1.2, cv = 0.1608793, smallest cv

video_function(x = mydta$Critic_Count) # bandwidth = 4.200000, cv = 0.164162

video_function(x = mydta$User_Score) # bandwidth = 0.2, cv = 0.2009796

##################### question 2 part a ######################
rm(list = ls())
library(randomForest)
N = 200
P = 4
set.seed(2)
X = matrix(rnorm(N * P, mean = 0, sd = 1), N, P) 

fittedy = matrix(0, nrow = N, ncol = 20) 
truey = matrix(0, nrow = N, ncol = 20) 

rf.fun = function(p, q){
for (i in 1: 20){
  truey1 = 1 + 0.5 * apply(X, 1, sum) + rnorm(N)
  rf.fit = randomForest(X, truey1, mtry = p, nodesize = q, xtest = X, ytest = truey1)
  fittedy[, i] = rf.fit$predicted
  truey[, i] = as.vector(truey1)
}
  return (list(fittedy = fittedy, truey = truey))
}

sample_cov = function(x, y){
  n_samples = nrow(X)
  scale_factor = (n_samples - 1) / n_samples
  scale_factor * cov(x, y)
}

# mtry = 1, nodesize = 20
sum(mapply(sample_cov, as.data.frame(rf.fun(p = 1, q = 20)$truey), as.data.frame(rf.fun(p = 1, q = 20)$fittedy)))
# mtry = 1, nodesize = 30
sum(mapply(sample_cov, as.data.frame(rf.fun(p = 1, q = 30)$truey), as.data.frame(rf.fun(p = 1, q = 30)$fittedy)))
# mtry = 1, nodesize = 40
sum(mapply(sample_cov, as.data.frame(rf.fun(p = 1, q = 40)$truey), as.data.frame(rf.fun(p = 1, q = 40)$fittedy)))



# mtry = 3, nodesize = 20
sum(mapply(sample_cov, as.data.frame(rf.fun(p = 3, q = 20)$truey), as.data.frame(rf.fun(p = 3, q = 20)$fittedy)))
# mtry = 3, nodesize = 30
sum(mapply(sample_cov, as.data.frame(rf.fun(p = 3, q = 30)$truey), as.data.frame(rf.fun(p = 3, q = 30)$fittedy)))
# mtry = 3, nodesize = 40
sum(mapply(sample_cov, as.data.frame(rf.fun(p = 3, q = 40)$truey), as.data.frame(rf.fun(p = 3, q = 40)$fittedy)))

# increase nodesize, decreases the df; increase the mtry, increases the df.
############ part b ###########
rm(list = ls())
N = 200
P = 4
set.seed(2)
X = matrix(rnorm(N * P, mean = 0, sd = 1), N, P) 

fittedy = matrix(0, nrow = N, ncol = 20) 
truey = matrix(0, nrow = N, ncol = 20) 

ntree.fun = function(z){
  for (i in 1: 20){
    truey1 = 1 + 0.5 * apply(X, 1, sum) + rnorm(N)
    rf.fit = randomForest(X, truey1, ntree = z, xtest = X, ytest = truey1)
    fittedy[, i] = rf.fit$predicted
    truey[, i] = as.vector(truey1)
  }
  return (list(fittedy = fittedy, truey = truey))
}


mean((ntree.fun(500)$fittedy - mean(ntree.fun(500)$fittedy))^2)
mean((ntree.fun(800)$fittedy - mean(ntree.fun(800)$fittedy))^2)
mean((ntree.fun(1000)$fittedy - mean(ntree.fun(1000)$fittedy))^2)

################  question 3 part a ###################
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


