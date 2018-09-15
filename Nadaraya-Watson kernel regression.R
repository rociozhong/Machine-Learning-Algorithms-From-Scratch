## Write my own code to fit a Nadaraya-Watson kernel regression estimator for a one dimensional problem.
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
