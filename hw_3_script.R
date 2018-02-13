rm(list = ls())
# install.packages("quadprog")
library(quadprog)

################################ quesiton 1 (a) #####################################
# generate a set of separable data
set.seed(1)
n = 40
p = 2
xpos = matrix(rnorm(n * p, mean = 0, sd = 1), n, p)
xneg = matrix(rnorm(n * p, mean = 4, sd = 1), n, p)
x = rbind(xpos, xneg)
y = matrix(c(rep(1, n), rep(-1, n)))

# primal form of linear separable SVM
# solve QP with quadprog and the perturbance hack
# problems of the form min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0.
# quadprog solver requires that the D matrix be symmetric positive definite
# we can perturb D by a small diagonal matrix and obtain positive definite matrix
# choose eps relatively small value for the diagonal perturbance
eps = 5e-4
Q = sapply(1: dim(x)[1], function(i) y[i] * t(x)[, i])


# beta^T * beta , where beta = (beta0, beta1, beta2)
D1 = rbind( c(eps, rep(0, p)), cbind(rep(0, p), diag(p)))

d1 = rep(0, p + 1)
A1 = t(cbind(y, t(Q)))
b01 = rep(1, 2*n)

results = solve.QP(D1, d1, A1, b01, meq = 0, factorized = FALSE)

# visualize the data
library(data.table)
library(ggplot2)
df = cbind(x, y)
df = as.data.frame(df)

plotMargin = function(w = 1*c(-1, 1), b = 1){
  x1 = seq(-20, 20, by = .01)
  x2 = (-w[1]*x1 + b)/w[2]
  l1 = (-w[1]*x1 + b + 1)/w[2]
  l2 = (-w[1]*x1 + b - 1)/w[2]
  dt = data.table(X1 = x1, X2 = x2, L1 = l1, L2 = l2)
  ggplot(dt) + geom_line(aes(x = X1, y = X2)) + geom_line(aes(x = X1, y = L1), color = "blue") + geom_line(aes(x = X1, y = L2), color = "green")+
    geom_hline(yintercept = 0, color = "red") + geom_vline(xintercept = 0, color = "red") + xlim(-5, 5) + ylim(-5, 5)+
    labs(title = paste0("w = (", w[1], ",", w[2], "), b = ", b))
}

plotMargin(w = c( -results$solution[2], -results$solution[3]), b = results$solution[1]) + 
  geom_point(data = df, aes(x = V1, y = V2, color = as.factor(V3))) + labs(title = "Using primal via solve.QP")

# intercept, and slope
primalline = c(- results$solution[1] / results$solution[3],  - results$solution[2]/results$solution[3])



################################ question 1 (b) ######################################

# dual form of linear separable SVM
# build the system matrices
D = t(Q) %*% Q

d = matrix(1, nrow = dim(x)[1])
b0 = rbind(matrix(0, 1, 1), matrix(0, 2*n, 1))
A = t(rbind(matrix(y, 1, 2*n), diag(nrow = 2*n)))
sol = solve.QP(D + eps * diag(2*n), d, A, b0, meq = 1, factorized = FALSE)
qpsol = matrix(sol$solution, nrow = 2*n)

findLine = function(a, y, X){
  nonzero = abs(a) > 1e-5
  W = rowSums(sapply(which(nonzero), function(i) a[i] * y[i] * X[i,]))
  b = mean(sapply(which(nonzero), function(i) X[i,] %*% W- y[i]))
  slope = -W[1] / W[2]
  intercept = b / W[2]
  return(c(intercept, slope))
}

findwb = function(a, y, X){
  nonzero = abs(a) > 1e-5
  W = rowSums(sapply(which(nonzero), function(i) a[i] * y[i] * X[i,]))
  b = mean(sapply(which(nonzero), function(i) X[i,] %*% W- y[i]))
  return(c(W, b))
}


qpline = findLine(qpsol, y, x)
dual_result = findwb(qpsol, y, x)

# plot of dual form via solve.QP
plotMargin(w = c(0.9332568, 0.3848944), b = -dual_result[3]) + geom_point(data = df, aes(x = V1, y = V2, color = as.factor(V3))) + 
  labs(title = "Using dual via solve.QP")


# use e1017 svmfit function 
library(e1071)
svm.fit = svm(y ~ ., data = data.frame(x, y), type ='C-classification', kernel='linear', scale = FALSE, cost = 10000)
w = t(svm.fit$coefs) %*% svm.fit$SV  # w = alpha_i*y_i*x_i
svmline = c(svm.fit$rho / w[1, 2], -w[1, 1] / w[1, 2]) 

# Visualize three lines together
plot(x,col=ifelse(y>0,"red", "blue"), pch = 19, cex = 1.2, lwd = 2, xlab = "X1", ylab = "X2", cex.lab = 1.5)
legend("bottomleft", c("Positive","Negative"),col=c("red", "blue"),pch=c(19, 19),text.col=c("red", "blue"), cex = 1.5)
abline(a = svmline[1], b = svmline[2], col = "black", lty = 3, lwd = 2)
abline(a = qpline[1], b = qpline[2], col = "red", lty = 3, lwd = 2)
abline(a = primalline[1], b = primalline[2], col = "yellow", lty = 3, lwd = 2)


# e1071 line
plotMargin(w = c(0.933874, 0.3844957), b = - svm.fit$rho) + geom_point(data = df, aes(x = V1, y = V2, color = as.factor(V3))) + 
labs(title = "Using e1071 package")

################################ question 1 (c) ######################################

# nonseparable SVM
# using solve.QP on dual form
rm(list = ls())
set.seed(2)
n = 40
p = 2

# Generate the positive and negative examples
xpos = matrix(rnorm(n * p, mean = 0, sd = 1), n, p)
xneg = matrix(rnorm(n * p, mean = 1.5, sd = 1), n, p)
x = rbind(xpos, xneg)
y = matrix(as.factor(c(rep(1, n), rep(-1, n))))
y1 = as.numeric(as.character(y))


eps = 5e-4
Q = sapply(1: dim(x)[1], function(i) y1[i] * t(x)[, i])

D = t(Q) %*% Q
C = 1
d = matrix(1, nrow = dim(x)[1])
b0 = rbind(matrix(0, 1, 1), matrix(0, 2*n, 1), matrix(-C, 2*n, 1))
A = t(rbind(matrix(y1, 1, 2*n), diag(nrow = 2*n), -diag(2*n)))
sol = solve.QP(D + eps * diag(2*n), d, A, b0, meq = 1, factorized = FALSE)
qpsol_dual = matrix(sol$solution, nrow = 2*n)

findw = function(a, y, X){
  nonzero = abs(a) > 1e-5
  W = rowSums(sapply(which(nonzero), function(i) a[i] * y[i] * X[i,]))
  return(W)
}

w_dual = findw(qpsol_dual, y1, x)
b_dual = median(y1 - c(w_dual[1], w_dual[2]) %*% t(x))

# using e1071 package
svm.fit = svm(y ~ ., data = data.frame(x, y), type = 'C-classification', kernel='linear', scale = FALSE, cost = 1)
w = t(svm.fit$coefs) %*% svm.fit$SV
b = -svm.fit$rho


# Visualize the data
plot(x,col=ifelse(y>0,"red", "blue"), pch = 19, cex = 1.2, lwd = 2, xlab = "X1", ylab = "X2", cex.lab = 1.5, main = "Decision Boundary via e1071 package and solve.QP")
legend("topright", c("Positive","Negative"),col=c("red", "blue"),pch=c(19, 19),text.col=c("red", "blue"), cex = 1.5)
abline(a= -b/w[1,2], b=-w[1,1]/w[1,2], col="black", lty=3, lwd = 2) # decision boundary from svm.fit
abline(a= -b_dual/w_dual[2], b=-w_dual[1]/w_dual[2], col="red", lty=3, lwd = 2) # decision boundary from dual form solve.QP

################################ question 1 (c) try primal (not required in HW) ######################################
# nonseperable SVM
# try using solve.QP on primal form

rm(list = ls())
set.seed(2)
n = 40
p = 2

# Generate the positive and negative examples
xpos = matrix(rnorm(n * p, mean = 0, sd = 1), n, p)
xneg = matrix(rnorm(n * p, mean = 1.5, sd = 1), n, p)
x = rbind(xpos, xneg)
y = matrix(as.factor(c(rep(1, n), rep(-1, n))))
y1 = as.numeric(as.character(y))


Q = sapply(1: dim(x)[1], function(i) y1[i] * t(x)[, i])
C = 1
d = c(0, rep(0, p), rep(-C, 2*n)) # -C * sum(zeta)
eps = 5e-4
D = diag(c(eps, rep(1, p), rep(eps, 2*n))) # beta^T * beta
In = diag(2*n) # 80*80 identity matrix

A_1 = cbind(matrix(0, ncol = p + 1, nrow = 2*n), In) # zeta_i > 0, for all i; 2n rows
A_2 = cbind(y1, t(Q), In) # zeta_i + beta_0 * y_i + beta^T * x_i * y_i >= 1, for all i; 2n rows
rownames(A_1) = NULL
rownames(A_2) = NULL
colnames(A_1) = NULL
colnames(A_2) = NULL

A = t(rbind(A_1, A_2))
b0 = c(rep(0, 2*n), rep(1, 2*n))

results_nonsep = solve.QP(D, d, A, b0)
boptim = results_nonsep$solution

beta0 = boptim[1]
beta = boptim[1 + (1: p)]
zeta = boptim[(p + 1) + (1: 2*n)]




################### question 1 (d) ########################

# for south african heart disease 
rm(list = ls())
set.seed(3)
library(quadprog)
library(ElemStatLearn)
str(SAheart)
dta = SAheart 
y = ifelse(dta$chd == 0, -1, 1)
xf1 = ifelse(dta$famhist == "Absent", 1, 0)
xf2 = ifelse(dta$famhist == "Present", 1, 0)
xf = cbind(xf1, xf2)
x = cbind(as.matrix(dta[ , - c(which(colnames(dta) == "famhist"), which(colnames(dta) == "chd"))]), xf)

p = ncol(x)
n = nrow(x)
mydta = cbind(x, y)

# set a series of cost values
C = c(exp(seq(-3, 0, 0.5)), seq(1.05, 2, 0.05))

# 10 fold cross validation to choosing the cost parameter c
nfold = 10
infold = sample(rep(1: nfold, length.out = n))

ny = dim(mydta)[2]
cv_error = matrix(NA, nrow = nfold, ncol = length(C))

for(i in 1: length(C)){
    c = C[i]
  for(l in 1: nfold){
    traindata = mydta[infold != l, ]
    trainx = traindata[, -ny]
    trainy = traindata[, ny]
    Q = sapply(1: dim(trainx)[1], function(i) trainy[i] * t(trainx)[, i])
    train_n = dim(trainx)[1]
    d = c(0, rep(0, p), rep(-c, train_n)) 
    eps = 5e-4
    D = diag(c(eps, rep(1, p), rep(eps, train_n))) 
    In = diag(train_n) 
    
    A_1 = cbind(matrix(0, ncol = p + 1, nrow = train_n), In) 
    A_2 = cbind(trainy, t(Q), In)
    rownames(A_1) = NULL
    rownames(A_2) = NULL
    colnames(A_1) = NULL
    colnames(A_2) = NULL
    A = t(rbind(A_1, A_2))
    b0 = c(rep(0, train_n), rep(1, train_n))
    results_nonsep <- solve.QP(D, d, A, b0)
    boptim = results_nonsep$solution
    
    beta0 = boptim[1]
    beta = boptim[1 + (1: p)]
    zeta = boptim[(p + 1) + (1: train_n)]
    
    testx = mydta[infold == l, -ny]
    testy = mydta[infold == l, ny]
    Y_pred = sign(apply(testx, 1, function(x) beta0 + sum(beta * as.vector(x))))
    cv_error[l, i] = sum(Y_pred != testy)
  }
}

which.min(apply(cv_error, 2, sum)) # give the smallest sum of predicting errors
C[4] # 0.2231302
# the sum of predicting errors for each C
# 132 131 131 129 131 130 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131 131
# the fourth gives the smallest sums of error
# thus we choose the cost = 0.2231302 based on 10 fold cross validation



# define a svm function
svm.function = function(x, y, c){
  Q = sapply(1: dim(x)[1], function(i) y[i] * t(x)[, i])
  n = dim(x)[1]
  p = dim(x)[2]
  d = c(0, rep(0, p), rep(-c, n)) 
  eps = 5e-4
  D = diag(c(eps, rep(1, p), rep(eps, n))) 
  In = diag(n) 
  
  A_1 = cbind(matrix(0, ncol = p + 1, nrow = n), In) 
  A_2 = cbind(y, t(Q), In)
  rownames(A_1) = NULL
  rownames(A_2) = NULL
  colnames(A_1) = NULL
  colnames(A_2) = NULL
  A = t(rbind(A_1, A_2))
  b0 = c(rep(0, n), rep(1, n))
  results = solve.QP(D, d, A, b0)
  boptim = results$solution
  
  beta0 = boptim[1]
  beta = boptim[1 + (1: p)]
  zeta = boptim[(p + 1) + (1: n)]
  resultlist = list(beta0, beta, zeta)
  return(resultlist)
}

# plug c  = 0.2231302 into my SVM code
mysvm_result = svm.function(x = x, y = y, c = 0.2231302)
mysvm_result[[1]]
mysvm_result[[2]]



# using e1071 package
library(e1071)
#y_factor = as.factor(y)
svm.fit_heart = svm(y ~ ., data = data.frame(x, y), type = 'C-classification', kernel='linear', scale = FALSE, cost = 0.2231302)
summary(svm.fit_heart)
pred = predict(svm.fit_heart, x)
table(pred, y_factor)

w = t(svm.fit_heart$coefs) %*% svm.fit_heart$SV
b = -svm.fit_heart$rho

# show results in table
library(xtable)
tbl_e1071 = xtable(t(w), digits = 9)



################## Question 3 ########################
rm(list = ls())
getwd()
setwd("/Users/rociozhong/Library/Mobile Documents/com~apple~CloudDocs/STAT_542")
shooting = read.csv("Mass Shootings Dataset Ver 2.csv", header = T, sep = ",", na.strings=c("","NA"))

# Any correlation with calendar dates? Do we have more deadly days, weeks or months on average
head(shooting)
dta = shooting[, c("Date", "Fatalities", "Injured", "Total.victims")]  #"Gender", "Race", "Gender", "Latitude", "Longitude")]

# dta[ dta == "Unknown"] <- NA
# dta[ dta == "unknown"] <- NA
# dta[ dta == "Unclear"] <- NA
# dta$Gender[dta$Gender == "M"] = "Male"
# dta$central = ifelse(dta$Longitude < -90 & dta$Longitude > -110, 1, 0)
# dta$male = ifelse(dta$Gender == "Male", 1, 0)
# dta$white = ifelse(grepl("White", dta$Race) == TRUE, 1, 0)


dta = na.omit(dta)
str(dta)


# add new variables
dta$time = as.Date(as.factor(dta$Date), format = "%m/%d/%Y")


library(lubridate)  # for "datetime"

dta$day= as.numeric(as.character(as.factor(day(dta$time))))
dta$year = as.numeric(as.character(as.factor(year(dta$time))))
dta$month = as.numeric(as.character(as.factor(month(dta$time))))

# data visualization
library(dplyr)
by_year = group_by(dta, year)
shot_year = summarise(by_year, shooting = n(), vic = sum(Total.victims), inj = sum(Injured), fat = sum(Fatalities))

library(reshape2)
year_melt = melt(shot_year, id = "year", measure = c("vic", "fat"))
ggplot(year_melt, aes(x = year, y = value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") # vic, inj, fat per year
ggplot(shot_year, aes(x = year, y = shooting)) + geom_bar(stat = "identity") + scale_x_continuous(breaks = seq(min(shot_year$year), max(shot_year$year), by = 3)) # shooting per year

# add a new variable weekday, weekday = 1, weekend = 0
# library(chron)
# dta$wdy = ifelse(is.weekend(dta$time) == FALSE, 1, 0)

# add a new variable that differentiates year before 2000 and year after 2000
dta$later2000 = ifelse(dta$year > 2000, 1, -1)
dta$later2000 = factor(dta$later2000)


# remove the outlier
dta = dta[-1, ]
plot(c(dta$Injured, dta$Fatalities), col=ifelse(as.numeric(as.character(dta$later2000)) > 0, "red", "blue"), 
     pch = 19, cex = 1.2, lwd = 2, xlab = "X1", ylab = "X2", cex.lab = 1.5)



x = cbind(dta$Injured, dta$Fatalities)
y = dta$later2000
mydta = data.frame(y = y, x = x)
colnames(mydta) = c("y", "x1", "x2")

svm.fit = svm(y ~., data = mydta, gamma=10, cost=10)


pred = predict(svm.fit, mydta[, -1])
table(pred, mydta[, 1])

# Compute accuracy
sum(pred==y)/length(y)






