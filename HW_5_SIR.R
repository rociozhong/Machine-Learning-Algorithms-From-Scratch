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


#################### question 3 ###########################
## predict the movie, "Star Wars: the Last Jedi" box office and whether its vote average is greater than 7.  
## the dataset is from https://www.kaggle.com/tmdb/tmdb-movie-metadata
                     
rm(list = ls())
getwd()
setwd("/Users/rociozhong/Library/Mobile Documents/com~apple~CloudDocs/STAT_542")
load_data = read.csv("tmdb_5000_movies.csv", header = T, sep = ",")
load_data = na.omit(load_data)


#### genre ####
temp = load_data$genres
library(stringr)
after_remove = str_replace_all(temp, "[[:punct:]]", " ")
after_remove = gsub("\\s+", " ", str_trim(after_remove))

stopwords = c("id", "name")
library(tm)
after_remove = removeWords(after_remove, stopwords)
after_remove = gsub('[0-9]+', '', after_remove)
after_remove = gsub("\\s+", " ", str_trim(after_remove))
load_data$temp_genre = after_remove


library(splitstackshape)
library(reshape2)
mydata = concat.split.expanded(load_data, split.col = "temp_genre", sep = " ", type = "character",
                      mode = "binary", fixed = TRUE, fill = 0)

#### production companies ####

company = load_data$production_companies
company = as.character(company)
company_upd = purrr::map(company, jsonlite::fromJSON) # or use: lapply(company, jsonlite::fromJSON)

# load company inf highest revenue production company information, grab first 50 companies' name
com_inf = read.csv("movie_company.csv", header = T, sep = ",")
com_inf = com_inf[, -c(5:7)]
com_inf$key = gsub("([A-Za-z]+).*", "\\1", com_inf$Company)

# find the movie id for those 50 production companies
list_com = lapply(company_upd, function(x) x$name) 

index = list()
movie_id = list()
for(i in 1:10){
  index[[i]] = grep(com_inf$key[i], list_com)
  movie_id[[i]] = load_data$id[index[[i]]]
}
# movies id produced by first 10 companies
id_10 = unique(unlist(movie_id)) 

# if movie is produced by first big 10, set as 1, otherwise 0.
mydata$large_prod = ifelse(mydata$id %in% id_10, 1, 0)


# add time new variables
mydata$time = as.Date(as.factor(mydata$release_date), format = "%Y-%m-%d")


library(lubridate)  # for "datetime"

mydata$day= as.numeric(as.character(as.factor(day(mydata$time))))
mydata$year = as.numeric(as.character(as.factor(year(mydata$time))))
mydata$month = as.numeric(as.character(as.factor(month(mydata$time))))

mydata = na.omit(mydata)
# create three indicators identifying Thanksgiving, or the Christmas season (dec 18 - 31) or summer season(Memorial - labor))

# during thanks giving break

library(chron)
thanks = holiday(unique(mydata$year), Holiday = "USThanksgivingDay")
thanks_year = as.numeric(as.character(as.factor(year(thanks))))
thanks_day = as.numeric(as.character(as.factor(day(thanks))))
range(thanks_day) # 22 - 28

# thanks giving season break  --- 1
mydata$thanksg = NA
mydata$thanksg[which(mydata$month == "11")] = 
  ifelse(mydata$day[which(mydata$month == "11")] <= 28 & mydata$day[which(mydata$month == "11")] >= 22, 1, 0)

# chrismas season dec 18-31 --- 2
mydata$chrism = NA
mydata$chrism[which(mydata$month == "12")] = 
  ifelse(mydata$day[which(mydata$month == "12")] <= 31 & mydata$day[which(mydata$month == "12")] >= 18, 2, 0)


# summer season 5 - 8 --- 3
mydata$summer = NA
mydata$summer = ifelse(mydata$month %in% c(5,6,7,8), 3, 0)


# rest date --- 0
mydata[, c("chrism", "thanksg")][is.na(mydata[, c("chrism", "thanksg")])] = 0

mydata$holiday = rowSums(mydata[, c("chrism", "thanksg", "summer")])
# mydata$holiday = as.factor(mydata$holiday)

#### production countries ####

countries = mydata$production_countries
countries = as.character(countries)
countries_upd = lapply(countries, jsonlite::fromJSON)
countries_name = lapply(countries_upd, function(x) x$name)


for(i in 1: length(countries_name)){
  mydata$us = ifelse("United States of America" %in% countries_name[[i]], 1, 0)
}

mydata$lang = ifelse(mydata$original_language == "en", 1, 0)

# mydata$lang = as.factor(mydata$lang)
# mydata$us = as.factor(mydata$us)

names(mydata)


cols = c("genres", "homepage", "keywords", "original_language", "original_title", "overview",
         "popularity", "production_companies", "production_countries", "release_date", "spoken_languages", 
         "title", "vote_average", "vote_count", "tagline", "time", "temp_genre", 
         "thanksg", "summer", "chrism", "temp_genre_Science") # since Science row == Fiction

mydata$budget = log(mydata$budget )
mydata$revenue = log(mydata$revenue )


mydata = mydata[-which(is.infinite(mydata$budget) == TRUE), ]
mydata = mydata[-which(is.infinite(mydata$revenue) == TRUE), ]


final_data = mydata[, - which(names(mydata) %in% cols)]


# odd id as train data
train_data = final_data[which(final_data$id %% 2 != 0), ]

# even id as test data
test_data  = final_data[which(final_data$id %% 2 == 0), ]


set.seed(3)
# variable selection, linear regression
# variables chosed by AIC
var_sel = step(lm(revenue ~., data =  train_data[, -2], direction = "both", trace = 1))
summary(var_sel) # 15 vars


# important vars: budget, runtime, genre: action, adventure, comedy, crime, drama, history,
# thriller, western, large_prod, holiday

col_imp = c("revenue", "budget", "runtime", "large_prod", "holiday",
            "temp_genre_Adventure", "temp_genre_Action", "temp_genre_Animation", 
            "temp_genre_Documentary", "temp_genre_Drama", "temp_genre_Foreign", "temp_genre_Horror",
            "temp_genre_Romance", "temp_genre_War", "temp_genre_Western")

fit_final = lm(revenue ~., data = train_data[, which(names(train_data) %in% col_imp)])
summary(fit_final)

library(ggplot2)
p = ggplot(aes(x = actual, y = pred), 
           data = data.frame(actual = test_data$revenue, 
                             pred = predict(fit_final, test_data[, which(names(test_data) %in% col_imp)])))
p + geom_point() + geom_abline(color = "red")

# model prediction MSE
sqrt(mean((predict(fit_final, newdata = test_data[, which(names(test_data) %in% col_imp)]) - 
             test_data$revenue)^2))

# revenue prediction
sqrt(mean((exp(predict(fit_final, newdata = test_data[, which(names(test_data) %in% col_imp)])) - 
  exp(test_data$revenue))^2))

# randomforest predicting vote_average
# predict classfication problem whether vote_average >7

cols_2 = c("genres", "homepage", "keywords", "original_language", "original_title", "overview",
           "popularity", "production_companies", "production_countries", "release_date", "spoken_languages", 
           "title", "revenue", "vote_count", "tagline", "time", "temp_genre", "revenue",
           "thanksg", "summer", "chrism", "temp_genre_Science", "us") # since Science row == Fiction


final_data_2 = mydata[, - which(names(mydata) %in% cols_2)]

# vote average > 7 is one class, otherwise another class
final_data_2$vote_average = ifelse(final_data_2$vote_average > 7, 1, 0)

# odd id as train data
train_data_2 = final_data_2[which(final_data_2$id %% 2 != 0), ]

# even id as test data
test_data_2  = final_data_2[which(final_data_2$id %% 2 == 0), ]



set.seed(333)
### 5 fold cross-validation on mtry
library(class)
library(randomForest)
library(miscTools)

nfold = 5
infold = sample(rep(1:nfold, length.out = dim(train_data_2)[1]))

K  = 15 #dim(train_data_2)[2] - 2 # maximum number of variables 
cverror = matrix(0, nrow = nfold, ncol = K)


for (l in 1:nfold)
{
  for (k in 1:K)
  {rf_temp = randomForest(as.factor(vote_average) ~., mtry = k, data = train_data_2[infold != l, -2], 
                             importance = TRUE)
  
  cv_predict = predict(rf_temp, train_data_2[infold == l, -2])
  
  cverror[l, k] = sum(cv_predict != train_data_2[infold == l, ]$vote_average)
  }
}

which.min(apply(cverror, 2, sum))  
# then mtry will be 3


rf = randomForest(as.factor(vote_average) ~., mtry = 3, data = train_data_2[, -2] , importance = TRUE)
plot(rf)
varImpPlot(rf, sort = T, main = "Variable Importance")

# variable importance table
var.imp = data.frame(importance(rf, type = 2))

var.imp$Variables = row.names(var.imp)
var.imp[order(var.imp$MeanDecreaseGini,decreasing = T),]

test_data_2$prediction_rf = predict(rf, test_data_2[, -2])

table(test_data_2$prediction_rf, test_data_2$vote_average)
(19 + 267) / dim(test_data_2)[1] #  error 0.1807838
(1259 + 27) / dim(test_data_2)[1] # accuracy  0.8128951


########## star war movie ###############
star_wars_2 = data.frame(
  budget = log(300000000),
  id = "NA",
  runtime = 152,
  status = "Released",
  vote_average = "NA",
  temp_genre_Action = 1,
  temp_genre_Adventure = 1,
  temp_genre_Animation = 0,
  temp_genre_Comedy = 0,
  temp_genre_Crime = 0,
  temp_genre_Documentary = 0,
  temp_genre_Drama = 0,
  temp_genre_Family = 0,
  temp_genre_Fantasy = 1,
  temp_genre_Fiction = 0,
  temp_genre_Foreign = 0,
  temp_genre_History = 0,
  temp_genre_Horror = 0,
  temp_genre_Movie = 0,
  temp_genre_Music = 0,
  temp_genre_Mystery = 0,
  temp_genre_Romance = 0,
  temp_genre_Thriller = 0,
  temp_genre_TV = 0,
  temp_genre_War = 0,
  temp_genre_Western = 0,
  large_prod = 1,
  day = 15,
  year = 2017,
  month = 12,
  holiday = 0,
  lang = 1)

newtrain2 = rbind(train_data_2, star_wars_2)
star_wars_predict = predict(rf, newtrain2[nrow(newtrain2), -2]) # vote average > 7



names(train_data[, which(names(train_data) %in% col_imp)])
star_wars_1 = data.frame(
  budget = log(300000000),
  revenue = "NA",
  runtime = 152,
  temp_genre_Action = 1,
  temp_genre_Adventure = 1,
  temp_genre_Animation = 0,
  temp_genre_Documentary = 0,
  temp_genre_Drama = 0,
  temp_genre_Foreign = 0,
  temp_genre_Horror = 0,
  temp_genre_Romance = 0,
  temp_genre_War = 0,
  temp_genre_Western = 0,
  large_prod = 1,
  holiday = 0)

star_rev_predict = predict(fit_final, newdata = star_wars_1)
exp(star_rev_predict)







