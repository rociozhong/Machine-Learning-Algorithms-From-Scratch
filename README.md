# Machine Learning Algorithms From Scratch

## Lasso

## Support Vector Machine
I install the "quadprog" package and utilize the function "solve.QP" to solve SVM. The "solve.QP" function is trying to perform the minimizing problem: 

![](http://latex.codecogs.com/gif.latex?%5Ctext%7Bminimize%7D%5Cquad%20%5Cfrac%7B1%7D%7B2%7D%20%5Cbeta%20%5ET%20D%20%5Cbeta%20-d%5ET%20%5Cbeta%20%5Cquad%20%5Ctext%7Bsubject%20to%7D%20%5Cquad%20A%5ET%5Cbeta%20%5Cgeq%20b_0)

I formulate both the primal and dual form of SVM. And compare with "e1071" package

## Nadaraya-Watson kernel regression estimator for a one dimensional problem
I consider just the Gaussian kernel function, and include a tuning parameter for the bandwidth. Then I download the “Video Game Sales with Ratings” dataset from Kaggle https:// www.kaggle.com/rush4ratio/video-game-sales-with-ratings. Perform my kernel estimator to estimate the log(1 + Global Sales) using each of the three scores: Critic Score, Critic Count and User Score. For this case, I ignore observations with any missing value (ignore text value too).  And I tune the bandwidth using cross-validation. 

## A one-dimensional adaboost using a stump model as the weak learner

Algorithm: A stump model, with weights
Input: A set of data 
1. Search for a splitting rule that will maximize the weighted reduction of Gini impurity.
2. Calculate the left and the right node weighted predictions respectively.
Output: The cutting point c, and the mode predictions.

## Sliced inverse regression
I use 10 as the number of slices,validate it by comparing to the "dr" package.
