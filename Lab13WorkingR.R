#########################################################################################################
# Lab 13 
library(tidyverse)
library(e1071)
finches_data <- read_csv("zebrafinches.csv")

# Question 1 Part A:
skew <- skewness(finches_data$further)
n <- 25
t <- t.test(finches_data$further, mu = 0, 
                    conf.level = 0.95, alternative = "two.sided")$statistic 
fz <- dnorm(t)
Fz <- pnorm(t)

potential_error <- Fz+(skew/sqrt(n))*(((2*t^2)+1)/6)*fz


# Question 1 Part B:




# Question 1 Part C:




# Question 2 Part A:





# Question 2 Part B:






# Question 2 Part C:





# Question 2 Part D:








# Question 3 Part A:





# Question 3 Part B:






# Question 3 Part C:






#########################################################################################################