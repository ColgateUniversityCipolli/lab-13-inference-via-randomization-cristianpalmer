#########################################################################################################
# Lab 13 
library(tidyverse)
library(e1071)
library(boot)
finches_data <- read_csv("zebrafinches.csv")

# Question 1 Part A:
skew <- skewness(finches_data$further)
n <- 25
t <- t.test(finches_data$further, mu = 0, 
                    conf.level = 0.95, alternative = "less")$statistic 
fz <- dnorm(t)
Fz <- pnorm(t)

(potential_error <- Fz+(skew/sqrt(n))*(((2*t^2)+1)/6)*fz)


# Question 1 Part B:
skew <- skewness(finches_data$further)
n <- 25
values <- seq(-10, 10, by = 0.001)
t_values <- tibble(t_stat = numeric(length(values)))
for (i in 1:length(values)) {
  fz <- dnorm(values[i])
  Fz <- pnorm(values[i])
  potential_error <- (skew/sqrt(n))*(((2*values[i]^2)+1)/6)*fz
  t_values$t_stat[i] <- potential_error
}

ggplot(t_values, aes(x = values, y = t_stat)) +
  geom_line(color = "blue") +  # Line plot
  labs(title = "Potential Error vs. t", x = "t", y = "Potential Error") +
  theme_minimal()


# Question 1 Part C:
skew <- skewness(finches_data$further)
a <- 0.05

t <- qnorm(0.05)
fz <- dnorm(t)

(n <- ((skew/(6*0.10*a))*((2*t^2)+1)*fz)^2)

###########################################################################

# Question 2 Part A:

# resamples.null.closer

n <- 25
R <- 10000
resamples.null.closer <- tibble(t = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finches_data$closer,
                          size = nrow(finches_data),
                          replace = T)
  
  resamples.null.closer$t[i] <- (mean(curr.resample))/(sd(finches_data$closer)/sqrt(n)) 
}
# Center 
resamples.null.closer <- (resamples.null.closer)$t
resamples.null.closer <- resamples.null.closer - mean(resamples.null.closer)

# resamples.null.further

n <- 25
R <- 10000
resamples.null.further <- tibble(t = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finches_data$further,
                          size = nrow(finches_data),
                          replace = T)
  
  resamples.null.further$t[i] <- (mean(curr.resample))/(sd(finches_data$further)/sqrt(n)) 
}
# Center 
resamples.null.further <- (resamples.null.further)$t
resamples.null.further <- resamples.null.further - mean(resamples.null.further)


# resamples.null.diff

n <- 25
R <- 10000
resamples.null.diff <- tibble(t = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finches_data$diff,
                          size = nrow(finches_data),
                          replace = T)
  
  resamples.null.diff$t[i] <- (mean(curr.resample))/(sd(finches_data$diff)/sqrt(n)) 
}
# Center 
resamples.null.diff <- (resamples.null.diff)$t
resamples.null.diff <- resamples.null.diff - mean(resamples.null.diff)


# Question 2 Part B:






# Question 2 Part C:





# Question 2 Part D:



# Question 3 Part A:





# Question 3 Part B:






# Question 3 Part C:






#########################################################################################################