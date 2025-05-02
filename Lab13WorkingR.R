#########################################################################################################
# Lab 13 
library(tidyverse)
library(e1071)
library(boot)
library(boot.pval)

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
  labs(title = "Potential Error vs. t from -10 to 10", x = "t", y = "Potential Error") +
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
resamples.null.closer.shifted <- resamples.null.closer - mean(resamples.null.closer)

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
resamples.null.further.shifted <- resamples.null.further - mean(resamples.null.further)


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
resamples.null.diff.shifted <- resamples.null.diff - mean(resamples.null.diff)


# Question 2 Part B:

# Closer
t_closer <- t.test(finches_data$closer, mu = 0, 
                   conf.level = 0.95, alternative = "two.sided")$statistic
(p_closer <- mean(resamples.null.closer.shifted >= t_closer))

# Further
t_further <- t.test(finches_data$further, mu = 0, 
                    conf.level = 0.95, alternative = "two.sided")$statistic
(p_further <- mean(resamples.null.further.shifted <= t_further))

# Diff
t_diff <- t.test(finches_data$diff, mu = 0, 
                       conf.level = 0.95, alternative = "two.sided")$statistic
low <- -(t_diff) 
high <- (t_diff)
p.low = mean(resamples.null.diff.shifted <= low)
p.high = mean(resamples.null.diff.shifted >= high)
(p_diff = p.low + p.high)

# Question 2 Part C:

# Closer
(firth_percentile_closer <- quantile(resamples.null.closer.shifted, .05))
(firth_percentile_closer_t <- qt(0.05, df=length(finches_data$closer-1)))

# Further
(firth_percentile_further <- quantile(resamples.null.further.shifted, .05))
(firth_percentile_further_t <- qt(0.05, df=length(finches_data$further-1)))

# Diff
(firth_percentile_diff <- quantile(resamples.null.diff.shifted, .05))
(firth_percentile_diff_t <- qt(0.05, df=length(finches_data$diff-1)))



# Question 2 Part D:

# Boot strap confidence intervals using resampling

# Closer

library(boot)
boot.t <- function(d, i){
  d[i]  
}
boots.closer <- boot(data = finches_data$closer,
              statistic = boot.t,
              R = R) 

boot.ci(boots.closer, type="bca")

# Further
boots.further <- boot(data = finches_data$further,
                     statistic = boot.t,
                     R = R) 

boot.ci(boots.further, type="bca")

# Diff
boots.diff <- boot(data = finches_data$diff,
                      statistic = boot.t,
                      R = R) 

boot.ci(boots.diff, type="bca")


# Question 3 Part A:

# Randomization Test Closer

rand.closer <- tibble(t = rep(NA, R))

for(i in 1:R){
  curr.rand <- resamples.null.closer.shifted *
    sample(x = c(-1, 1),
           size = length(resamples.null.closer.shifted),
           replace = T)
  
  rand.closer$t[i] <- mean(curr.rand)
}
# Thinking is hard
rand.closer.shifted <- rand.closer |>
  mutate(t = t + mean(resamples.null.closer)) # shifting back

# Randomization Test Further

rand.further <- tibble(t = rep(NA, R))

for(i in 1:R){
  curr.rand <- resamples.null.further.shifted *
    sample(x = c(-1, 1),
           size = length(resamples.null.further.shifted),
           replace = T)
  
  rand.further$t[i] <- mean(curr.rand)
}
# Thinking is hard
rand.further.shifted <- rand.further |>
  mutate(t = t + mean(resamples.null.further)) # shifting back

# Randomization Test Diff

rand.diff <- tibble(t = rep(NA, R))

for(i in 1:R){
  curr.rand <- resamples.null.diff.shifted *
    sample(x = c(-1, 1),
           size = length(resamples.null.diff.shifted),
           replace = T)
  
  rand.diff$t[i] <- mean(curr.rand)
}
# Thinking is hard
rand.diff.shifted <- rand.diff |>
  mutate(t = t + mean(resamples.null.diff)) # shifting back


# Question 3 Part B:

# closer randomization p-value
(p_closer_randomization <- mean(rand.closer >= t_closer))

# further randomization p-value
(p_further_randomization <- mean(rand.further <= t_further))

# diff randomization p-value
low <- -(t_diff) 
high <- (t_diff)
p.low.rand = mean(rand.diff <= low)
p.high.rand = mean(rand.diff >= high)
(p_diff_rand = p.low.rand + p.high.rand)


# Question 3 Part C:

