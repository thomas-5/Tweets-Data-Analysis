set.seed(20880634)
library(dplyr)
library(rootSolve)
library(moments)
library(ggplot2)
library(MASS)

# Read in dataset
dataset <- read.csv('./dataset.csv')

a1 <- dataset$urls.binary

thetahat <- mean(dataset$urls.binary)

# The number of samples
n <- length(a1)

# The observed data
y <- thetahat * n

BinomRLF <- function(theta, n, thetahat, y){
  (theta/thetahat)^y * ((1-theta)/(1-thetahat))^(n-y)
}

theta <- seq(0.34, 0.46, 0.0005)

plot(theta, BinomRLF(theta, n, thetahat, y), xlab = expression(theta),
     ylab = expression(paste("R(", theta, ")")), type = "l", lwd = "2",
     main = "Binomial relative likelihood function", las = 1,
     col = "navy")
abline(h = 0.15, col = "red", lwd = 2)

# The 15% likelihood interval for the data
Interval_15 <- uniroot.all(function(x) BinomRLF(x, n, thetahat, y) - 0.15, lower = 0.34, upper = 0.46)

#Observed value of likelihood test statistics
lambda <- -2*log(((0.19^y)*(1-0.19)^(n-y))/((thetahat^y)*(1-thetahat)^(n-y)))

a1_p <- 2*(1-pnorm(sqrt(lambda)))

tod_hour_O <- dataset$time.of.day[dataset$username == "@ONThealth"]/3600
tod_hour_A <- dataset$time.of.day[dataset$username == "@GoAHealth"]/3600

nA <- length(tod_hour_A)
nO <- length(tod_hour_O)

# qq-plot for analysis 2
qqnorm(tod_hour_O, pch = 1, frame = FALSE, main = "Normal Q-Q Plot for @ONThealth")
qqline(tod_hour_O, col = "steelblue", lwd = 2)

qqnorm(tod_hour_A, pch = 1, frame = FALSE, main = "Normal Q-Q Plot for @GoAHealth")
qqline(tod_hour_A, col = "steelblue", lwd = 2)

# 99% confidence interval for the data of @GoAHealth
CI99_A <- c(sqrt((nA-1)*var(tod_hour_A)/qchisq(0.995, df=nA-1)),
            sqrt((nA-1)*var(tod_hour_A)/qchisq(0.005, df=nA-1)))

# 99% confidence interval for the data of @ONThealth
CI99_O <- c(sqrt((nO-1)*var(tod_hour_O)/qchisq(0.995, df=nA-1)),
            sqrt((nO-1)*var(tod_hour_O)/qchisq(0.005, df=nA-1)))

thetahat_OA <- mean(tod_hour_O) - mean(tod_hour_A)

# 95% confidence interval for Testing H0: mu_A = mu_O
A2_CI95 <- c(thetahat_OA-1.96*sqrt(var(tod_hour_O)/nO+var(tod_hour_A)/nA),
          thetahat_OA+1.96*sqrt(var(tod_hour_O)/nO+var(tod_hour_A)/nA))

# p-value for Testing H0: mu_A = mu_O
p_2 <- 2*(1-pnorm(thetahat_OA/sqrt(var(tod_hour_O)/nO+var(tod_hour_A)/nA)))

x <- log(dataset$retweets+1)
y <- log(dataset$likes+1)
mod <- lm(y~x)
betahat <- 0.9810
stdres <- rstandard(mod)

plot(x,y, xlab = "retweets.log", ylab = "likes.log", 
     main = "Scatterplot of retweets.log and likes.log",
     pch = 1, cex = 0.5, col = "navy", las = 1, lwd = 1)
abline(coef(mod), lwd = 2, lty = 2, col = "red")

# standardized residuals versus the fitted values 
plot(fitted(mod), stdres, main = "Std residuals vs. fitted values", xlab = "Fitted values",
     ylab = "Standardized Residuals", pch = 1, col = "navy", cex = 0.5)
abline(h = 0, lty = 2, col = "red", lwd = 2)
x
# aa-plot of standard residual
qqnorm(stdres, main = "qqplot of std residuals", xlab = "G(0, 1) Quantiles", ylab = "Standardized Residuals",
       pch = 1, col = "navy", cex = 0.5)
qqline(stdres, lty = 2, col = "red", lwd = 2)

# 95% confidence interval for beta in analysis 3
Beta_CI95 <- c(coef(summary(mod))[2, 1] - qt(0.975, 923) * coef(summary(mod))[2, 2],
               coef(summary(mod))[2, 1] + qt(0.975, 923) * coef(summary(mod))[2, 2])

pi <- predict(mod, data.frame(x = log(30 + 1)), interval = "prediction", level = 0.90)
predictedy30x <- coef(mod)[1] + coef(mod)[2] * log(30 + 1)

# R code for Anaysis 4
longwords <- dataset$long.words
A4thetahat <- mean(longwords)

n_A4 <- length(longwords)
A4_CI95 <- c(A4thetahat-1.96*sqrt(A4thetahat/n_A4),
             A4thetahat+1.96*sqrt(A4thetahat/n_A4))
dataset$long.words[dataset$long.words >= 7] <- "7+"
observed_frequency <- table(dataset$long.words)
expected_frequency <- dpois(0:6, A4thetahat) * n_A4
expected_frequency <- append(expected_frequency, (n-sum(expected_frequency)))

testStatistic <- 2*(124*log(124/82.04443)+181*log(181/198.57773)+207*log(207/240.31563)+184*log(184/193.88412)+126*log(126/117.31774)+53*log(53/56.79043)+31*log(31/22.90896)+17*log(17/11.16096))
A4_p <- 1-pchisq(testStatistic,df=6)

barplot(rbind(observed_frequency, expected_frequency),beside=T, main = "Observed value v.s. Expected value", legend = TRUE)

# R code for Analysis 5
f1 <- sum(dataset$hashtags.binary)
f2 <- sum(dataset$media.binary)

A5_lambda <- 2*(162*log(162/121.948)+172*log(172/212.052)+175*log(175/215.02)+414*log(414/373.948))

A5_p <- 2*(1-pnorm(sqrt(A5_lambda)))

# R code for Analysis 6

mle <- dataset %>%
  group_by(username) %>%
  summarise(proportion_retweets = mean(is.retweet))

MLE <- mean(dataset$is.retweet)

nretweets <- dataset %>%
  group_by(username) %>%
  summarise(num_retweets = sum(is.retweet)*(1-MLE))

