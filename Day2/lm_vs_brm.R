#Simulating data for linear regression

cookies <- c(0, 3, 5, 7, 9, 10, 13, 16, 18, 20)
gain <- c(-0.1, 0.4, -0.7,  1.7, 1.3, 0.2, 0.8, 2.5, 1.0, 2.8)
df_cookies <- data.frame(gain, cookies)

plot(cookies, gain, pch = 16, col = "blue", cex = 1.5)
abline(lm(gain ~ cookies))

#OLS regression

m_ols <- lm(gain ~ cookies, data = df_cookies)
summary(m_ols)

confint(m_ols)

library(effects)
plot(effect("cookies", m_ols))

df_cookies$sports <- c("Yes", "No", "Yes", "No", "No", "Yes", "No", "No", "Yes", "No")

m_ols1 <- update(m_ols, . ~ . + sports)
anova(m_ols, m_ols1)
AIC(m_ols)
AIC(m_ols1)

plot(m_ols1)

library(car)
vif(m_ols1)

fitted(m_ols1)[1:10]


#Bayesian linear model

library(brms)
m_brm <- brm(gain ~ cookies, data = df_cookies, chains = 1, iter = 20)
summary(m_brm)

plot(m_brm)

as_draws_df(m_brm)$b_Intercept
as_draws_df(m_brm)$b_cookies
as_draws_df(m_brm)$sigma

m_brm1 <- brm(gain ~ cookies, data = df_cookies, chains = 1, iter = 2000)
plot(m_brm1)

m_brm2 <- brm(gain ~ cookies, data = df_cookies, chains = 4, iter = 2000)
plot(m_brm2)

posterior_cookies <- as_draws_df(m_brm2)$b_cookies
mean(posterior_cookies)
median(posterior_cookies)
sd(posterior_cookies) 

quantile(posterior_cookies, c(0.025, 0.975)) 
library(bayestestR)
hdi(posterior_cookies, 0.95)

plot(density(posterior_cookies))

mean(posterior_cookies > 0)

m_brm3 <- update(m_brm2, formula = . ~ . + sports, newdata = df_cookies)

m_brm1 <- add_criterion(m_brm1, criterion = "waic")
m_brm3 <- add_criterion(m_brm3, criterion = "waic")
loo_compare(m_brm1, m_brm3, criterion="waic")

m_brm1 <- add_criterion(m_brm1, criterion = "loo")
m_brm3 <- add_criterion(m_brm3, criterion = "loo")
loo_compare(m_brm1, m_brm3, criterion="loo")

bayes_R2(m_brm1) 
bayes_R2(m_brm3)

fitted(m_brm3)

pp_check(m_brm1)

res <- residuals(m_brm)
head(res)

get_prior(gain ~ cookies + sports, data = df_cookies)

median(gain)

prior_summary(m_brm3)

x <- seq(-10, 10, by = 0.1)
plot(x, dnorm(x, mean = 0, sd = 1),type="l")
lines(x, dnorm(x, mean = 0, sd = 2), col = "blue")
lines(x, dnorm(x, mean = 1, sd = 3), col = "red")

plot(x, dcauchy(x, location = 0, scale = 1), type = "l")
lines(x, dcauchy(x, location = -1, scale = 2), col = "blue")
lines(x, dcauchy(x, location = 2, scale = 3), col = "red")

y <- rt(100000, df = 3)*2.5 + 0.9 
#the values from our model
plot(density(y), xlim = c(-10, 10), ylim = c(0, 0.2))
y <- rt(100000, df = 3)*2 + 0
lines(density(y), col = "blue")
y <- rt(100000, df = 4)*5 - 1
lines(density(y), col = "red")

summary(m_brm3)

my_priors1 <- prior(normal(0, 1), class = "b")
m_brm4 <- brm(gain ~ cookies + sports, data = df_cookies, prior = my_priors1)
summary(m_brm4)

my_priors2 <- prior(normal(3, 1), class = "b")
m_brm5 <- brm(gain ~ cookies + sports, data = df_cookies, prior = my_priors2)
summary(m_brm5)

my_priors3 <- prior(cauchy(3, 2), class = "b")
m_brm6 <- brm(gain ~ cookies + sports, data = df_cookies, prior = my_priors3)
summary(m_brm6)

my_priors4 <- c(prior(normal(0, 1), class = "b", coef = "cookies"), prior(cauchy(2, 1), class = "b", coef = "sportsYes"))     
m_brm7 <- brm(gain ~ cookies + sports, data = df_cookies, prior = my_priors4)
summary(m_brm7)
