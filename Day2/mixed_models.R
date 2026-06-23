###This file contains code that will teach you how to create simulation data. 

#Background: inspired by Jongman et al. (2021)
#Language production data from 70 university and 30 vocational college students.
#The response variable is the number of words they produced in 1 minute.
#The task was, for example, to tell about how they spent their weekend.
#Below you find instructions how to simulate data for different situations,
#and how to fit the "true" models.

#Every student produces a response on ten topics.
#This means that the data are no longer independent.
#We need random intercepts for Student and Topic.

set.seed(21)
stud_intercepts <- rnorm(n = 100, mean = 0, sd = 10)
hist(stud_intercepts)
head(stud_intercepts)

set.seed(66)
top_intercepts <- rnorm(n = 10, mean = 0, sd = 7)
top_intercepts

#This time we have 1000 values in each of the vectors

stud_intercepts_long <- rep(stud_intercepts, times = 10)
stud_intercepts_long[c(1, 2, 101, 102, 201, 202)]
top_intercepts_long <- rep(top_intercepts, each = 100)
top_intercepts_long[1:5]

stud_id <- 1:100
top_id <- 1:10

stud_id_long <- rep(stud_id, times = 10)
stud_id_long[c(1, 2, 101, 102, 201, 202)]

top_id_long <- rep(top_id, each = 100)
top_id_long[1:5]

intercept <- 125


set.seed(154)
vocab_score <- sample(1:99, 100, replace = TRUE)
hist(vocab_score)
vocab_score_long <- rep(vocab_score, times = 10)

university <- c(rep(1, times = 70), rep(0, times = 30))
university_long <- rep(university, times = 10)

set.seed(77)
noise <- rnorm(n = 1000, mean = 0, sd = 5)

nwords <- intercept + stud_intercepts_long + top_intercepts_long + 0.15*vocab_score_long + 7*university_long + noise
nwords <- round(nwords)
hist(nwords)

df <- data.frame(Student_ID = stud_id_long, Topic_ID = top_id_long, 
                 Nwords = nwords, Vocab_Score = vocab_score_long,
                 University = university_long)

df$Student_ID <- as.factor(df$Student_ID)
df$Topic_ID <- as.factor(df$Topic_ID)
df$University <- as.factor(df$University)
str(df)        

#Restricted Maximum Likelihood model
        
library(lme4)

m_lmer <- lmer(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df)
summary(m_lmer)

#To get p-values:
library(lmerTest)

mm_lmer <- lmer(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df)
summary(mm_lmer)
#Are the residuals, random effects and fixed effects similar to what we created?


#Bayesian model
#default priors, default MCMC
library(brms)
mm_brm <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df, cores = 4)
summary(mm_brm)

#Compare the results

fixef(mm_lmer)
fixef(mm_brm)

ranef(mm_lmer)
ranef(mm_brm)

get_prior(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df)

#Exercise
conditional_effects(mm_brm)

pp_check(mm_brm)

bayes_R2(mm_brm)

mm_brm1 <- update(mm_brm, formula = . ~ . - University, newdata = df, cores = 4)

mm_brm <- add_criterion(mm_brm, criterion = "waic")
mm_brm1 <- add_criterion(mm_brm1, criterion = "waic")
loo_compare(mm_brm, mm_brm1, criterion="waic")

#specify: my_priors_mm <- c(prior(normal(..., ...), class = "b", coef = "University1"))
#mm_brm2 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df, prior = my_priors_mm, cores = 4)

