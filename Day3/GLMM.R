help_df <- read.table("help_big.txt", header = T, stringsAsFactors = TRUE)
str(help_df)

###Fitting a Maximum Likelihood GLMM 

library(lme4)

help_glmer <- glmer(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + (1|Verb), data = help_df, family = binomial)
#convergence issues
summary(help_glmer)

help_glmer <- glmer(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + (1|Verb), data = help_df, family = binomial, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(help_glmer)

glmer_predict <- predict(help_glmer, type = "response")
length(glmer_predict)

library(Hmisc)
somers2(glmer_predict, as.numeric(help_df$Response) - 1)

###Fitting a Bayesian GLMM
library(brms)

my_formula <- bf(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + (1|Verb))

#weakly informative priors

get_prior(formula = my_formula, data = help_df, family = bernoulli)

my_priors <- prior(cauchy(0, 2), class = b) + 
  prior(normal(0, 1), class = sd)

help_brm <- brm(formula = my_formula, data = help_df, prior = my_priors, family = bernoulli, cores = 4)

summary(help_brm)

cbind(glmer = fixef(help_glmer), brm = fixef(help_brm)[, 1])


pp_check(help_brm)
pp_check(help_brm, type = "error_binned")

preds <- posterior_linpred(help_brm, newdata = help_df, transform = TRUE)
dim(preds)
preds[1:5, 1:5]

c_scores <- numeric(4000)

library(pROC)

for (i in 1:4000) {
  c_obj <- roc(response = help_df$Response, predictor = preds[i, ], quiet = TRUE)
  c_scores[i] <- auc(c_obj)
}

mean(c_scores)
quantile(c_scores, probs = c(0.025, 0.975))


###Small sample

help_df_small <- read.table("help_small.txt", header = TRUE, stringsAsFactors = TRUE)
 
str(help_df_small)

#Small sample: GLMER model

help_glmer_small <- glmer(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + Stress + (1|Verb), data = help_df_small, family = binomial, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))

summary(help_glmer_small)

#Small sample: Bayesian model with flat (default) priors

help_brm_flat <- brm(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + Stress + (1|Verb), data = help_df_small, family = bernoulli, cores = 4)
summary(help_brm_flat)

#Small sample: Bayesian model with weakly informative priors

help_brm_weak <- brm(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + Stress + (1|Verb), data = help_df_small, family = bernoulli, prior = my_priors, cores = 4)
summary(help_brm_weak)

#Small sample: Bayesian model with strong informative priors for the problematic interaction, weakly informative for others

my_priors_strong <- prior(cauchy(0, 2), class = "b") +  prior(normal(4.5, 1.5), class = "b", coef = "HorrorYes:Distance_log")
help_brm_strong <- brm(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + Stress + (1|Verb), data = help_df_small, family = bernoulli, prior = my_priors_strong, cores = 4)
summary(help_brm_strong)

#Compare the estimates

cbind(big_glmer = fixef(help_glmer), 
      big_brms = fixef(help_brm)[, 1],
      small_glmer = fixef(help_glmer_small)[-c(9:10)],
      small_brm_flat = fixef(help_brm_flat)[-c(9:10), 1], 
      small_brm_weak = fixef(help_brm_weak)[-c(9:10), 1], 
      small_brm_strong = fixef(help_brm_strong)[-c(9:10), 1])