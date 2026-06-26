library(brms)

#flat priors for cookies
m_brm_flat0 <- brm(gain ~ 1, data = df_cookies, cores = 4, save_pars = save_pars(all = TRUE))
m_brm_flat1 <- brm(gain ~ cookies, data = df_cookies, cores = 4, save_pars = save_pars(all = TRUE))

logml_flat0 <- bridge_sampler(m_brm_flat0)
logml_flat0

logml_flat1 <- bridge_sampler(m_brm_flat1)
logml_flat1

exp(logml_flat0$logml)/exp(logml_flat1$logml)
bayes_factor(logml_flat0, logml_flat1)

exp(logml_flat1$logml)/exp(logml_flat0$logml)
bayes_factor(logml_flat1, logml_flat0)

#weakly informative priors on cookies
my_priors_weak <- prior(normal(0, 10), class = "b")
m_brm_weak0 <- brm(gain ~ 1, data = df_cookies, cores = 4, save_pars = save_pars(all = TRUE))
m_brm_weak1 <- brm(gain ~ cookies, data = df_cookies, cores = 4, prior = my_priors_weak, save_pars = save_pars(all = TRUE))
bayes_factor(m_brm_weak1, m_brm_weak0)
bayes_factor(m_brm_weak0, m_brm_weak1)

#strong specific prior on cookies
my_priors_strong <- prior(normal(0.10, 0.1), class = "b")
m_brm_strong0 <- brm(gain ~ 1, data = df_cookies, cores = 4, save_pars = save_pars(all = TRUE))
m_brm_strong1 <- brm(gain ~ cookies, data = df_cookies, cores = 4, prior = my_priors_strong, save_pars = save_pars(all = TRUE))
bayes_factor(m_brm_strong0, m_brm_strong1)
bayes_factor(m_brm_strong1, m_brm_strong0)

###mixed model

#flat priors on University and Vocab_Score
mm_brm0 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score, data = df, cores = 4, save_pars = save_pars(all = TRUE))
mm_brm1 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df, cores = 4, save_pars = save_pars(all = TRUE))
bayes_factor(mm_brm0, mm_brm1)
bayes_factor(mm_brm1, mm_brm0)

#weak generic priors on University and Vocab_Score
my_prior_university <- prior(normal(0, 10), class = "b", coef = "University1")
my_prior_vocab <- prior(normal(0, 10), class = "b", coef = "Vocab_Score")
mm_brm_weak0 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score, data = df, cores = 4, prior = my_prior_vocab, save_pars = save_pars(all = TRUE))
mm_brm_weak1 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df, cores = 4, prior = my_prior_vocab + my_prior_university, save_pars = save_pars(all = TRUE))
bayes_factor(mm_brm_weak0, mm_brm_weak1)
bayes_factor(mm_brm_weak1, mm_brm_weak0)

#strong specific priors on University and Vocab_Score
my_prior_university <- prior(normal(8, 1), class = "b", coef = "University1")
my_prior_vocab <- prior(normal(0.1, 0.1), class = "b", coef = "Vocab_Score")
mm_brm_strong0 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score, data = df, cores = 4, prior = my_prior_vocab, save_pars = save_pars(all = TRUE))
mm_brm_strong1 <- brm(Nwords ~ (1|Student_ID) + (1|Topic_ID) + Vocab_Score + University, data = df, cores = 4, prior = my_prior_vocab + my_prior_university, save_pars = save_pars(all = TRUE))
bayes_factor(mm_brm_strong0, mm_brm_strong1)
bayes_factor(mm_brm_strong1, mm_brm_strong0)
