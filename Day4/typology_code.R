#Packages we'll need for this session
library(brms)
library(ape) #install it first if needed
remotes::install_github("erichround/glottoTrees",
                        dependencies = TRUE)

library(glottoTrees)

data <- read.table("data_typology.txt", header = T, sep = "\t")

#Is there a correlation between the proportion
#of overlapping S and O forms in a language, 
#and Shannon’s entropy of SO order?

plot(data$SO_Overlap, data$Entropy)
#The plot suggests a negative correlation. But we need a proper statistical test to be sure.

###Naive, anti-Galtonian model###

#We first fit a Bayesian model without any genealogical or geographic information. 
#Before doing so, we define weakly informative priors.

my_priors_fixed <- set_prior("normal(0, 5)", class = "b")

brm_naive <- brm(Entropy ~ SO_Overlap, data = data, cores = 4, iter = 4000, prior = my_priors_fixed)
summary(brm_naive)

conditional_effects(brm_naive)

bayes_R2(brm_naive)

pp_check(brm_naive)

#What's the problem?

#Entropy values for a binomial outcome (SO or OS) 
#cannot be greater than 1 and less than 0.
#We need to switch to beta regression.
#But we first transform the Entropy variable 
#because beta regression cannot deal with values 
#that are strictly equal to 0 or 1.

summary(data$Entropy)
data$Entropy1 <- data$Entropy
data$Entropy1[data$Entropy == 0] <- 0.001
summary(data$Entropy1)

#Let us fit a Beta model. 
#We will use the same priors as above, 
#but note that in many cases you will need to adjust the priors.

brm_naive_beta <- brm(Entropy1 ~ SO_Overlap, data = data, cores = 4, iter = 4000, 
                 prior = my_priors_fixed, family = Beta)
summary(brm_naive_beta)

conditional_effects(brm_naive_beta)

#What has changed?
pp_check(brm_naive_beta)
bayes_R2(brm_naive_beta)

###Mixed-effects model with Genus###

#Let us try to take into account the phylogenetic information by adding random intercepts for every Genus. 
#It is also good to provide priors for the random effects.
#Our priors are quite generous.

my_priors_mixed <- c(set_prior("normal(0, 5)", class = "b"), 
                         set_prior("normal(0, 2)", class = "sd"))

brm_genus <- brm(Entropy1 ~ SO_Overlap + (1|genus), data = data, family = Beta, cores = 4, iter = 4000, 
                 prior = my_priors_mixed)
summary(brm_genus)
plot(brm_genus)
pp_check(brm_genus)
bayes_R2(brm_genus)

#Which model is better? Let us use the LOO information criterion for model comparison.

brm_naive_beta <- add_criterion(brm_naive_beta, "loo")
brm_genus <- add_criterion(brm_genus, "loo")

#To fix the problem indicated in the warning message, 
#we refit the model brm_genus with save_pars = save_pars(all = TRUE) 
#and then improve the LOO estimation. 

brm_genus <- brm(Entropy1 ~ SO_Overlap + (1|genus), 
                 family = Beta,
                 data = data, cores = 4, iter = 4000, 
                 prior = my_priors_mixed, 
                 save_pars = save_pars(all = TRUE))

brm_genus <- add_criterion(brm_genus, "loo", moment_match = TRUE)
#This takes a while...

#Now we are ready to compare the models again:
loo_compare(brm_naive_beta, brm_genus, criterion = "loo")

#Which model is better? Why?

###Phylogenetic model###

#Our first phylogenetic model will include a variance-covariance matrix 
#with covariances (correlations) between all languages in the sample. 
#To create the matrix, we will use some functions from the package glottoTrees:

my_supertree <- assemble_supertree(macro_groups = NULL)
my_supertree <- abridge_labels(my_supertree)
my_tree <- keep_as_tip(my_supertree, label = data$glottocode)
plot_glotto(my_tree)

#see more information in this tutorial: 
#https://ladal.edu.au/tutorials/phylogenetic_methods/phylogenetic_methods.html


my_cor <- ape::vcv.phylo(my_tree, corr = TRUE)
print(my_cor[1:10, 1:10])


brm_phylo <- brm(Entropy1 ~ SO_Overlap + 
                 (1 |gr(glottocode, cov = my_cor)),
                 data = data, 
                 data2 =list(my_cor = my_cor), 
                 family = Beta, cores = 4, 
                 prior = my_priors_mixed,
                 iter = 4000, control = list(adapt_delta = 0.9999, max_treedepth = 20), 
                 save_pars = save_pars(all = TRUE))

summary(brm_phylo)
plot(brm_phylo)

pp_check(brm_phylo)
bayes_R2(brm_phylo)

brm_phylo <- add_criterion(brm_phylo, "loo", moment_match = TRUE)
loo_compare(brm_naive_beta, brm_genus, brm_phylo, criterion = "loo")

#There are seem to be some numerical problems with brm_phylo.
#It can be useful to use higher-level units.
#We'll take Genera. The difference from brm_genus is that
#we also take into account the tree information above the Genus level.

my_tree1 <- keep_as_tip(my_supertree, label = data$glottocode_genus)
plot_glotto(my_tree1)
my_cor1 <- ape::vcv.phylo(my_tree1, corr = TRUE)

brm_phylo1 <- brm(Entropy1 ~ SO_Overlap + 
                 (1 |gr(glottocode_genus, cov = my_cor1)),
                 data = data, 
                 data2 =list(my_cor1 = my_cor1), cores = 4, 
                 family = Beta,
                 prior = my_priors_mixed,
                 iter = 4000, control = list(adapt_delta = 0.99), 
                 save_pars = save_pars(all = TRUE))

summary(brm_phylo1)
plot(brm_phylo1)
pp_check(brm_phylo1)
bayes_R2(brm_phylo1)

#But what if we exclude SO_Overlap?

my_priors_random <- set_prior("normal(0, 2)", class = "sd")

brm_phylo_only <- brm(Entropy1 ~ 
                    (1 |gr(glottocode, cov = my_cor)),
                  data = data, 
                  data2 = list(my_cor = my_cor), cores = 4, 
                  family = Beta,
                  prior = my_priors_random,
                  iter = 4000, control = list(adapt_delta = 0.99), 
                  save_pars = save_pars(all = TRUE))

summary(brm_phylo_only)
plot(brm_phylo_only)
pp_check(brm_phylo_only)
bayes_R2(brm_phylo_only)

#Which of the five models is the best, according to LOOIC?

brm_phylo1 <- add_criterion(brm_phylo1, "loo", moment_match = TRUE) #this can take a while
brm_phylo_only <- add_criterion(brm_phylo_only, "loo", moment_match = TRUE) #this can take a while


loo_compare(brm_naive_beta, brm_genus, brm_phylo, brm_phylo1, brm_phylo_only, criterion = "loo")



##########################################################
###Useful  code for exploring gamma (or any other) prior distributions
#Modified from https://solomonkurz.netlify.app/blog/2023-06-25-causal-inference-with-beta-regression/

library(ggdist)

df <- rbind(parse_dist(prior(gamma(0.01, 0.01))),  
            parse_dist(prior(gamma(2, 0.1))),
            parse_dist(prior(gamma(4, 0.5))))  

ggplot(df, aes(xdist = .dist_obj, y = prior)) + 
  stat_halfeye(.width = c(.5, .99), p_limits = c(.0001, .9999)) +
  scale_x_continuous(expression(italic(p)(phi)), breaks = 0:4 * 25) +
  scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
  labs(title = "Prior distributions for phi parameter") +
  coord_cartesian(xlim = c(0, 75))

#Normal

df <- parse_dist(prior(normal(0, 5)))  

ggplot(df, aes(xdist = .dist_obj, y = prior)) + 
  stat_halfeye(.width = c(.5, .99), p_limits = c(.0001, .9999)) +
  scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
  labs(title = "My priors") + xlim(-15, 15)


#Student t

df <- parse_dist(prior(student_t(3, 0.9, 2.5)))  

ggplot(df, aes(xdist = .dist_obj, y = prior)) + 
  stat_halfeye(.width = c(.5, .99), p_limits = c(.0001, .9999)) +
  scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
  xlim(-20, 20) + labs(title = "My priors")


#Cauchy

df <- parse_dist(prior(cauchy(0, 2.5)))  

ggplot(df, aes(xdist = .dist_obj, y = prior)) + 
  stat_halfeye(.width = c(.5, .99), p_limits = c(.0001, .9999)) +
  scale_y_discrete(NULL, expand = expansion(add = 0.1)) +
  labs(title = "My priors") + xlim(-20, 20)
