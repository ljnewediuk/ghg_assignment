

# This script has four parts:

# 1 - Setting up the environment (loading libraries, functions, etc.)
# 2 - Fitting the models. We are fitting the same causal model but using
#     WAIC to select the parameters for the ARMA model.
# 3 - Performing posterior checks and goodness of fit tests on the final model
# 4 - Predicting from the model and plotting predictions from the model plus
#     the conditional effects of land use accounting for other variables in the
#     model

#################################
# 1 - Set up the environment ====
#################################

# Load libraries
library(brms)
library(tidyverse)

# Also required but not loaded: bayesplot, bayestestR, performance, cowplot

# Load cleaned data
ghg_ts <- readRDS('data/ghg_ts_data.rds')

################################################
# 2 - Fit models and select ARMA parameters ====
################################################

# Each model is fit with the following terms:
#   - log-transformed methane emissions (main effects)
#   - year to account for annual variation (main effects)
#   - site nested in province to account for repeated measures (random effects)
#   - ARMA term to account for temporal autocorrelation

# Each model is fit with the following specs:
#   - Weakly informative priors (ß ~ N(0, 10) and sigma ~ C(0, 1))
#   - Gaussian family with identity link (i.e., linear regression)
#   - 2000 iterations with 1000 warmup iterations
#   - 4 chains and 4 cores

# First I am performing model selection to select parameters "p" and "q" for the
# ARMA term. "p" is the order of the auto-regressive (AR) part of the model,
# and "q" is the order of the moving average (MA) part of the model. It is 
# difficult to choose these parameters, so I am using WAIC (a Bayesian analog
# of AIC) to compare model fit with all combinations of p and q from 0--5, while 
# keeping the rest of the causal model the same.

# Use expand_grid() to get all combos of p and q from 0 to 5 and names m1 to mn. 
# The first row must be removed where both p and q are zero.
m_tbl <- expand_grid(p = 0:5, q = 0:5)[-1,]
m_tbl$name <- paste0('m', 1:nrow(m_tbl))

# Now, loop through fitting models with the parameter combinations in each row 
# of "m_tbl". WARNING: This is an extremely slow step. Fitting all 35 models
# took approximately 3 hours on a 3.1 GHz Dual-Core Intel Core i5 MacBook Pro
# with 8 GB memory.
m_list <- list()
for(i in 1:nrow(m_tbl)) {
  # Specify p and q from row i
  p_i <- m_tbl[i,]$p
  q_i <- m_tbl[i,]$q
  # Fit model
  m <- brm(data = ghg_ts,
           family = gaussian,
           # Specifying terms (land use, year, site ID/province, and ARMA)
           bf(log(methane) ~ Land_use + Year_f + (Land_use | Site_ID/Province) +
                arma(time = Time, gr = Site_ID, p = p_i, q = q_i)),
           # Setting priors
           prior=c(prior(normal(0,10), class=b), 
                   prior(cauchy(0,1), class=sigma)),
           # Model specs
           iter=2000, warmup=1000, chains=4, cores=4,
           backend = 'cmdstanr')
  # Add the model to the model list and name it for reference
  m_list[[i]] <- m
  names(m_list)[[i]] <- m_tbl[i,]$name
}

# Save model list
saveRDS(m_list, 'output/ARMA_param_select_mods.rds')

# Compute WAICs of the model based on the posterior likelihood in list
waic_list <- lapply(m_list, function(x) waic(x))

# Save WAIC list
saveRDS(waic_list, 'output/ARMA_WAICs.rds')

# loo_compare() for best combination of p and q. The best-fitting model is the
# top model, but large SEs and differences ≥ 4 are not really substantial.
loo_compare(waic_list)

# The results suggest m1 and m2 are comparable, but m1 is lowest, so we will go 
# with m1. This means p = 0 and q = 1.

# Save m1 as the final model
saveRDS(m_list$m1, 'output/final_ghg_model.rds')

#########################################################
# 3 - Perform model checks and goodness of fit tests ====
#########################################################

# First, we will just pass the model to the plot function. This draws "trace
# plots", which just show the model exploration of the posterior. There should
# be a "fuzzy caterpillar" pattern with all chains overlapping and indistinct.
# These trace plots look fine. The function also shows us the parameter
# posteriors.
plot(m_list$m1)

# Posterior predictive check, indicating whether simulations from the model 
# posterior (yrep) predict the data (y). Spaghetti lines should overlap with y 
# if the model predicts the data well. 

# Extracting the data and predicting from the posterior
y <- m_list$m1$data[,1]
yrep <- posterior_predict(m_list$m1, ndraws = 100)

# Setting the colour scheme
bayesplot::color_scheme_set('brightblue')

# Plotting the posterior predictive check
pp_plot <- bayesplot::ppc_dens_overlay(y, yrep) + 
  theme(plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black',
                                        linewidth = 1),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -4),
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.position = 'inside',
        legend.position.inside = c(0.8, 0.8))

# Save the plot
ggsave('pp_check.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, 
       height = 12, width = 15, units = 'cm', bg = 'white')

# The probability of a directional effect, analogous to the frequentist p-value.
# p_direction ≥ 95% can be thought of as a "significant" p-value ≤ 0.05. For
# land use the value is 98%, suggesting a relatively strong relationship between
# land use and methane emissions.
bayestestR::p_direction(m_list$m1)

# The Bayesian R-squared. The interpretation is again similar to the that of a
# frequentist R-squared, but instead of variance ratio between the predicted 
# values and data, the Bayesian version is the variance ratio between the 
# predicted values and the variance of predicted values plus the expected error. 
# The model explains a decent amount of the variation in the data.
performance::r2_bayes(m_list$m1)

###########################################
# 4 - Predict from the model and plot ====
##########################################

# To summarize the model, we need to predict "draws" from the posterior. Here we 
# are taking 4000 draws.
p1 <- posterior_epred(m1, summary = F, ndraws = 4000, re_formula = NA)

# Now we can summarize the posterior for plotting. We are computing the mean of 
# the posterior predictive draws and upper and lower 95% credible intervals and 
# adding them to the data frame.
ghg_ts$model_mean <- apply(p1, 2, mean)
ghg_ts$model_lower <- apply(p1, 2, function(x) quantile(x, probs = 0.025))
ghg_ts$model_upper <- apply(p1, 2, function(x) quantile(x, probs = 0.975))

# Plot the time series
ts_plot <- ghg_ts %>%
  ggplot(aes(x = Day, y = log(methane), colour = Land_use)) + 
  geom_point(size = 1.5) +
  scale_colour_manual(values = c('#a4e996', '#ffc878')) +
  scale_fill_manual(values = c('#a4e996', '#ffc878')) +
  geom_ribbon(aes(x = Day, ymin = model_lower, ymax = model_upper, fill = Land_use),
              alpha = 0.4, colour = NA) +
  geom_line(aes(x = Day, y = model_mean), linewidth = .5) + 
  facet_wrap(~Year, scales = 'free', ncol = 1) +
  ylim(0,10) +
  labs(y = 'Methane concentration (log ppm)', x = 'Ordinal day') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black',
                                        linewidth = 1),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 4),
        strip.background = element_rect(fill = 'white', colour = 'black',
                                        linewidth = 1),
        strip.text = element_text(size = 15, colour = 'black'),
        legend.position = 'none')

# Next, we will pull out the conditional effects of land use on methane 
# emissions accounting for the other variables in the model.
ce1 <- data.frame(Land_use = conditional_effects(m_list$m1)$Land_use$Land_use,
                  beta = conditional_effects(m_list$m1)$Land_use$`estimate__`,
                  lower = conditional_effects(m_list$m1)$Land_use$`lower__`,
                  upper = conditional_effects(m_list$m1)$Land_use$`upper__`)

# Now we can plot the conditional effects
ci_plot <- ggplot(ce1, aes(x = Land_use)) +
  geom_segment(aes(x = Land_use, xend = Land_use, y = lower, yend = upper),
               linewidth = 2, lineend = 'round') +
  geom_point(aes(x = Land_use, y = beta), size = 5) +
  labs(y = 'Methane emissions (log ppm)', x = 'Land use') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black',
                                        linewidth = 1),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -4),
        axis.title.y = element_blank())

# Make a panel plot
cowplot::plot_grid(ts_plot, ci_plot, labels = "AUTO", label_size = 18, rel_widths = c(1, 1))

# Save the plot
ggsave('methane_effects.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, 
       height = 15, width = 28, units = 'cm', bg = 'white')

