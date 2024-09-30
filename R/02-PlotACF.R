
#######################################
# 02 - PLOT AUTOCORRELATION FUNCTION #
######################################

# This script has four parts:

# 1 - Setting up the environment (loading libraries, functions, etc.)
# 2 - Computing the confidence intervals for autocorrelation, which will be used
#     to plot boundaries to visualize when the temporal autocorrelation is 
#     significant or not
# 3 - Computing the autocorrelation function by group
# 4 - Creating and saving the plot

#################################
# 1 - Set up the environment ====
#################################

# Load libraries
library(tidyverse)

# Load cleaned data
ghg_ts <- readRDS('data/ghg_ts_data.rds')

# Load time lags ~ days (for translating)
lag_data <- readRDS('data/lag_summary.rds')

####################################################
# 2 - Compute confidence intervals for ACF plot ====
####################################################

# Get max length of site time series in days. This is a parameter needed to
# compute the confidence intervals for the ACF plot below.
max_n <- ghg_ts %>%
  group_by(Site_ID) %>%
  summarize(N = n()) %>%
  pull(N) %>%
  max()

# Compute the widest confidence intervals for autocorrelation. The formula to 
# calculate the confidence intervals for an ACF plot time series is Var(rk)∼1/N,
# where r is the correlation at lag k and N is the number of samples. Because 
# the data are grouped by site, each site only has a fraction of the total
# number of samples. We will calculate the CIs for a "worst case scenario", the
# maximum number of samples taken at a single site.
max_ci <- qnorm((1 + .95)/2)/sqrt(max_n)

######################################
# 3 - Build the ACF plot by group ====
######################################

# Manipulating the data so we can compute ACF at each lag for each group 
# (Site_ID) in the data frame.
grouped_acf_values <- ghg_ts %>%
  # Log-transforming methane
  mutate(methane = log(methane)) %>%
  # Arrange the time series in order
  arrange(Time) %>%
  # Create a list column of data frames in our original data frame for computing
  # acf by group
  nest(-Site_ID) %>%
  mutate(acf_results = map(data, ~ acf(.x$methane, lag.max = 15, plot = F)),
                ACF_values = map(acf_results, ~ drop(.x$acf))) %>%
  # Then unnest back into a column
  unnest(ACF_values) %>%
  group_by(Site_ID) %>%
  # Compute the lag
  mutate(Lag = seq(0, n() - 1)) %>%
  # Clean up the data frame
  select(! c(data, acf_results))

# We really only need to plot the max and min ACF values for each lag to show a 
# "worst case scenario" of autocorrelation. So, we can just get the max and min
# ACF by lag.
max_acf_lag <- grouped_acf_values %>%
  group_by(Lag) %>%
  summarize(max_ACF = max(ACF_values),
            min_ACF = min(ACF_values)) %>%
  # Now joining the lag data to translate lags into days
  left_join(lag_data)

# Set the first day to zero and minimum ACF = -1
max_acf_lag[1,]$Days <- 0
max_acf_lag[1,]$min_ACF <- -1

##########################
# 4 - Create the plot ====
##########################

# Plot the ACF plot
max_acf_lag %>% 
  ggplot(aes(x = factor(Days))) + 
  geom_rect(aes(ymin = -max_ci, ymax = max_ci,
                xmin = -Inf, xmax = Inf), fill = 'grey90') +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = max_ci, linetype = 'dashed') + 
  geom_hline(yintercept = -max_ci, linetype = 'dashed') +
  geom_segment(aes(x = factor(Days), y = min_ACF, 
                   xend = factor(Days), yend = max_ACF),
               linewidth = 3, lineend = 'round', colour = '#1560bd') +
  labs(y = 'Autocorrelation in methane (log ppm)', x = 'Days') +
  theme(plot.margin = unit(c(0.5, 0.5, 1, 1), 'cm'),
        panel.background = element_rect(fill = 'white', colour = 'black',
                                        linewidth = 1),
        panel.grid = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -4),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 4))

# Save the plot
ggsave('acf_plot.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, 
       height = 12, width = 20, units = 'cm', bg = 'white')

