
#############################################
# 01 - CLEAN DATA, ORGANIZE DATES AND TIMES #
#############################################

# This script has four parts:

# 1 - Setting up the environment (loading libraries, functions, etc.)
# 2 - Converting the "Date" column in the raw data to a proper Datetime object,
#     as we'll need this to understand the temporal autocorrelation in the data
# 3 - Adding a time series to the data from the beginning to the end of the 
#     sampling period
# 4 - Translating "lags" in the data into days to understand the autocorrelation

#################################
# 1 - Set up the environment ====
#################################

# Load libraries
library(tidyverse)

# Load raw data
ghg_data <- read.csv('data/ghg_data.csv')

##############################################################
# 2 - Convert the Date variable to proper datetime format ====
##############################################################

# Datetime formats make it simpler to set ordinal day, month, and calculate days
# between samples for the time series

ghg_data_dates <- ghg_data %>%
  # Adding "-" between year, month, day for date formatting, then casting as a 
  # datetime type
  mutate(Date = as.Date(paste(substr(Date, 1, 4),
                              substr(Date, 5, 6),
                              substr(Date, 7, 8), sep = '-')),
         # Adding columns for month and day from datetime
         Month = month(Date),
         Day = yday(Date)) %>%
  # Relocating so all datetime-related columns are together
  relocate(c(Month, Day), .after = Year)

# In a few instances, there were two methane measurements at the same site on a 
# single day. The time series model we will use later expects a single 
# observation per time point per group, so we need to take the mean of those two
# measurements made on the same day.
ghg_summ <- ghg_data_dates %>% 
  group_by(Site_ID, Province, Date, Day, Month, Year, Land_use) %>% 
  summarize(methane = mean(pCH4_ppm))

########################################
# 3 - Convert days to a time series ====
########################################

# Days are currently ordinal, but we need to create a variable for "Time"
# beginning at the first sample of the time series and ending at the last.

# Minimum date (earliest sample) in time series
min_date <- as.Date(min(ghg_summ$Date))

# Compute the time difference (in days) between rows of the data frame and the 
# minimum date
ghg_summ$DateDiff <- apply(ghg_summ, 1, function(x) difftime(x[3], min_date))

# Set min date to date 1, and all other dates min + 1
ghg_ts <- ghg_summ %>%
  arrange(Date) %>%
  mutate(DateDiff = round(as.numeric(DateDiff), 0)) %>%
  mutate(Time = ifelse(Date == as.Date('2021-05-10'), 1, DateDiff + 1)) %>%
  # Factor the year variable (for modelling later)
  mutate(Year_f = factor(Year))

# A quick look at the distribution of the data... They likely follow a lognormal 
# distribution. We can either specify a lognormal distribution in the model
# with a log link function for sigma, or we can transform and use a simple 
# linear regression.
ghg_ts %>%
  ggplot(aes(x = Land_use, y = methane)) +
  geom_boxplot()

# Save the cleaned time series data
saveRDS(ghg_ts, 'data/ghg_ts_data.rds')

################################################
# 4 - Calculate days ~ lags between samples ====
################################################

# We need to compute what a "lag" -- two adjacent samples -- means in terms of 
# days. When we fit the time series, it tell us how autocorrelation changes
# with unitless lags. It will be helpful to know how each lag translates to days 
# so we can interpret how autocorrelation changes with days.

# First we need to split the data into groups by Site_ID
site_grps <- ghg_summ %>%
  group_by(Site_ID) %>%
  group_split() 

# Defining a function to calculate days between rows
calc_days_btwn <- function(x) {
  diff_vec <- c()
  for(i in 1:nrow(x)) {
    if(i == 1) t <- NA
    if(i > 1) t <- as.numeric(x[3][i,] - x[3][i-1,])
    diff_vec <- c(diff_vec, t)
  }
  diff_df <- bind_cols(x, 'DaysBtwnSamples' = diff_vec)
  return(diff_df)
} 

# Now we can use the function to calculate days between lags
site_days_btwn <- lapply(site_grps, function(x) calc_days_btwn(x)) %>%
  bind_rows() 

# Summarize
lag_summary <- site_days_btwn %>%
  summarize(Days = unique(DaysBtwnSamples)) %>%
  arrange(Days) %>%
  rownames_to_column('Lag') %>%
  # Convert the lag column to numeric
  mutate(Lag = as.numeric(Lag))

# How many sample intervals were < 16 days? (in the next script we'll see this 
# is approximately the point where temporal autocorrelation ceases to be 
# significant)
site_days_btwn %>%
  filter(DaysBtwnSamples < 16) %>%
  summarize(n())

# Save the lag summary data
saveRDS(lag_summary, 'data/lag_summary.rds')
