#------------------------------------------------------------------------------
# Script uses data reported in Steve Smith's manuscript and recreates his
# power analysis. 
#
# Author: Ryan N. Kinzer
# Created: 02/20/2018
#------------------------------------------------------------------------------
# load packages
#------------------------------------------------------------------------------
library(tidyverse)
library(lme4)
library(MASS)

#------------------------------------------------------------------------------
# load fish and cohort spread factors
#------------------------------------------------------------------------------
load('./data/data.rda')

#------------------------------------------------------------------------------
# select species and origin group of interest
#------------------------------------------------------------------------------
species <- c('Chinook', 'Steelhead')
origin <- c('both', 'hatchery', 'wild')

spp <- species[1]
org <- origin[1]

#------------------------------------------------------------------------------
# get retrospective data based on spp and org
#------------------------------------------------------------------------------
retro_df <- fish_list %>%
  filter(species == spp & origin == org) %>%
  unnest() %>%
  rename(A1 = A) %>%
  mutate(A0 = J - A1,
         gas_cap = as.factor(rep(0, n())),
         obs = as.factor(1:n()),
         period = factor(period))

#------------------------------------------------------------------------------
# get cohort-spread factors based on spp
#------------------------------------------------------------------------------

tmp <- spread_list %>% # probably a better/cleaner way to do this.
  filter(species == spp) %>%
  unnest()

spread_factors <- list('x0011' = as.tibble(tmp$data[1][[1]]),
                'x1100' = as.tibble(tmp$data[2][[1]]))


#------------------------------------------------------------------------------
# set power analysis parameters
#------------------------------------------------------------------------------
delta <- 1.50 # percent increase in SAR only takes one number (need to loop)
Yp <- c(2,4,6) # number of prospective years (takes a vector of years)
cj <- 1.0 # sample size multiplier (only takes one number)
iter <- 100

#------------------------------------------------------------------------------
# run power analysis 
#------------------------------------------------------------------------------
sim_results <- spill_power(delta, Yp, cj, retro_df, spread_factors, iter)


# combine all results into one data frame
sim_df <- bind_rows(sim_results, .id = 'sim_data') %>% 
  gather(key = 'year', value = 'value', -sim_data) %>%
  separate(sim_data, into = c('design', 'metric'), sep = '_')

# pull out only p-values
p_df <- sim_df %>%
  filter(metric == 'pvalue') %>%
  mutate(delta = delta,
         alpha_05 = ifelse(value <= .05, 1, 0), 
         alpha_10 = ifelse(value <= .1, 1, 0),
         year = as.numeric(year))

# estimated power
p_df %>%
  group_by(delta, design, year) %>%
  summarise(power_05 = sum(alpha_05)/n(),
            power_10 = sum(alpha_10)/n())

        