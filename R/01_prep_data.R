#------------------------------------------------------------------------------
# Script loads two .csv files.  The first file contains juvenile and adult 
# return data transcribed from Steve Smith's power analysis document (1/18/18).
# The second .csv file loads Steve's reported cohort spread factors.  After
# loading, the data is transformed into formats that are accessible for the
# simulation and power analysis functions.  Data is then stored in lists and
# saved to a .rda file.
#
# Author: Ryan N. Kinzer
# Created: 02/20/2018
#------------------------------------------------------------------------------

# load packages
library(tidyverse)

# load data
fish_data <- read_csv('./data/sar_data.csv')
spread_data <- read_csv('./data/spread_factors.csv')


# transform fish data into lists 
fish_list <- fish_data %>%
  rename(period = cohort) %>%
  group_by(species, origin) %>%
  nest()
  

# transform cohort spread factors into lists 
spread_list <- spread_data %>%
  group_by(species, origin, spill_schedule) %>%
  nest() %>%
  group_by(species, origin) %>%
  nest()

save(fish_list, spread_list, file = './data/data.rda')
  
