library(tidyverse)
library(sbtools)
library(zip)
library(ncdf4)

# get helper functions
source('code/file_functions.R')

# download SB forecasts
temp <- tempfile()
sbtools::item_file_download(sb_id = '5f6a28a782ce38aaa2449137', 
                            names = 'forecast[2021-04-16_2021-07-16]_files.zip', destinations = temp)

zip::unzip(zipfile = temp, exdir = 'in/forecasts')
unlink(temp)

# NOTE # 
# I didn't take the time to figure out automated recursive unzip. 
# The commented out code below works for unzipping multiple files 
# but takes much longer than navigating to these zipped files
# select all > right click > 7zip > Extract all here (on Windows)
# such that all files are unzipped to in/forecasts

# unzip many zip files
#file_names <- list.files("in/forecasts", full.names = TRUE)
#walk(file_names, ~ unzip(zipfile = .x, 
#                         exdir = 'in/forecasts'))

# download SB temperature data
temp <- tempfile()
sbtools::item_file_download(sb_id = '5f6a287382ce38aaa2449131', 
                            names = 'temperature_observations_forecast_sites.zip', destinations = temp)

zip::unzip(zipfile = temp, exdir = 'in')
unlink(temp)


# get list of files you want to process
forecast_files <- return_da_files('in/forecasts', model = 'DA')

# extract predictions from each file and bind together
forecasts <- extract_and_bind(forecast_files)

# pair with observation data
obs <- readr::read_csv('in/temperature_observations_forecast_sites.csv') %>%
  mutate(seg_id_nat = as.character(seg_id_nat)) %>%
  select(seg_id_nat, time = date, obs_max_temp_c = max_temp_c)

sites <- tibble(site_name = c('DR @ Lordville', 'EBDR @ Harvard','WBDR @ Hancock','WBDR @ Hale Eddy','NR @ Bridgeville'),
                seg_id_nat = c('1573', '1450', '1571', '1565', '1641'))

# now obs are ready to merge with forecasts
# may want to do this after you summarize the ensembles
# since this creates repeated obs per ensemble
# but including this code to show obs are ready to merge
forecasts_with_obs <- forecasts %>%
  left_join(obs) %>%
  left_join(sites)


