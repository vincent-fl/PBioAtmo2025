# load functions
source("Code/HelperFunctions.R")
source("Code/PlotFunctions.R")

library(dplyr)
library(lubridate)
library(dygraphs)

# load full output data
full_out_2605 = load_ec_data("Data/eddypro_advanced_01_full_output_2025-05-26T212140_adv.csv")
full_out_2805 = load_ec_data("Data/eddypro_advanced_01_full_output_2025-05-28T182940_adv.csv")
meteo_full = load_meteo_data("Data/MeteoFull_PBioAtmo25_30May.csv")

# combine datasets
full_data <- bind_rows(full_out_2605, full_out_2805)

# transform date time
full_data$date_time <- ymd_hms(full_data$date_time)

# filter all values above -999
full_data <- full_data %>% 
  filter(co2_flux > -999)

# hourly
co2_hourly <- full_data %>%
  mutate(date_time = floor_date(date_time, "hour")) %>%
  group_by(date_time) %>%
  summarise(co2_sum = sum(co2_flux, na.rm = TRUE))
  #mutate(co2_cum = cumsum(co2_sum))

# daily
co2_daily <- full_data %>%
  mutate(date_time = floor_date(date_time, "day")) %>%
  group_by(date_time) %>%
  summarise(co2_sum = sum(co2_flux, na.rm = TRUE))
  #mutate(co2_cum = cumsum(co2_sum))

# weekly
co2_weekly <- full_data %>%
  mutate(date_time = floor_date(date_time, "week")) %>%
  group_by(date_time) %>%
  summarise(co2_sum = sum(co2_flux, na.rm = TRUE))
  #mutate(co2_cum = cumsum(co2_sum))

# cleaning NA-values
co2_hourly <- co2_hourly %>%
  filter(complete.cases(.))
co2_daily <- co2_daily %>%
  filter(complete.cases(.))
co2_weekly <- co2_weekly %>%
  filter(complete.cases(.))

# create interactive dygraphs
co2_sum_hourly <- xts::xts(co2_hourly$co2_sum, order.by = co2_hourly$date_time)
co2_sum_daily <- xts::xts(co2_daily$co2_sum, order.by = co2_daily$date_time)
co2_sum_weekly <- xts::xts(co2_weekly$co2_sum, order.by = co2_weekly$date_time)

interactive_dygraph_simple(co2_sum_hourly, co2_hourly$date_time)
interactive_dygraph_simple(co2_sum_daily, co2_daily$date_time)
interactive_dygraph_simple(co2_sum_weekly, co2_weekly$date_time)

#######
# combine with meteo data
# convert date time
meteo_full$date_time <- ymd_hms(meteo_full$date_time)

# cleaning data: NA-, NaN-, Inf-values 
meteo_full <- meteo_full %>%
  filter(complete.cases(.)) %>%
  filter(is.finite(Ta) & is.finite(RH) & is.finite(Pa) & is.finite(SWin))

# daily sum
meteo_daily <- meteo_full %>%
  mutate(date_time = floor_date(date_time, "day")) %>%
  group_by(date_time) %>%
  summarise(Ta = mean(Ta, na.rm = TRUE),
            RH = mean(RH, na.rm = TRUE),
            Pa = mean(Pa, na.rm = TRUE),
            SWin = mean(SWin, na.rm = TRUE),
            SWout = mean(SWout, na.rm = TRUE),
            LWin = mean(LWin, na.rm = TRUE),
            LWout = mean(LWout, na.rm = TRUE),
            SHF_mean = mean(SHF_mean, na.rm = TRUE)) %>%
  ungroup()

# combine CO2 data with meteo data
co2_meteo_daily <- left_join(co2_daily, meteo_daily, by = c("date_time" = "date_time"))

# create interactive dygraphs
co2_meteo_daily <- xts::xts(co2_meteo_daily, order.by = co2_meteo_daily$date_time)
dygraph(co2_meteo_daily, main = "Aggregated CO₂-Flux (sum) and Meteo Data") 

