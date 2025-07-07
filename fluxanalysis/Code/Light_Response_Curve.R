# Elisas Skript zum Ausprobieren

# 3. Plot light response curve ---------------------------------------------------
# CO2 flux in response to SWin (for the time before and after mowing)
# T1: before mowing 
# T2: after Mowing 
# pick only fluxed during NEP prime time :) --> 11-15 h 
# convert co2-fluxes (NEE) to NEP      NEP = -NEE

# load functions
source("Code/HelperFunctions.R")
source("Code/PlotFunctions.R")

library(dplyr)
library(ggplot2)
library(minpack.lm)

# load in all the data 
T1_ec <- load_ec_data("Data/eddypro_T1_full_output_2025-06-19T034008_adv.csv")
T2_ec <- load_ec_data("Data/eddypro_T2_full_output_2025-06-18T165605_adv.csv")
T3_ec <- load_ec_data("Data/eddypro_T3_full_output_2025-06-18T003008_adv.csv")
meteo <- load_meteo_data("Data/MeteoFull_PBioAtmo25_16Jun.csv")

# after mowing is T2 and T3
T2_T3_ec = rbind(T2_ec, T3_ec)

# aggregate 10 minute meteo data to 30 minute data 
meteo_30 <- tenMin_to_30minData(meteo)

hours <- c("11:00", "11:30", "12:00", "12:30", "13:00", "13:30", "14:00", "14:30", "15:00")

# data before mowing
T1 <- semi_join(meteo_30, T1_ec, by = "date_time") %>% 
  select(date_time, SWin) %>% 
  mutate(co2_flux = T1_ec$co2_flux, 
         qc_co2_flux = T1_ec$qc_co2_flux, 
         NEP = - co2_flux, # add colum with NEP (-NEE, also hier -co2_flux)
         time = T1_ec$time) %>% 
  filter(time %in% hours,
         qc_co2_flux < 7) # remove qualityflags >= 7

# T1 but not filtered by hours but by SWin > 20
T1_test <- semi_join(meteo_30, T1_ec, by = "date_time") %>% 
  select(date_time, SWin) %>% 
  mutate(co2_flux = T1_ec$co2_flux, 
         qc_co2_flux = T1_ec$qc_co2_flux, 
         NEP = - co2_flux, # add colum with NEP (-NEE, also hier -co2_flux)
         time = T1_ec$time) %>% 
  filter(SWin > 20,
         qc_co2_flux < 7)


# data after mowing 
T2 <- semi_join(meteo_30, T2_T3_ec, by = "date_time") %>% 
  select(date_time, SWin) %>% 
  mutate(co2_flux = T2_T3_ec$co2_flux, 
         qc_co2_flux = T2_T3_ec$qc_co2_flux, 
         NEP = - co2_flux,
         time = T2_T3_ec$time) %>% 
  filter(time %in% hours) %>% 
  filter(qc_co2_flux < 7) # remove qualityflags >= 7


T2_test <- semi_join(meteo_30, T2_ec, by = "date_time") %>% 
  select(date_time, SWin) %>% 
  mutate(co2_flux = T2_ec$co2_flux, 
         qc_co2_flux = T2_ec$qc_co2_flux, 
         NEP = - co2_flux, # add colum with NEP (-NEE, also hier -co2_flux)
         time = T2_ec$time) %>% 
  filter(SWin > 20,
         qc_co2_flux < 7)

# create scatterplot with regression line to look at data distribution

ggplot(T2, aes(x = SWin, y = NEP)) + # bzw T2 einsetzen
  geom_point() +
  geom_smooth()

# fit light response curve with Levenberg-Marquardt Algorithm 

# LRC-Modellformel
lrc_model <- NEP ~ (alpha * SWin * A_max) / (alpha * SWin + A_max) + R_d

# Nichtlineare Regression mit Levenberg–Marquardt
fit <- nlsLM(lrc_model,
             data = T2_test,
             start = list(alpha = 0.03, A_max = 15, R_d = -2))

# Zusammenfassung der Ergebnisse
summary(fit)

#### Plot
# Vorhersage auf feinem Raster
SWin_pred <- seq(0, 2000, length.out = 200)
NEP_pred <- predict(fit, newdata = data.frame(SWin = SWin_pred))

# Erstelle DataFrame mit Vorhersagen
pred_df <- data.frame(SWin = SWin_pred, NEP = NEP_pred)


# Plot mit korrekter Zuordnung
ggplot(T2, aes(x = SWin, y = NEP)) +
  geom_point(color = "blue", alpha = 0.5) +
  geom_line(data = pred_df, aes(x = SWin, y = NEP), color = "red", linewidth = 1) +
  labs(title = "Light Response Curve (LRC) after mowing ",
       x = "SWin (µmol m⁻² s⁻¹)",
       y = "NEP (µmol CO₂ m⁻² s⁻¹)") +
  theme_minimal()

# Residuen
rss <- sum((T1$NEP - NEP_pred)^2)              # Residual Sum of Squares
tss <- sum((T1$NEP - mean(T1$NEP))^2)   # Total Sum of Squares

r_squared <- 1 - (rss / tss)
r_squared



