getwd()
library(lattice)
library(fitdistrplus)
library(MASS)
library(survival)
library(RColorBrewer)
library(lattice)
library(permute)
library(vegan)
library(nlme)
library(splines)
library(carData)
library(effects)
library(car)
library(mgcv)
library(nlme)
library(MuMIn)
library(Matrix)
library(lme4)
library(MASS)
library(effects)
library(multcomp)
library(mvtnorm)
library(TH.data)
library(FSA)
library(rcompanion)
library(ggplot2)
library(ggsignif)
library(reshape)
library(ggpubr)
library(GGally)
library(psych)
library(tidyverse)
library(lmtest)
library(lsr)
library(rstatix)
library(dplyr)
library(Hmisc)
library(DHARMa)
library(brant)
library(emmeans)
library(dplyr)
library(lubridate)
library(dygraphs)
library(patchwork)


EC_full <- read.csv("ec_data_fullfull.csv", header=TRUE, sep= ",", dec = ".")
meteo_full <- read.csv("MeteoFull_PBioAtmo25_16Jun.csv", header=TRUE, sep= ",", dec = ".")


# load functions
source("Code/HelperFunctions.R")
source("Code/PlotFunctions.R")

####Plot CO2, H2O, Ta, VPD & SWin####

##2.1 transform date time & quality filter EC data

#solving midnight problem of missing time in EC file
EC_full$date_time <- ifelse(
  grepl("^\\d{4}-\\d{2}-\\d{2}$", EC_full$date_time),
  paste0(EC_full$date_time, " 00:00:00"),
  EC_full$date_time)


# transform date time
EC_full$date_time <- ymd_hms(EC_full$date_time)


# filter all values above -999 (1824 -> 1814 rows, 0.33 % loss)
EC_full <- EC_full %>% 
  filter(co2_flux >- 999, h2o_flux > -999)


#filter all quality classes above 6 - gq = good quality (1814 -> 1497 rows, 17.47 % loss)
EC_full_gq <- EC_full %>%
  filter(qc_co2_flux < 7, qc_h2o_flux < 7)



##2.2 transform date time & quality filter meteo data

# transform date time (without seconds!)
meteo_full$date_time <- ymd_hm(meteo_full$date_time)


#cleaning data: NA-, NaN-, Inf-values (5695 -> 5179 rows, 9.06 % loss)
meteo_full_gq <- meteo_full %>%
  filter(complete.cases(.)) %>%
  filter(is.finite(Ta) & is.finite(SWin))


###3. 2 Periods: before & after mowing
##3.1 set mowing periods
T1_start <- ymd_hms("2025-05-09 00:00:00")
T1_end  <- ymd_hms("2025-05-15 16:00:00")
T2_start <- ymd_hms("2025-05-15 17:00:00")
T2_end  <- ymd_hms("2025-06-16 00:00:00")


##3.2 combine daytime with mowing periods
EC_2p <- EC_full_gq %>%
  filter((date_time >= T1_start & date_time <= T1_end) |
           (date_time >= T2_start & date_time <= T2_end)) %>%
  mutate(Period = case_when(
    date_time >= T1_start & date_time <= T1_end ~ "T1",
    date_time >= T2_start & date_time <= T2_end ~ "T2"),
    Time = format(date_time, "%H:%M")) %>%
  filter(!is.na(Period))


meteo_2p <- meteo_full_gq %>%
  filter((date_time >= T1_start & date_time <= T1_end) |
           (date_time >= T2_start & date_time <= T2_end)) %>%
  mutate(Period = case_when(
    date_time >= T1_start & date_time <= T1_end ~ "T1",
    date_time >= T2_start & date_time <= T2_end ~ "T2"),
    Time = format(date_time, "%H:%M")) %>%
  filter(!is.na(Period))



###4. Calculate median & SE per daytime for CO2, H2O, Ta & SWin
EC_median <- EC_2p %>%
  group_by(Period, Time) %>%
  summarise(median_flux_co2 = median(co2_flux, na.rm = TRUE),
            sd_co2 = sd(co2_flux, na.rm = TRUE),
            n_co2 = sum(!is.na(co2_flux)),
            se_co2 = sd_co2 / sqrt(n_co2),
            
            median_flux_h2o = median(h2o_flux, na.rm = TRUE),
            sd_h2o = sd(h2o_flux, na.rm = TRUE),
            n_h2o = sum(!is.na(h2o_flux)),
            se_h2o = sd_h2o / sqrt(n_h2o),
            
            .groups = "drop")



meteo_median <- meteo_2p %>%
  group_by(Period, Time) %>%
  summarise(median_swin = median(SWin, na.rm = TRUE),
            sd_swin = sd(SWin, na.rm = TRUE),
            n_swin = sum(!is.na(SWin)),
            se_swin = sd_swin / sqrt(n_swin),
            
            median_ta = median(Ta, na.rm = TRUE),
            sd_ta = sd(Ta, na.rm = TRUE),
            n_ta = sum(!is.na(Ta)),
            se_ta = sd_ta / sqrt(n_ta),
            
            .groups = "drop")

###5. Plot

##5.1 Plot CO2 Fluxes
#scale y-axis in 2-steps starting with 0, 2, 4 etc.
y_min_co2 <- floor(min(EC_median$median_flux_co2 - EC_median$se_co2, na.rm = TRUE) / 2) * 2 
y_max_co2 <- ceiling(max(EC_median$median_flux_co2 + EC_median$se_co2, na.rm = TRUE) / 2) * 2


#plot, y-axis in 2 steps
CO2_plot <- ggplot(EC_median, aes(x = Time, y = median_flux_co2, color = Period, 
                                  group = Period)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median_flux_co2 - se_co2,
                  ymax = median_flux_co2 + se_co2,
                  fill = Period),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  
  annotate("segment", x = "12:00", xend = "12:00",y = 0, yend = 5,
           arrow = arrow(length = unit(0.6, "cm")), color = "black", size = 1.7) +
  
  annotate("text", x = "12:00", y = 5.6, label = "Source", vjust = 0, size = 13,
           fontface = "bold") +
  
  annotate("segment", x = "12:00", xend = "12:00", y = 0, yend = -5,
           arrow = arrow(length = unit(0.6, "cm")), color = "black", size = 1.7) +
 
  annotate("text", x = "12:00", y = -5.6, label = "Sink", vjust = 1, size = 13,
           fontface = "bold") +
  
  scale_x_discrete(breaks = format(seq(ISOdatetime(2000, 1, 1, 0, 0, 0),
                                       by = "2 hour", length.out = 12), "%H:%M"),
                   guide = guide_axis(angle = 45), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(y_min_co2, y_max_co2, by = 2), 
                     expand = c(0, 0)) +
  scale_color_manual(values = c("T1" = "#1f77b4", "T2" = "#ff7f0e"),
                     labels = c("T1" = "before mowing", "T2" = "after mowing")) +
  scale_fill_manual(values = c("T1" = "#1f77b4", "T2" = "#ff7f0e"),
                    labels = c("T1" = "before mowing", "T2" = "after mowing")) +
  labs(title = " ", x = " ", y = " ") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.03, size = 20),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major.x = element_line(color = "gray90"),
        legend.position = "right",
        legend.justification = "center",
        legend.background = element_rect(fill = "#FFFFFFB3", color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(size = 23))

CO2_plot


#plot, y-axis in 4 steps
CO2_plot <- ggplot(EC_median, aes(x = Time, y = median_flux_co2, color = Period, 
                                  group = Period)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median_flux_co2 - se_co2,
                  ymax = median_flux_co2 + se_co2,
                  fill = Period),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  
  annotate("segment", x = "12:00", xend = "12:00",y = 0, yend = 5,
           arrow = arrow(length = unit(0.7, "cm")), color = "black", size = 1.6) +
  
  annotate("text", x = "12:00", y = 5.6, label = "Source", vjust = 0, size = 13,
           fontface = "bold") +
  
  annotate("segment", x = "12:00", xend = "12:00", y = 0, yend = -5,
           arrow = arrow(length = unit(0.7, "cm")), color = "black", size = 1.6) +
  
  annotate("text", x = "12:00", y = -5.6, label = "Sink", vjust = 1, size = 13,
           fontface = "bold") +
  
  scale_x_discrete(breaks = format(seq(ISOdatetime(2000, 1, 1, 0, 0, 0),
                                       by = "2 hour", length.out = 12), "%H:%M"),
                   guide = guide_axis(angle = 45), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(floor(y_min_co2 / 4) * 4, 
                                  ceiling(y_max_co2 / 4) * 4, by = 4), 
                     expand = c(0, 0)) +
  scale_color_manual(values = c("T1" = "#1f77b4", "T2" = "#ff7f0e"),
                     labels = c("T1" = "before mowing", "T2" = "after mowing")) +
  scale_fill_manual(values = c("T1" = "#1f77b4", "T2" = "#ff7f0e"),
                    labels = c("T1" = "before mowing", "T2" = "after mowing")) +
  labs(title = " ", x = " ", y = " ") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.03, size = 20),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.grid.major.x = element_line(color = "gray90"),
        legend.position = "right",
        legend.justification = "center",
        legend.background = element_rect(fill = "#FFFFFFB3", color = "black"),
        legend.title = element_blank(),
        axis.text = element_text(size = 34))

CO2_plot




##5.2 Plot H2O Fluxes
#scale y-axis in 0.5-steps starting with 0, 0.5, 1.0 etc.
y_min_h2o <- floor(min(EC_median$median_flux_h2o - EC_median$se_h2o, 
                       na.rm = TRUE) / 0.5) * 0.5
y_max_h2o <- ceiling(max(EC_median$median_flux_h2o + EC_median$se_h2o, 
                         na.rm = TRUE) / 0.5) * 0.5


#plot
ggplot(EC_median, aes(x = Time, y = median_flux_h2o, color = Period, group = Period)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median_flux_h2o - se_h2o,
                  ymax = median_flux_h2o + se_h2o,
                  fill = Period),
              alpha = 0.2, color = NA) +
  scale_x_discrete(breaks = format(seq(ISOdatetime(2000, 1, 1, 0, 0, 0),
                                       by = "2 hours", length.out = 12), "%H:%M"),
                   guide = guide_axis(angle = 45)) +
  scale_y_continuous(breaks = seq(y_min_h2o, y_max_h2o, by = 0.5), expand = c(0, 0)) +
  labs(title = "Diurnal median H2O-Fluxes before & after mowing (± SE)",
       x = "Time of Day",
       y = expression(paste("H"[2], "O-Flux [", mu, "mol ", m^{-2}, s^{-1}, "]")),
       color = "Period",
       fill = "Period") +
  theme_minimal(base_size = 14)


##5.3 Plot H2O, SWin & Ta
#join H2O with meteo data
meteo_H2O_median <- left_join(EC_median %>% 
                                select(Period, Time, median_flux_h2o, sd_h2o,
                                       n_h2o, se_h2o), meteo_median,
                              by = c("Period", "Time"))


#rename column
meteo_H2O_median <- meteo_H2O_median %>%
  rename(median_h2o = median_flux_h2o)


#long & wide format to plot
meteo_H2O_median_wider <- meteo_H2O_median %>%
  pivot_longer(cols = starts_with("median") | starts_with("se"),
               names_to = c("stat", "variable"),
               names_sep = "_",
               values_to = "value") %>%
  pivot_wider(names_from = stat,
              values_from = value)


#label each headline
headlines <- c(h2o = " ",
               swin = " ",
               ta = " ")


#plot
Rest_plot <- ggplot(meteo_H2O_median_wider, aes(x = Time, y = median, 
                                                group = Period, color = Period, 
                                                fill = Period)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = median - se, ymax = median + se), alpha = 0.4, color = NA) +
  facet_wrap(~ variable, scales = "free_y", ncol = 1,
             labeller = labeller(variable = function(x) "")) +
  scale_x_discrete(breaks = sprintf("%02d:00", seq(0, 23, by = 2)),
                   guide = guide_axis(angle = 45), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("T1" = "#1f77b4", "T2" = "#ff7f0e"),
                     labels = c("T1" = "before mowing", "T2" = "after mowing")) +
  scale_fill_manual(values = c("T1" = "#1f77b4", "T2" = "#ff7f0e"),
                    labels = c("T1" = "before mowing", "T2" = "after mowing")) +
  expand_limits(y = 0) +
  labs(title = " ", x = " ", y = "", color = "Period", fill = "Period") +
  theme_minimal(base_size = 14) +
  theme(strip.text = element_text(size = 20, hjust = 0.03, face = "plain"), 
        legend.position = "right",
        legend.justification = "center",
        legend.background = element_rect(fill = "#FFFFFFB3", color = "black"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing = unit(-1, "lines"),
        plot.margin = margin(t = - 25, r = 10, b = 2, l = 2),
        axis.text = element_text(size = 23))
Rest_plot

