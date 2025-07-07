source("Code/helperfunctions.R")

library(dplyr)
library(lubridate)

T1 <- load_ec_data("Data/eddypro_T1_full_output_2025-06-19T034008_adv.csv")
T2 <- load_ec_data("Data/eddypro_T2_full_output_2025-06-18T165605_adv.csv")
T3 <- load_ec_data("Data/eddypro_T3_full_output_2025-06-18T003008_adv.csv")

met <- load_meteo_data("Data/MeteoFull_PBioAtmo25_16Jun.csv")

start_T1 <- as.POSIXct("2025-05-09 00:00:00")
end_T1 <- as.POSIXct("2025-05-15 16:00:00")

start_T2 <- as.POSIXct("2025-05-15 17:00:00")
end_T2 <- as.POSIXct("2025-05-19 13:00:00")

start_T3 <- as.POSIXct("2025-05-19 14:00:00")
end_T3 <- as.POSIXct("2025-06-16 00:00:00")


start_met <- as.POSIXct("2025-05-09 00:00:00")
end_met <- as.POSIXct("2025-06-16 00:00:00")

T1 <- T1[T1$date_time >= start_T1 & T1$date_time <= end_T1, ]
T2 <- T2[T2$date_time >= start_T2 & T2$date_time <= end_T2, ]
T3 <- T3[T3$date_time >= start_T3 & T3$date_time <= end_T3, ]

T_all <- bind_rows(T1, T2, T3)
full_time <- data.frame(
  date_time = seq(from = min(start_T1),
                  to = max(end_T3),
                  by = "30 min")
)
T_all <- full_time %>%
  left_join(T_all, by = "date_time")
names(T_all)[names(T_all) == "RH"] <- "RH_EddyPro"

met <- met[met$date_time >= start_met & met$date_time <= end_met, ]

# time stemp at end of period: ceiling_date(); at beginning: floor_date()
met <- met %>%
  mutate(date_time = floor_date(date_time, unit = "30 minutes"))

met <- met %>%
  group_by(date_time) %>%
  summarise(
    Ta       = mean(Ta, na.rm = TRUE),
    RH       = mean(RH, na.rm = TRUE),
    Pa       = mean(Pa, na.rm = TRUE),
    SWin     = mean(SWin, na.rm = TRUE),
    SWout    = mean(SWout, na.rm = TRUE),
    LWin     = mean(LWin, na.rm = TRUE),
    LWout    = mean(LWout, na.rm = TRUE),
    SHF_1    = mean(SHF_1, na.rm = TRUE),
    SHF_2    = mean(SHF_2, na.rm = TRUE),
    SHF_3    = mean(SHF_3, na.rm = TRUE),
    Rn       = mean(Rn, na.rm = TRUE),
    albedo   = mean(albedo, na.rm = TRUE),
    P        = sum(P, na.rm = TRUE),
    SHF_mean = mean(SHF_mean, na.rm = TRUE)
  )

combined_df <- merge(T_all, met, by = "date_time", all = TRUE)
combined_df$date_time <- format(combined_df$date_time, "%Y-%m-%d %H:%M:%S")

write.csv(combined_df, "Data/ec_meteo_combined.csv", row.names = F)
