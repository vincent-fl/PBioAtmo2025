
data <- read.csv("Data/ec_meteo_combined_gapfilled.csv")
data$date_time <- as.POSIXct(data$date_time)
test <- data[,c("date_time", "co2_flux_f")]
test$cumsum_co2 <- cumsum(test$co2_flux_f)
plot(x=test$date_time, y=test$cumsum_co2, type="l")
test[nrow(test),]
#date_time co2_flux_f cumsum_co2
#1825 2025-06-16    8.84515  -3752.078

# from μmol m^-2 to g m^-2 (time = our measure period)
-3752.078 * 1e-6 * 44.01
#-0.165129
