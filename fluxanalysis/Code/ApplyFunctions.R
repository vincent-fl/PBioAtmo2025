# Beispiel-Skript zur Anwendung von Vincents Funktionen
# Stand 01.06.2025
# zuletzt bearbeitet von Elisa zu Erstellung der Plots für CO2
# H von Flo

# load functions
source("Code/HelperFunctions.R")
source("Code/PlotFunctions.R")

# load full output data
full_out_2605 = load_ec_data("Data/eddypro_advanced_01_full_output_2025-05-26T212140_adv.csv")
full_out_2805 = load_ec_data("Data/eddypro_advanced_01_full_output_2025-05-28T182940_adv.csv")
meteo_full = load_meteo_data("Data/MeteoFull_PBioAtmo25_30May.csv")

# Available arguments and default values:
# y_col, QC_col (only _qc plots), x_col, max_val = NULL, min_val = NULL, save = FALSE, show_NA = FALSE

# 1.1 CO2 -----------------------------------------------------------------------

# create non-interactive ggplots (simple ant with qc values)
noninteractive_ggplot_simple(full_out_2605$co2_flux, full_out_2605$date_time)
noninteractive_ggplot_qc(full_out_2605$co2_flux,full_out_2605$qc_co2_flux, 
                         full_out_2605$date_time)

noninteractive_ggplot_simple(full_out_2805$co2_flux, full_out_2805$date_time)
noninteractive_ggplot_qc(full_out_2805$co2_flux,full_out_2805$qc_co2_flux, 
                         full_out_2805$date_time)

#create interactive dygraphs (simple and with qc values)
interactive_dygraph_simple(full_out_2605$co2_flux, full_out_2605$date_time)
interactive_dygraph_qc(full_out_2605$co2_flux, full_out_2605$qc_co2_flux,
                       full_out_2605$date_time)

interactive_dygraph_simple(full_out_2805$co2_flux, full_out_2805$date_time)
interactive_dygraph_qc(full_out_2805$co2_flux, full_out_2805$qc_co2_flux, 
                           full_out_2805$date_time)

# 1.2 H -----------------------------------------------------------------------
# create non-interactive ggplots (simple ant with qc values)
noninteractive_ggplot_simple(full_out_2605$H, full_out_2605$date_time)
noninteractive_ggplot_qc(full_out_2605$H,full_out_2605$qc_H, 
                         full_out_2605$date_time)

noninteractive_ggplot_simple(full_out_2805$H, full_out_2805$date_time)
noninteractive_ggplot_qc(full_out_2805$H,full_out_2805$qc_H, 
                         full_out_2805$date_time)


#create interactive dygraphs (simple and with qc values)
interactive_dygraph_simple(full_out_2605$H, full_out_2605$date_time)
interactive_dygraph_qc(full_out_2605$H, full_out_2605$qc_H,
                       full_out_2605$date_time)

interactive_dygraph_simple(full_out_2805$H, full_out_2805$date_time)
interactive_dygraph_qc(full_out_2805$H, full_out_2805$qc_H, 
                       full_out_2805$date_time)

# 1.3 LE ------------------------------------------------------------------------
# create non-interactive ggplots (simple ant with qc values)
noninteractive_ggplot_simple(full_out_2605$LE, full_out_2605$date_time)
noninteractive_ggplot_qc(full_out_2605$LE,full_out_2605$qc_LE, 
                         full_out_2605$date_time, save=T)

noninteractive_ggplot_simple(full_out_2805$LE, full_out_2805$date_time)
noninteractive_ggplot_qc(full_out_2805$LE,full_out_2805$qc_LE, 
                         full_out_2805$date_time, save=T)


#create interactive dygraphs (simple and with qc values)
interactive_dygraph_simple(full_out_2605$LE, full_out_2605$date_time, show_NA=T)
interactive_dygraph_qc(full_out_2605$LE, full_out_2605$qc_LE,
                       full_out_2605$date_time, show_NA=T, save=T)

interactive_dygraph_simple(full_out_2805$LE, full_out_2805$date_time, show_NA=T)
interactive_dygraph_qc(full_out_2805$LE, full_out_2805$qc_LE, 
                       full_out_2805$date_time, show_NA=T, save=T)


# 2. Meteo Data
# 2.1. Radiative Balance -------------------------------------------------------

ggplot_radiative_balance(meteo_full, x_col = date_time,
                       y_col1 = SWin,
                       y_col2 = SWout,
                       y_col3 = LWin,
                       y_col4 = LWout,
                       min_date = "2025-05-25", max_date = "2025-05-28")

dygraph_radiative_balance(meteo_full, x_col = date_time,
                               y_col1 = SWin,
                               y_col2 = SWout,
                               y_col3 = LWin,
                               y_col4 = LWout)


#=======
# 3) Energy balance ------------------------------------------------------------
plot_energybalance_ggplot(ec_df=full_out_2605, meteo_df=meteo_data, save=T)
plot_energybalance_ggplot(ec_df=full_out_2805, meteo_df=meteo_data, save=T)


# 2.2. Energy balance ----------------------------------------------------------

plot_energybalance_ggplot(ec_df=full_out_2605, meteo_df=meteo_full, save=T,
                          max_qc=7)
plot_energybalance_ggplot(ec_df=full_out_2805, meteo_df=meteo_full, save=T,
                          max_qc=7)

#=======
# 3) poster plots --------------------------------------------------------------


