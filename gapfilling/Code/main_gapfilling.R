library(patchwork)
source("Code/Functions_gapfilling.R")

# ============================================================================ #
##### DATA #####

# load full output data
EC <- read.csv("Data/ec_meteo_combined.csv")
EC$date_time <- as.POSIXct(EC$date_time, format="%Y-%m-%d %H:%M:%S")


###### CO2 ######
df_co2 <- EC[EC$qc_co2_flux < 7, ]
df_co2 <- df_co2[, c('date_time', 'co2_flux', 'Ta', 'SWin', 'VPD', 'RH', 'Pa', 'SHF_mean')]
# , 'SWin', 'SWout', 'LWin', 'LWout', , 'VPD', 'Rn'

df_co2 <- lowNr_to_NA(df_co2)
df_co2 <- na.omit(df_co2)

df_co2 <- df_co2[df_co2$co2_flux < 20 & df_co2$co2_flux > -25, ]  # Unnecessary?


###### H2O ######
df_h2o <- EC[EC$qc_h2o_flux < 7, ]
df_h2o <- df_h2o[, c('date_time', 'h2o_flux', 'Ta', 'SWin', 'VPD', 'RH', 'Pa', 'SHF_mean')]

df_h2o <- lowNr_to_NA(df_h2o)
df_h2o <- na.omit(df_h2o)


# ============================================================================ #
##### WORKFLOW #####

###### gapsize wise ######

gapsize_wise <- function(df, y, set_seed) {
  main_gapfill(
    df = df,
    y = y,
    methods = c("lasso", "lasso_weighted", "GLS", "MDS", "xgboost"),
    gap_sizes = c("30 min", "6 h", "12 h", "1 d", "5 d"),
    repetitions = 100,
    gap_percents = 20,
    day_time = NULL,
    time_col = "date_time",
    use_time_col = T,
    set_seed = set_seed,  # if NULL random seed still gets stored in the result list
    MDS_V1="SWin", MDS_V2="VPD", MDS_V3="Ta"
  )
}

###### CO2 ######
reslist_co2 <- gapsize_wise(df_co2, "co2_flux", set_seed = NULL)
model_res_co2 <- reslist_co2$results  # View(model_res_co2)
hists_co2 <- reslist_co2$histograms  # plot with e.g. hists[[1]]
seeds <- reslist_co2$seeds

saveRDS(model_res_co2, file = "Results/model_res_co2.rds")
saveRDS(hists_co2, file = "Results/hists_co2.rds")
saveRDS(seeds, file = "Results/seeds.rds")

###### H2O ######
reslist_h2o <- gapsize_wise(df_h2o, "h2o_flux", set_seed = seeds)
model_res_h2o <- reslist_h2o$results  # View(model_res)
hists_h2o <- reslist_h2o$histograms  # plot with e.g. hists[[1]]

saveRDS(model_res_h2o, file = "Results/model_res_h2o.rds")
saveRDS(hists_h2o, file = "Results/hists_h2o.rds")

###### by day time ######

daytime_wise <- function(df, y, day_time, set_seed=seeds) {
  main_gapfill(
    df = df,
    y = y,
    methods = c("lasso", "lasso_weighted", "GLS", "MDS", "xgboost"),
    gap_sizes = c("30min"),
    repetitions = 30,
    gap_percents = 5,
    day_time = day_time,
    time_col = "date_time",
    use_time_col = T,
    set_seed = set_seed,
    MDS_V1="SWin", MDS_V2="VPD", MDS_V3="Ta"
  )
}

run_daytime_gapfill <- function(df, y, set_seed = seeds) {
  # Define the time slots
  timeslots <- list(
    c("00:00", "02:30"),
    c("03:00", "05:30"),
    c("06:00", "08:30"),
    c("09:00", "11:30"),
    c("12:00", "14:30"),
    c("15:00", "17:30"),
    c("18:00", "20:30"),
    c("21:00", "23:30")
  )
  
  # Initialize empty list to collect results
  result_list <- list()
  hists_list <- list()
  
  # Loop through each time slot
  for (slot in timeslots) {
    time_label <- paste(slot[1], slot[2], sep = "-")
    
    result <- daytime_wise(df = df, y = y, day_time = slot, set_seed = set_seed)
    
    res_df <- result$results
    res_df$time_of_day <- time_label
    
    result_list[[time_label]] <- res_df
    
    hists_list[[time_label]] <- result$histograms
  }
  
  # Combine all into one data frame
  final_results <- do.call(rbind, result_list)
  
  return(list(results = final_results, histograms = hists_list))
}

###### CO2 ######
reslist_co2_2 <- run_daytime_gapfill(df=df_co2, y="co2_flux", set_seed=seeds)
slot_results_co2 <- reslist_co2_2$results
slot_hists_co2 <- reslist_co2_2$histograms
saveRDS(slot_results_co2, file = "Results/slot_results_co2.rds")
saveRDS(slot_hists_co2, file = "Results/slot_hists_co2.rds")

###### H2O ######
reslist_h2o_2 <- run_daytime_gapfill(df=df_h2o, y="h2o_flux", set_seed=seeds)
slot_results_h2o <- reslist_h2o_2$results
slot_hists_h2o <- reslist_h2o_2$histograms
saveRDS(slot_results_h2o, file = "Results/slot_results_h2o.rds")
saveRDS(slot_hists_h2o, file = "Results/slot_hists_h2o.rds")

# ============================================================================ #
# readRDS(c("Results/hists_co2.rds", "Results/hists_h2o.rds", 
#           "Results/model_res_co2.rds", "Results/model_res_h2o.rds", 
#           "Results/slot_results_co2.rds", "Results/slot_results_h2o.rds",
#           "Results/seeds.rds"))
#seeds <- readRDS("Results/seeds.rds")
model_res_co2 <- readRDS("Results/model_res_co2.rds")
model_res_h2o <- readRDS("Results/model_res_h2o.rds")
slot_results_co2 <- readRDS("Results/slot_results_co2.rds")
slot_results_h2o <- readRDS("Results/slot_results_h2o.rds")

hists1 <- readRDS("Results/hists_co2.rds")
hists1[[1]]

##### PLOTS #####

###### CO2 ######
p1 <- boxplot_metric(model_res_co2, x="gap_size", y="RMSE", title="CO[2]")
ggsave("Plots/co2_gapsize_rmse.png", plot = p1, width = 8, height = 3, units = "in", dpi = 500)
p1
p2 <- boxplot_metric(model_res_co2, x="gap_size", y="R2", title="CO[2]")
ggsave("Plots/co2_gapsize_r2.png", plot = p2, width = 8, height = 3, units = "in", dpi = 500)
p2

p3 <- boxplot_metric(slot_results_co2, x="time_of_day", y="RMSE", title="CO[2]")
ggsave("Plots/co2_daytime_rmse.png", plot = p3, width = 8, height = 3, units = "in", dpi = 500)
p3
p4 <- boxplot_metric(slot_results_co2, x="time_of_day", y="R2", title="CO[2]")
ggsave("Plots/co2_daytime_r2.png", plot = p4, width = 8, height = 3, units = "in", dpi = 500)
p4

###### H2O ######
p5 <- boxplot_metric(model_res_h2o, x="gap_size", y="RMSE", title="H[2]*O")
ggsave("Plots/h2o_gapsize_rmse.png", plot = p5, width = 8, height = 3, units = "in", dpi = 500)
p5
p6 <- boxplot_metric(model_res_h2o, x="gap_size", y="R2", title="H[2]*O")
ggsave("Plots/h2o_gapsize_r2.png", plot = p6, width = 8, height = 3, units = "in", dpi = 500)
p6

p7 <- boxplot_metric(slot_results_h2o, x="time_of_day", y="RMSE", title="H[2]*O")
ggsave("Plots/h2o_daytime_rmse.png", plot = p7, width = 8, height = 3, units = "in", dpi = 500)
p7
p8 <- boxplot_metric(slot_results_h2o, x="time_of_day", y="R2", title="H[2]*O")
ggsave("Plots/h2o_daytime_r2.png", plot = p8, width = 8, height = 3, units = "in", dpi = 500)
p8

###### combined ######
# pA
p1_1 <- p1 + labs(title = NULL)
p2_1 <- p2 + labs(x = NULL) +
  theme(axis.text.x = element_blank())
p5_1 <- p5 + labs(title = NULL, y = NULL)
p6_1 <- p6 + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())

pA <- ((p2_1 / p1_1) | (p6_1 / p5_1)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("Plots/r2_rmse_combined_gapsize_co2_h2o.png", plot = pA, width = 8, height = 6, units = "in", dpi = 500)
pA

# pB
p3_1 <- p3 + labs(title = NULL) #+
  #theme(axis.text.x = element_text(size = 6))
p4_1 <- p4 + labs(x = NULL) +
  theme(axis.text.x = element_blank())
p7_1 <- p7 + labs(title = NULL, y = NULL) #+
  #theme(axis.text.x = element_text(size = 6))
p8_1 <- p8 + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank())

pB <- ((p4_1 / p3_1) | (p8_1 / p7_1)) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave("Plots/r2_rmse_combined_daytime_co2_h2o.png", plot = pB, width = 8, height = 6, units = "in", dpi = 500)
pB
