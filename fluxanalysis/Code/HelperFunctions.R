

load_ec_data <- function(path) {
  df <- read.csv(path, skip=1, header=TRUE)
  df <- df[-1, ]
  df$date_time <- as.POSIXct(paste(df$date, df$time), format = "%Y-%m-%d %H:%M")
  df <- df[, c(1, ncol(df), 2:(ncol(df)-1))]  # date_time column to position 2
  df[, 5:ncol(df)] <- lapply(df[, 5:ncol(df)], as.numeric)
  return(df)
}


load_biomet_data <- function(path) {
  df <- read.csv(path, header=TRUE)
  df <- df[-1, ]
  df$date_time <- as.POSIXct(paste(df$date, df$time), format = "%Y-%m-%d %H:%M")
  df <- df[, c(ncol(df), 1:(ncol(df)-1))]  # date_time column to position 2
  df[, 4:ncol(df)] <- lapply(df[, 4:ncol(df)], as.numeric)
  return(df)
}


calc_SHFmean <- function(SHF1, SHF2, SHF3) {
  G <- mapply(function(x, y, z) mean(c(x, y, z), na.rm = TRUE),
              SHF1, SHF2, SHF3)
  return(G)
}


# not for biomet file, but for MeteoFull file
load_meteo_data <- function(path) {
  df <- read.csv(path, header=TRUE)
  df$date_time <- as.POSIXct(df$date_time, format = "%Y-%m-%d %H:%M")
  df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], as.numeric)
  df$SHF_mean <- calc_SHFmean(df$SHF_1, df$SHF_2, df$SHF_3)
  return(df) 
}


group_qc <- function(qc_col) {
  # Map 1-3 -> 0, 4-6 -> 1, 7-9 -> 2
  cut(qc_col,
      breaks = c(-Inf, 3, 6, Inf),
      labels = c(0, 1, 2),
      right = TRUE)
}


filter_qc <- function(df, data_col, qc_col, max_qc = 7) {
  df[[data_col]][df[[qc_col]] > max_qc] <- NA
  return(df)
}


lowNr_to_NA <- function(df_col, threshold=-9999) {
  replace(df_col, df_col <= threshold, NA)
}

tenMin_to_30minData <- function(df, date_time = "date_time", 
                                mean_cols = NULL, sum_cols = NULL) {
  
  # Time rounding helper: floor to 00 or 30 minutes
  round_to_30min <- function(t) {
    mins <- as.numeric(format(t, "%M"))
    secs <- as.numeric(format(t, "%S"))
    t <- t - mins * 60 - secs              # to full hour
    t + ifelse(mins >= 30, 1800, 0)        # add 30 min if needed
  }
  
  # Extract and round time column
  time_vals <- df[[date_time]]
  df$grp <- round_to_30min(time_vals)
  
  # Determine columns to average/sum
  all_cols <- setdiff(names(df), c(date_time, "grp"))
  
  if (is.null(mean_cols)) {
    mean_cols <- setdiff(all_cols, "P")
  }
  
  if (is.null(sum_cols)) {
    sum_cols <- intersect(all_cols, "P")  # only sum "P" if it's there
  }
  
  # Aggregate means
  result_mean <- aggregate(df[mean_cols], by = list(date_time = df$grp), 
                           FUN = function(x) mean(x, na.rm = TRUE))
  
  # Aggregate sums (if any)
  if (length(sum_cols) > 0) {
    result_sum <- aggregate(df[sum_cols], by = list(date_time = df$grp), 
                            FUN = function(x) sum(x, na.rm = TRUE))
    # Merge both
    result <- merge(result_mean, result_sum, by = "date_time")
  } else {
    result <- result_mean
  }
  
  return(result)
}

