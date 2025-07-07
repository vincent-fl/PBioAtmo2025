library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)

source("Code/MDS_standalone.R")

##### HELPER FUNCTIONS #####
# already defined in HelperFunctions.R in main project
calc_SHFmean <- function(SHF1, SHF2, SHF3) {
  G <- mapply(function(x, y, z) mean(c(x, y, z), na.rm = TRUE),
              SHF1, SHF2, SHF3)
  return(G)
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

# already defined in HelperFunctions.R in main project
load_meteo_data <- function(path) {
  df <- read.csv(path, header=TRUE)
  df$date_time <- as.POSIXct(df$date_time, format = "%Y-%m-%d %H:%M")
  df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df)], as.numeric)
  df$SHF_mean <- calc_SHFmean(df$SHF_1, df$SHF_2, df$SHF_3)
  return(df) 
}

# already defined in HelperFunctions.R in main project
lowNr_to_NA <- function(df_col, threshold=-9999) {
  replace(df_col, df_col <= threshold, NA)
}

combine_data_gapFilling <- function(ec_df, met_df, y_col, x_cols=NULL,
                                    met_mean_cols=NULL, met_sum_cols=NULL) {
  
  # Aggregate meteorological data from 10-minute to 30-minute intervals
  met_df <- tenMin_to_30minData(met_df, mean_cols = met_mean_cols, sum_cols = met_sum_cols)
  
  # Identify overlapping column names (excluding the key)
  common_cols <- intersect(setdiff(names(ec_df), "date_time"), names(met_df))
  # Drop those columns from ec_df so met_df values are preferred
  ec_df <- ec_df[, !names(ec_df) %in% common_cols]
  
  # Merge ecosystem and meteorological data by 'date_time' (full outer join)
  combined_df <- merge(ec_df, met_df, by = "date_time", all = TRUE)
  
  # Decide which columns to select
  if (is.null(x_cols)) {
    cols <- c('date_time', y_col, 'Ta', 'Rn', 'VPD', 'RH', 'Pa', 'SWin', 'SWout', 'LWin', 
              'LWout', 'SHF_mean')
  } else {
    cols <- c('date_time', y_col, x_cols)
  }
  
  # Check if all requested columns exist in the merged data frame
  missing_cols <- setdiff(cols, colnames(combined_df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are missing in the combined data frame:", 
               paste(missing_cols, collapse = ", ")))
  }
  
  # Subset combined data to only the selected columns
  new_df <- combined_df[, cols, drop = FALSE]
  
  return(new_df)
}


##### PLOT FUNCTIONS #####
## 1.1. Plot-Function to display the gap distribution 

# INPUT data for all functions is the Output from the splitTrainTest_byGapsize()
# function (list with 2 dataframes)

# This function plots the entire gap --> all data that is missing 
plot_gap_distribution_by_halfhour <- function(gap_result, gap_size, time_col = "date_time") {
  test_df <- gap_result$test_df
  
  # time in 30-min steps
  test_df$time_of_day <- format(test_df[[time_col]], "%H:%M")
  test_df$time_of_day <- factor(test_df$time_of_day, levels = format(seq(
    as.POSIXct("00:00", format = "%H:%M"),
    as.POSIXct("23:30", format = "%H:%M"),
    by = "30 mins"
  ), "%H:%M"))
  
  ggplot(test_df, aes(x = time_of_day)) +
    geom_bar(fill = "steelblue") +
    scale_x_discrete(
      breaks = c("00:00", "03:00", "06:00", "09:00", "12:00", "15:00", "18:00", "21:00")
    ) +
    labs(
      title = "Frequency of Gaps",
      subtitle = paste("Gap Size:", gap_size),
      x = "Time",
      y = "Number of gaps"
    ) +
    theme_minimal(base_size = 14) 
}

# This function only plots the frequency of the start time of the gap 
plot_gap_start_times <- function(gap_result, gap_size, expected_step_mins = 30, time_col = "date_time") {
  train_df <- gap_result$train_df
  time_vec <- as.POSIXct(train_df[[time_col]])
  
  # calculate time diff
  time_diff <- diff(time_vec)
  expected_diff <- as.difftime(expected_step_mins, units = "mins")
  gap_start_indices <- which(time_diff > expected_diff)
  gap_start_times <- time_vec[gap_start_indices + 1]
  
  # time in 30 moin steps
  all_times <- format(seq(
    as.POSIXct("00:00", format = "%H:%M"),
    as.POSIXct("23:30", format = "%H:%M"),
    by = "30 mins"
  ), "%H:%M")
  
  gap_start_df <- data.frame(start_time = format(gap_start_times, "%H:%M")) |>
    mutate(start_time = factor(start_time, levels = all_times)) |>
    count(start_time) |>
    complete(start_time = factor(all_times, levels = all_times), fill = list(n = 0))
  
  # Plot
  ggplot(gap_start_df, aes(x = start_time, y = n)) +
    geom_col(fill = "grey20", width = 0.9) +
    labs(
      title = "Frequency of gap start times",
      subtitle = paste("Gap Size:", gap_size),
      x = "Time",
      y = "Number of gaps"
    ) +
    scale_x_discrete(
      breaks = all_times[seq(1, length(all_times), by = 6)]  
    ) +
    theme_minimal(base_size = 14) 
}

# This function combines the ones above: It plots the frequency of all gaps and 
# marks the startin time of the gap
plot_gap_dist_combined <- function(gap_result, gap_size, 
                                   expected_step_mins = 30, 
                                   time_col = "date_time") {
  
  test_df <- gap_result$test_df
  train_df <- gap_result$train_df
  time_vec <- as.POSIXct(train_df[[time_col]])
  
  # 30-min steps as time factor
  all_times <- format(seq(
    as.POSIXct("00:00", format = "%H:%M"),
    as.POSIXct("23:30", format = "%H:%M"),
    by = "30 mins"), "%H:%M")
  
  # gaps time of day
  test_df$time_of_day <- factor(format(test_df[[time_col]], "%H:%M"), levels = all_times)
  
  # extract time of gap start
  time_diff <- diff(time_vec)
  expected_diff <- as.difftime(expected_step_mins, units = "mins")
  gap_start_indices <- which(time_diff > expected_diff)
  gap_start_times <- time_vec[gap_start_indices + 1]
  
  gap_start_df <- data.frame(time_of_day = factor(format(gap_start_times, "%H:%M"), levels = all_times))
  
  # count
  gap_counts <- test_df %>% count(time_of_day) %>% mutate(type = "All gaps")
  start_counts <- gap_start_df %>% count(time_of_day) %>% mutate(type = "Gap starts")
  
  # combine
  combined_counts <- bind_rows(gap_counts, start_counts)
  
  # set x axis ticks (hier alle 3 Stunden)
  x_ticks <- all_times[seq(1, length(all_times), by = 6)]
  
  # Plot 
  plot_list <- list()  # Leere Liste erstellen
  
  for (i in 1:5) {
    p <- ggplot(combined_counts, aes(x = time_of_day, y = n, fill = type)) +
      geom_col(position = "identity", width = 0.9, alpha = 0.7, color = "black") +
      scale_fill_manual(name = " ",
                        values = c("All gaps" = "grey85", "Gap starts" = "grey30"),
                        labels = c("All gaps" = "Artificial Gaps", "Gap starts" = "Start of gap"))+
      scale_x_discrete(breaks = x_ticks) +
      labs(
        title = "Gap Distribution and Start Times",
        subtitle = paste("Gap Size:", gap_size),
        x = "Time",
        y = "Count") +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.position = "top")
    
    plot_list[[i]] <- p
  }
  
  return(plot_list)
}


# 1.2. Plot function R^2 and RMSE per gap sizeAdd commentMore actions

## Plot metric as boxplot per daytime or per gap size
boxplot_metric <- function(results_df, x = c("gap_size", "time_of_day"), y, title) {
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  
  x <- match.arg(x)
  
  # Gap size als Faktor ordnen
  if (x == "gap_size"){
    ordered_gap_levels <- c("30 min", "6 h", "12 h", "1 d", "5 d")
    results_df$gap_size <- factor(results_df$gap_size, levels = ordered_gap_levels, ordered = TRUE)
  }
  
  if (!x %in% names(results_df)) {
    stop(paste("Column", x, "not found in data."))
  }
  
  # Filtern nach Methode und Metrik
  df_filtered <- results_df %>%
    select(method, gap_size, repetition, metric_value = all_of(y), x_var = all_of(x))
  
  # Falls kein Ergebnis übrig bleibt: Warnung ausgeben
  if (nrow(df_filtered) == 0) {
    stop(paste("No data for metric =", y))
  }
  
  # Formatierung der Labels
  x_label <- if(x == "gap_size") "gap size" else "Time of day"
  y_label <- if (y == "R2") expression(R^2) else paste("RMSE")
  title_label <- title
  
  # Plot 
  p <- ggplot(df_filtered, aes(x = x_var, y = metric_value, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.8), notch=TRUE) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = if (x == "time_of_day") element_text(angle = 20, hjust = 1) else element_text(),
      panel.border = element_rect(colour = "grey", fill = NA, size = 0.7)
    ) +
    theme(legend.position = "bottom") +
    scale_fill_brewer(palette = "Set2",
                      labels = c("glm" = "GLM",
                                 "lasso" = "Lasso",
                                 "lasso_weighted" = "Lasso weighted",
                                 "mds" = "MDS",
                                 "xgboost" = "XGBoost")) +
    labs(
      title = parse(text = title_label),
      x = x_label,
      y = y_label,
      fill = NULL
      #fill = "Method"
    )
  
  if (y == "R2") {
    p <- p + scale_y_continuous(limits = c(0, 1))
  }
  
  # # Kleiner Inset-Plot mit Mittelwertlinien
  # k <- df_filtered %>%
  #   group_by(x_var, method) %>%
  #   summarise(mean_val = mean(metric_value), .groups = "drop") %>%
  #   ggplot(aes(x = x_var, y = mean_val, color = method, group = method)) +
  #   geom_line() +
  #   # geom_point() +
  #   theme_minimal(base_size = 10) +
  #   labs(x = x_label,
  #        y = y_label,
  #        fill = "Method") +
  #   theme(legend.position = "none",
  #         axis.title = element_text(size = 8))+
  #   scale_color_brewer(palette="Set2")+
  #   if (y == "R2") scale_y_continuous(limits = c(0, 1)) else NULL
  # 
  # # Kombination: Inset oben rechts platzieren
  # final_plot <- ggdraw() +
  #   draw_plot(p, y = 0, height = 1) + 
  #   draw_plot(k, x = 0.7, y = 0.67, width = 0.3, height = 0.25)  # oben rechts
  # 
  return(p)  # return(final_plot)
  
}


##### STATISTIC FUNCTIONS #####
# Installiere das Paket, wenn es noch nicht installiert ist

# Lade die Bibliothek
library(xgboost)

xgboost_fill <- function(df_train, y) {
  y_train <- df_train[[y]]
  x_train <- df_train[, setdiff(names(df_train), y), drop = FALSE]
  
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x_train), label = y_train)
  
  # Trainiere das XGBoost-Modell
  params <- list(
    objective = "reg:squarederror", 
    booster = "gbtree",
    eta = 0.1,
    max_depth = 6
  )
  
  model <- xgboost::xgb.train(params, dtrain, nrounds = 100)
  
  return(model)
}

linear_lasso_standard <- function(df_train, y) {
  y_train <- df_train[[y]]
  x_vars <- setdiff(names(df_train), y)
  
  x_train <- df_train[, x_vars, drop = FALSE]
  
  y_train <- data.matrix(y_train)
  x_train <- data.matrix(x_train)
  
  # Perform k-fold cross-validation to find optimal lambda value
  model <- glmnet::cv.glmnet(x_train, y_train, alpha=1)
  # Find optimal lambda value that minimizes test MSE
  best_lambda <- model$lambda.min
  
  # Final model
  best_model <- glmnet::glmnet(x_train, y_train, alpha=1, lambda=best_lambda)
  
  return(list(best_model=best_model, best_lambda=best_lambda))
}

linear_GLS <- function(df_train, y) {
  x_vars <- setdiff(names(df_train), y)
  
  formula_str <- paste(y, "~", paste(x_vars, collapse = " + "))
  model_formula <- as.formula(formula_str)
  
  model <- nlme::gls(model_formula, data = df_train)
  
  return(model)
}

linear_lasso_weighted <- function(df_train, y) {
  y_train <- df_train[[y]]
  x_vars <- setdiff(names(df_train), y)
  
  x_train <- df_train[, x_vars, drop = FALSE]
  
  y_train <- data.matrix(y_train)
  x_train <- data.matrix(x_train)
  
  # Fit initial OLS model to estimate residual variance
  ols_model <- lm(y_train ~ ., data = as.data.frame(x_train))
  resid_sq <- residuals(ols_model)^2
  
  # Avoid extreme weights: floor at 10th percentile or small epsilon
  threshold <- max(quantile(resid_sq, 0.1), 1e-6)
  resid_sq <- pmax(resid_sq, threshold)
  
  # Compute weights as inverse of variance
  weights <- 1 / resid_sq
  sqrt_weights <- sqrt(weights)
  
  # Transform predictors and response using weights
  x_train_w <- x_train * sqrt_weights
  y_train_w <- y_train * sqrt_weights
  
  # Cross-validation on weighted data
  model_cv <- glmnet::cv.glmnet(x_train_w, y_train_w, alpha = 1)
  best_lambda <- model_cv$lambda.min
  
  # Fit final weighted Lasso model
  best_model <- glmnet::glmnet(x_train_w, y_train_w, alpha = 1, lambda = best_lambda)
  
  return(list(best_model=best_model, best_lambda=best_lambda))
}


##### WORK FLOW FUNCTIONS #####
calculate_n_rows <- function(window_size, time_per_row = '30min') {
  
  # Helper function to extract numeric value and unit from a string like "30min"
  extract_quantity_unit <- function(x) {
    # Remove space characters
    x <- gsub(" ", "", x)
    # Use regular expression to separate number and unit
    matches <- regmatches(x, regexec("^(\\d+(\\.\\d+)?)([a-zA-Z]+)$", x))[[1]]
    
    # If the pattern doesn't match properly (i.e., no number or unit), stop with error
    if (length(matches) < 4) stop("Invalid format: ", x)
    
    # Convert number to integer and unit to lowercase string
    quantity <- as.numeric(matches[2])
    unit <- tolower(matches[4])
    
    return(list(quantity = quantity, unit = unit))
  }
  
  # Helper function to convert a quantity with unit to minutes
  to_minutes <- function(quantity, unit) {
    # Define valid units grouped by their time base
    minute_units <- c('min', 'minute', 'm', 'minutes', 'mins', 'ms', 'minuten')
    hour_units   <- c('h', 'hour', 'stunde', 'hr', 'hs', 'hours', 'stunden', 'hrs')
    day_units    <- c('d', 'day', 'tag', 'ds', 'days', 'tage')
    
    # Convert to minutes based on unit category
    if (unit %in% minute_units) {
      return(quantity)
    } else if (unit %in% hour_units) {
      return(quantity * 60)
    } else if (unit %in% day_units) {
      return(quantity * 1440)  # 24 * 60
    } else {
      stop("Unknown unit: ", unit)
    }
  }
  
  # Parse and convert both window size and time per row to minutes
  win <- extract_quantity_unit(window_size)
  row <- extract_quantity_unit(time_per_row)
  
  qt_win_min <- to_minutes(win$quantity, win$unit)
  qt_row_min <- to_minutes(row$quantity, row$unit)
  
  # Calculate how many rows fit into the window
  n_rows <- round(qt_win_min / qt_row_min)  # maybe floor()?
  
  # Return result if valid, otherwise stop with error
  if (n_rows >= 1) {
    return(n_rows)
  } else {
    stop("time_per_row cannot be larger than window_size!")
  }
}



#' Create Artificial Gaps in Time Series Data
#'
#' @description
#' This function splits a time series dataframe into a training and test set by 
#' removing artificial gaps of a specified size and percentage.
#' The gaps can be defined either by number of rows (if evenly spaced) or 
#' based on actual time duration using a POSIXct timestamp column.
#'
#' @param df            A data frame containing time series data.
#' @param gap_size      A string indicating the gap duration (e.g. "3h", "90min", "2days").
#'                      Accepts units like "min", "m", "minutes", "h", "hr", "hours", "d", "day", etc.
#' @param gap_percents  Percentage of the data to be removed as gaps (e.g., 30 for 30%).
#' @param day_time      Optional vector of two strings (e.g. c("06:00", "18:00")) to restrict gaps
#'                      to certain times of the day.
#' @param time_col      Name of the column containing POSIXct timestamps (default "date_time").
#' @param time_per_row  Estimated spacing between rows if not using time_col (default "30min").
#' @param use_time_col  Boolean; if TRUE, uses actual time stamps for gap sizing.
#'
#' @return A list with:
#'   - train_df: the training data excluding gaps
#'   - test_df: the data in the gaps
#'   - gap_rows: row indices used for the test set

splitTrainTest_byGapsize <- function(df, gap_size, gap_percents = 30, 
                                     day_time = NULL, time_col = "date_time",
                                     time_per_row = '30min', use_time_col = TRUE) {
  df <- df[complete.cases(df), ]
  
  # Parser for time strings like '2h', '30min', '1day'
  parse_gap_size <- function(gap_size) {
    # Remove space characters
    gap_size <- gsub(" ", "", gap_size)
    matches <- regmatches(gap_size, regexec("^(\\d+(\\.\\d+)?)([a-zA-Z]+)$", gap_size))[[1]]
    if (length(matches) < 4) stop("Invalid gap_size format.")
    quantity <- as.numeric(matches[2])
    unit <- tolower(matches[4])
    
    unit <- switch(unit,
                   min = "minutes", m = "minutes", minute = "minutes", minutes = "minutes", mins = "minutes", ms = "minutes", minuten = "minutes",
                   h = "hours", hr = "hours", hrs = "hours", hour = "hours", hours = "hours", stunde = "hours", stunden = "hours", hs = "hours",
                   d = "days", day = "days", days = "days", tag = "days", tage = "days", ds = "days",
                   stop("Unknown or unsupported time unit: ", unit)
    )
    
    return(lubridate::duration(quantity, units = unit))
  }
  
  gap_duration <- parse_gap_size(gap_size)
  total_rows <- nrow(df)
  total_gap_rows_target <- floor(total_rows * (gap_percents / 100))
  
  all_gap_rows <- integer(0)
  gap_starts <- integer(0)
  candidate_indices <- 1:(total_rows - 1)
  
  if (!is.null(day_time)) {
    start_time <- as.POSIXct(day_time[1], format = "%H:%M", tz = "UTC")
    end_time   <- as.POSIXct(day_time[2], format = "%H:%M", tz = "UTC")
    time_of_day <- format(df[[time_col]], "%H:%M")
    time_of_day <- as.POSIXct(time_of_day, format = "%H:%M", tz = "UTC")
  }
  
  if (use_time_col) {
    if (!inherits(df[[time_col]], "POSIXct")) {
      stop("The time_col must contain POSIXct timestamps.")
    }
    
    while (length(all_gap_rows) < total_gap_rows_target && length(candidate_indices) > 0) {
      i <- sample(candidate_indices, 1)
      start_time_val <- df[[time_col]][i]
      end_time_val <- start_time_val + gap_duration
      
      gap_end <- which(df[[time_col]] >= end_time_val)[1]
      if (is.na(gap_end)) break
      gap_indices <- i:(gap_end - 1)
      
      if (length(gap_indices) < 2) {
        candidate_indices <- setdiff(candidate_indices, i)
        next
      }
      
      # Optional time-of-day filtering
      if (!is.null(day_time)) {
        gap_times <- time_of_day[gap_indices]
        if (start_time <= end_time) {
          if (!all(gap_times >= start_time & gap_times <= end_time)) {
            candidate_indices <- setdiff(candidate_indices, i)
            next
          }
        } else {
          if (!all(gap_times >= start_time | gap_times <= end_time)) {
            candidate_indices <- setdiff(candidate_indices, i)
            next
          }
        }
      }
      
      gap_starts <- c(gap_starts, i)
      all_gap_rows <- c(all_gap_rows, gap_indices)
      forbidden <- seq(max(1, i - 1), min(total_rows, gap_end + 1))
      candidate_indices <- setdiff(candidate_indices, forbidden)
    }
  }
  
  # Fallback to row-based estimation
  if (!use_time_col || length(all_gap_rows) == 0) {
    rows_per_gap <- calculate_n_rows(window_size = gap_size, time_per_row = time_per_row)
    n_gaps <- floor(total_rows * (gap_percents / 100) / rows_per_gap)
    if (n_gaps == 0) stop("0 gaps possible with chosen gap size and data frame.")
    
    candidate_indices <- 2:(total_rows - rows_per_gap)
    gap_starts <- c()
    all_gap_rows <- integer(0)
    
    if (!is.null(day_time)) {
      start_time <- as.POSIXct(day_time[1], format = "%H:%M", tz = "UTC")
      end_time   <- as.POSIXct(day_time[2], format = "%H:%M", tz = "UTC")
      time_of_day <- format(df[[time_col]], "%H:%M")
      time_of_day <- as.POSIXct(time_of_day, format = "%H:%M", tz = "UTC")
      
      valid_indices <- c()
      for (i in candidate_indices) {
        gap_indices <- i:(i + rows_per_gap - 1)
        if (length(gap_indices) < rows_per_gap) next
        gap_times <- time_of_day[gap_indices]
        
        if (start_time <= end_time) {
          if (all(gap_times >= start_time & gap_times <= end_time)) {
            valid_indices <- c(valid_indices, i)
          }
        } else {
          if (all(gap_times >= start_time | gap_times <= end_time)) {
            valid_indices <- c(valid_indices, i)
          }
        }
      }
      candidate_indices <- valid_indices
    }
    
    while (length(gap_starts) < n_gaps && length(candidate_indices) > 0) {
      i <- sample(candidate_indices, 1)
      gap_starts <- c(gap_starts, i)
      gap_rows <- i:(i + rows_per_gap - 1)
      all_gap_rows <- c(all_gap_rows, gap_rows)
      forbidden <- seq(max(1, i - 1), min(total_rows, i + rows_per_gap + 1))
      candidate_indices <- setdiff(candidate_indices, forbidden)
    }
  }
  
  df_train <- df[-all_gap_rows, ]
  df_test <- df[all_gap_rows, ]
  
  return(list(train_df = df_train, test_df = df_test, gap_rows = all_gap_rows))
}

# splitTrainTest_byGapsize <- function(df, gap_size, gap_percents = 30, 
#                                      day_time = NULL, time_col = "date_time") {
#   # Function to create artificial gaps in a time series dataframe by removing blocks of rows.
#   #
#   # INPUT:
#   #   df           : DataFrame containing the time series data (assumed ordered by time)
#   #   gap_size     : String specifying gap length with units (e.g. "30hours", "12h")
#   #   gap_percents : Percentage of the dataframe rows to be removed as gaps (default 30%)
#   #   day_time     : Optional vector of two strings (e.g. c("06:00", "18:00")) to restrict gaps to daytime hours
#   #   time_col     : Name of the column that contains POSIXct timestamps (default "time")
#   #
#   # OUTPUT:
#   #   A list with two dataframes:
#   #     train_df : DataFrame excluding rows within the gaps
#   #     test_df  : DataFrame containing only rows within the gaps
#   
#   # Remove rows with any missing values (complete cases)
#   df <- df[complete.cases(df), ]
#   
#   # Calculate the number of rows that correspond to the gap size
#   # This function returns a rounded result
#   rows_per_gap <- calculate_n_rows(window_size = gap_size, time_per_row = '30min')
#   
#   # Calculate how many gaps to insert based on the gap percentage
#   n_gaps <- floor(nrow(df) * (gap_percents/100) / rows_per_gap)
#   if (n_gaps == 0) stop("0 gaps possible with chosen gap size and data frame.")
#   
#   # Candidate indices for gap start:
#   # Exclude first row and last 'rows_per_gap' rows to avoid invalid gaps
#   candidate_indices <- 2:(nrow(df) - rows_per_gap)
#   
#   # If day_time is provided, filter candidates based on time-of-day
#   if (!is.null(day_time)) {
#     start_time <- as.POSIXct(day_time[1], format = "%H:%M", tz = "UTC")
#     end_time   <- as.POSIXct(day_time[2], format = "%H:%M", tz = "UTC")
#     
#     # Extract time-of-day for all rows as POSIXct (ignoring the date)
#     time_of_day <- format(df[[time_col]], format = "%H:%M")
#     time_of_day <- as.POSIXct(time_of_day, format = "%H:%M", tz = "UTC")
#     
#     # Prepare matrix of time values for each candidate gap (rows_per_gap wide)
#     valid_indices <- c()
#     for (i in candidate_indices) {
#       gap_indices <- i:(i + rows_per_gap - 1)
#       gap_times <- time_of_day[gap_indices]
#       
#       if (length(gap_times) < rows_per_gap) next  # skip if gap overflows
#       
#       # Check each time in the gap
#       if (start_time <= end_time) {
#         # Standard window
#         if (all(gap_times >= start_time & gap_times <= end_time)) {
#           valid_indices <- c(valid_indices, i)
#         }
#       } else {
#         # Window wraps over midnight
#         if (all(gap_times >= start_time | gap_times <= end_time)) {
#           valid_indices <- c(valid_indices, i)
#         }
#       }
#     }
#     
#     candidate_indices <- valid_indices
#     
#     if (length(candidate_indices) == 0) {
#       stop("No valid gaps fit fully within the given time window and gap size.")
#     }
#   }
#   
#   # Vector to store chosen gap start indices
#   gap_starts <- c()
#   
#   # Randomly select gap starts while respecting minimum distance between gaps
#   while (length(gap_starts) < n_gaps && length(candidate_indices) > 0) {
#     i <- sample(candidate_indices, 1)
#     gap_starts <- c(gap_starts, i)
#     
#     # Exclude indices around the current gap start to avoid overlapping gaps
#     forbidden <- seq(max(1, i - 1), min(nrow(df), i + rows_per_gap))
#     candidate_indices <- setdiff(candidate_indices, forbidden)
#   }
#   
#   # Collect all row indices belonging to the gaps
#   all_gap_rows <- c()
#   for (start in gap_starts) {
#     gap_rows <- start:(start + rows_per_gap-1)
#     all_gap_rows <- c(all_gap_rows, gap_rows)
#   }
#   
#   # Create train and test datasets by excluding and including gap rows, respectively
#   df_train <- df[-all_gap_rows, ]
#   df_test <- df[all_gap_rows, ]
#   
#   return(list(train_df = df_train, test_df = df_test, gap_rows = all_gap_rows))
# }


#' Perform Gap-Filling Evaluation Using Multiple Methods
#'
#' This function evaluates the performance of various gap-filling methods (e.g., LASSO, GLS, MDS,XGBoost)
#' across different artificial gap durations and multiple repetitions. It introduces gaps into 
#' the dataset, trains or applies the specified models, and evaluates their predictive performance.
#'
#' @param df Data frame containing the input data.
#' @param y Character string; name of the target (response) variable.
#' @param methods Character vector of method names to evaluate (must match keys in `method_dict`).
#' @param gap_sizes Character vector of gap durations (e.g., c("30 mins", "1 hour")).
#' @param repetitions Integer; number of repetitions per gap size (default = 10).
#' @param gap_percents Numeric; percentage of data to remove (as gaps) in each test set (default = 30).
#' @param day_time Optional vector; used to restrict gapping to certain times of day (default = NULL).
#' @param time_col Character string; name of the datetime column in `df` (default = "date_time").
#' @param use_time_col Boolean; if TRUE, time_col will be used for gap_size index (more robust).
#'        If FALSE, a fixed number of rows will be calculated based on `time_per_row`, assuming
#'        regular time intervals between observations.
#' @param set_seed Optional; either a single integer to set the random seed for reproducibility,
#' or a named list of `.Random.seed` objects (e.g., from a previous run) to reproduce exact randomness.
#' If NULL (default), seeds will still be stored but random each time.
#' @param MDS_V1, MDS_V2, MDS_V3 Character strings; optional predictor variable names used by MDS-type methods.
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A data frame with evaluation metrics (R², RMSE) for each method, repetition, and gap size.}
#'   \item{histograms}{A list of ggplot2 objects showing the gap distributions used in each repetition.}
#'   \item{seeds}{A list of random seeds used in each repetition, for reproducibility.}
#' }
#'
#' @examples
#' \dontrun{
#' results <- main_gapfill(df = my_data, y = "NEE", methods = c("lasso", "GLS", "MDS"),
#'                         gap_sizes = c("30 mins", "1 hour"), repetitions = 5)
#' }
main_gapfill <- function(df, y, methods, gap_sizes, repetitions = 10,
                         gap_percents = 30, day_time = NULL,
                         time_col = "date_time", use_time_col = TRUE, set_seed = NULL,
                         MDS_V1 = "SWin", MDS_V2 = "VPD", MDS_V3 = "Ta") {
  
  # Validate input dataframe
  if (!is.data.frame(df)) stop("Input 'df' must be a data frame.")
  if (!y %in% names(df)) stop(paste("Response variable", y, "not found in data."))
  if (!time_col %in% names(df)) stop(paste("Time column", time_col, "not found in data."))
  if (!is.character(gap_sizes)) {
    stop("'gap_sizes' must be a character vector, e.g., c('30 mins', '1 hour').")}
  if (!is.numeric(gap_percents) || gap_percents <= 0 || gap_percents > 100) stop("'gap_percents' must be between 0 and 100.")
  if (!is.null(set_seed)) {
    if (is.numeric(set_seed) && length(set_seed) == 1 && set_seed %% 1 == 0) {
      use_seed_type <- "integer"
    } else if (is.list(set_seed)) {
      use_seed_type <- "list"
    } else {
      stop("'set_seed' must be either NULL, a single integer, or a list of .Random.seed values.")
    }
  } else {
    use_seed_type <- "none"
  }
  if (!all(sapply(list(MDS_V1, MDS_V2, MDS_V3), is.character))) {
    stop("MDS_V1, MDS_V2, and MDS_V3 must all be character strings.")
  }
  
  df <- na.omit(df)  # Remove rows with NA
  
  # Mapping of method names to model functions
  method_dict <- list(
    lasso          = linear_lasso_standard,
    lasso_weighted = linear_lasso_weighted,
    GLS            = linear_GLS,
    MDS            = MDSGapFill,
    xgboost        = xgboost_fill 
  )
  
  # Check that all methods are defined
  unknown_methods <- setdiff(methods, names(method_dict))
  if (length(unknown_methods) > 0) {
    stop(paste("Unknown method(s):", paste(unknown_methods, collapse = ", ")))
  }
  
  # Initialize outputs
  results_list <- list()
  histograms   <- list()
  seeds_list   <- list()
  
  for (gap_size in gap_sizes) {
    for (i in seq_len(repetitions)) {
      
      # Construct name used to identify seed for this run
      seed_name <- paste0("gap_", gap_size, "_rep_", i)
      
      # Seed handling
      if (use_seed_type == "integer") {
        # Use fixed numeric seed across all repetitions
        set.seed(set_seed)
      } else if (use_seed_type == "list") {
        # Use specific pre-stored .Random.seed object if available
        if (!is.null(set_seed[[seed_name]])) {
          assign(".Random.seed", set_seed[[seed_name]], envir = .GlobalEnv)
        } else {
          warning(paste("Seed for", seed_name, "not found in provided seed list."))
          if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
        }
      } else {
        # use_seed_type == "none" → No seeding, use current RNG
        if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
      }
      
      # Store the seed used in this run for reproducibility
      seeds_list[[seed_name]] <- .Random.seed
      
      # Create gapped train/test datasets
      res <- tryCatch({
        splitTrainTest_byGapsize(
          df = df,
          gap_size = gap_size,
          gap_percents = gap_percents,
          day_time = day_time,
          time_col = time_col
        )
      }, error = function(e) {
        warning(sprintf("Skipping gap_size = %s, repetition = %s due to error in gap creation: %s",
                        gap_size, i, e$message))
        return(NULL)
      })
      
      if (is.null(res)) next  # Skip this repetition on failure
      
      # Add NA values to the target column at the gap rows
      df_with_gaps <- df
      df_with_gaps[res$gap_rows, y] <- NA
      
      # Save gap histogram
      plot_obj <- tryCatch({
        plot_gap_distribution_by_halfhour(res, gap_size)
      }, error = function(e) {
        warning("Failed to generate gap histogram.")
        NULL
      })
      histograms[[length(histograms) + 1]] <- plot_obj
      
      train_df <- res$train_df
      test_df  <- res$test_df
      
      # Drop time column
      train_data <- train_df[ , !(names(train_df) %in% time_col)]
      test_data  <- test_df[ , !(names(test_df) %in% time_col)]
      
      x_test <- test_data[ , !(names(test_data) %in% y)]
      y_true <- test_data[[y]]
      
      for (method_name in methods) {
        method_func <- method_dict[[method_name]]
        
        # Try to train model or get direct predictions
        model_or_pred <- tryCatch({
          if (method_name == "MDS") {
            method_func(df_with_gaps, y, V1 = MDS_V1, V2 = MDS_V2, V3 = MDS_V3)
          } else if (method_name == "xgboost") {
            method_func(train_data, y)
          } else {
            method_func(train_data, y)
          }
        }, error = function(e) {
          warning(sprintf("Error in method '%s' at gap size %s, rep %s: %s", method_name, gap_size, i, e$message))
          return(NULL)
        })
        
        if (is.null(model_or_pred)) next
        
        # Handle prediction
        y_pred <- tryCatch({
          if (method_name == "MDS") {
            model_or_pred[res$gap_rows, paste0(y, "_f")]
          } else if (method_name == "xgboost") {
            x_test_dm <- xgboost::xgb.DMatrix(data = as.matrix(x_test))
            predict(model_or_pred, newdata = x_test_dm)
          } else if (inherits(model_or_pred, "list") && inherits(model_or_pred$best_model, "glmnet")) {
            predict(model_or_pred$best_model, newx = as.matrix(x_test), s = model_or_pred$best_lambda) |> as.vector()
          } else if (inherits(model_or_pred, "gls")) {
            predict(model_or_pred, newdata = x_test)
          } else {
            stop("Unknown model type – cannot predict.")
          }
        }, error = function(e) {
          warning(sprintf("Prediction failed for method '%s': %s", method_name, e$message))
          return(NULL)
        })
        
        if (is.null(y_pred)) next
        
        # Performance metrics
        compute_metrics <- function(y_true, y_pred) {
          if (var(y_true, na.rm = TRUE) == 0) return(c(R2 = NA, RMSE = NA))
          r2   <- 1 - sum((y_true - y_pred)^2, na.rm = TRUE) / sum((y_true - mean(y_true, na.rm = TRUE))^2, na.rm = TRUE)
          rmse <- sqrt(mean((y_true - y_pred)^2, na.rm = TRUE))
          return(c(R2 = r2, RMSE = rmse))
        }
        
        metrics <- compute_metrics(y_true, y_pred)
        r2_val <- metrics["R2"]
        rmse_val <- metrics["RMSE"]
        
        # Store result
        results_list[[length(results_list) + 1]] <- data.frame(
          method     = method_name,
          gap_size   = gap_size,
          repetition = i,
          R2         = r2_val,
          RMSE       = rmse_val,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine and return
  result_df <- if (length(results_list) > 0) {
    data.table::rbindlist(results_list)
  } else {
    warning("No successful runs to summarize.")
    data.frame()
  }
  
  return(list(
    results    = result_df,
    histograms = histograms,
    seeds      = seeds_list
  ))
}
