source("Code/HelperFunctions.R")


recurrence_plot <- function(time_series, date_time = NULL, dimension = 2, delay = 1, epsilon = NULL, 
                            method = "euclidean", save = FALSE) {
  # Check for required packages
  if (!requireNamespace("tseriesChaos", quietly = TRUE) || 
      !requireNamespace("stats", quietly = TRUE) || 
      !requireNamespace("graphics", quietly = TRUE)) {
    stop("The package '", pkg, "' is required but not installed. Please install it first.")
  }
  
  # Extract name of y_col for plot title and saving
  time_series_name_full <- deparse(substitute(time_series))
  time_series_name <- if (grepl("\\$", time_series_name_full)) {
    sub(".*\\$(.*)", "\\1", time_series_name_full)
  } else {
    time_series_name_full
  }
  
  # Check if date_time is provided and valid
  if (!is.null(date_time)) {
    if (length(date_time) != length(time_series)) {
      stop("Length of 'date_time' must match length of 'time_series'")
    }
    # If date_time is a vector of class Date or POSIXt, convert to character for plotting
    date_labels <- as.character(date_time)
  }
  
  # Phase space reconstruction
  embedded <- tseriesChaos::embedd(time_series, m = dimension, d = delay)
  
  # Distance matrix
  dist_matrix <- as.matrix(stats::dist(embedded, method = method))
  
  # Set epsilon if not given (e.g., 10% quantile of all distances)
  if (is.null(epsilon)) {
    epsilon <- stats::quantile(dist_matrix, 0.1)
  }
  
  # Recurrence matrix
  R <- ifelse(dist_matrix <= epsilon, 1, 0)
  
  # Plot with date_time on axes if provided, else numeric indices
  if (!is.null(date_time)) {
    # Subsample date labels if too many to avoid overcrowding (optional)
    n <- nrow(R)
    max_labels <- 10
    step <- max(1, floor(n / max_labels))
    tick_positions <- seq(1, n, by = step)
    tick_labels <- date_labels[tick_positions]
    
    graphics::image(1:n, 1:n, R,
                    col = c("white", "black"),
                    xaxt = "n", yaxt = "n",
                    xlab = "Time i", ylab = "Time j",
                    main = paste("Recurrence plot of", time_series_name))
    # X-axis
    axis(1, at = tick_positions, labels = FALSE)
    text(x = tick_positions, y = par("usr")[3] - 10.5,  # slightly below the plot
         labels = tick_labels, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
    
    # Y-axis
    axis(2, at = tick_positions, labels = FALSE)
    text(x = par("usr")[1] - 5.5, y = tick_positions,  # slightly left of the plot
         labels = tick_labels, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
  } else {
    graphics::image(1:nrow(R), 1:ncol(R), R,
                    col = c("white", "black"),
                    xlab = "time i", ylab = "time j",
                    main = paste("Recurrence plot of", time_series_name))
  }
  
  if (save) {
    file_name <- paste0("Plots/RecurrencePlot_", time_series_name, "_graphics.png")
    
    png(filename = file_name, width = 800, height = 800)
    if (!is.null(date_time)) {
      graphics::image(1:n, 1:n, R,
                      col = c("white", "black"),
                      xaxt = "n", yaxt = "n",
                      xlab = "Time i", ylab = "Time j",
                      main = paste("Recurrence plot of", time_series_name))
      # X-axis
      axis(1, at = tick_positions, labels = FALSE)
      text(x = tick_positions, y = par("usr")[3] - 0.5,  # slightly below the plot
           labels = tick_labels, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
      
      # Y-axis
      axis(2, at = tick_positions, labels = FALSE)
      text(x = par("usr")[1] - 0.5, y = tick_positions,  # slightly left of the plot
           labels = tick_labels, srt = 45, adj = 1, xpd = TRUE, cex = 0.7)
    } else {
      graphics::image(1:nrow(R), 1:ncol(R), R,
                      col = c("white", "black"),
                      xlab = "time i", ylab = "time j",
                      main = paste("Recurrence plot of", time_series_name))
    }
    dev.off()
  }
  
  return(invisible(R))
}



interactive_dygraph_multi <- function(y_cols, x_col, max_val = NULL, min_val = NULL, 
                                      save = FALSE, show_NA = FALSE, normalize = FALSE) {
  # Check required packages
  required_pkgs <- c("dygraphs", "xts", "zoo", "htmlwidgets")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("The package '", pkg, "' is required but not installed. Please install it first."))
    }
  }
  
  # Check that y_cols is a data.frame, matrix or a single vector
  if (is.vector(y_cols)) {
    y_cols <- data.frame(Series = y_cols)
  } else {
    y_cols <- as.data.frame(y_cols)
  }
  
  # Optional normalization to [0,1]
  if (normalize) {
    y_cols <- as.data.frame(lapply(y_cols, function(col) {
      rng <- range(col, na.rm = TRUE)
      if (diff(rng) == 0) return(rep(0.5, length(col))) # constant vector
      (col - rng[1]) / diff(rng)
    }))
  }
  
  # Apply lowNr_to_NA to each column
  y_cols <- as.data.frame(lapply(y_cols, lowNr_to_NA))
  
  # Create xts object
  xts_val <- xts::xts(y_cols, order.by = x_col)
  
  # Identify timestamps with NA values (optional)
  na_times <- if (show_NA) {
    unique(zoo::index(xts_val)[apply(is.na(xts_val), 1, any)])
  } else {
    NULL
  }
  
  # Title based on column names
  title <- paste("Plot of", paste(names(y_cols), collapse = ", "))
  
  # Create dygraph plot
  dygraph_plot <- dygraphs::dygraph(xts_val, main = title)
  dygraph_plot <- dygraphs::dyRangeSelector(dygraph_plot)
  
  # Add vertical lines for NA timestamps (if any)
  if (!is.null(na_times) && length(na_times) > 0) {
    for (t in na_times) {
      dygraph_plot <- dygraphs::dyEvent(
        dygraph_plot, t,
        label = "", color = "red", strokePattern = "dashed"
      )
    }
  }
  
  # Add max and min limit lines (if provided)
  if (!is.null(max_val)) {
    dygraph_plot <- dygraphs::dyLimit(
      dygraph_plot, max_val,
      label = "max", color = "red", strokePattern = "dashed"
    )
  }
  if (!is.null(min_val)) {
    dygraph_plot <- dygraphs::dyLimit(
      dygraph_plot, min_val,
      label = "min", color = "red", strokePattern = "dashed"
    )
  }
  
  # Save widget as HTML file (if requested)
  if (save) {
    start_date <- format(min(x_col, na.rm = TRUE), "%d%m")
    end_date <- format(max(x_col, na.rm = TRUE), "%d%m")
    file_name <- paste0(
      "Plots/dygraph_", paste(names(y_cols), collapse = "_"),
      "_", start_date, "_", end_date, ".html"
    )
    
    htmlwidgets::saveWidget(
      dygraph_plot,
      file = file_name,
      selfcontained = TRUE
    )
  }
  
  # Return dygraph object
  return(dygraph_plot)
}



interactive_dygraph_simple <- function(y_col, x_col, max_val = NULL, min_val = NULL, save = FALSE, show_NA = FALSE) {
  # Check required packages
  required_pkgs <- c("dygraphs", "xts", "zoo", "htmlwidgets")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("The package '", pkg, "' is required but not installed. Please install it first."))
    }
  }
  
  # Extract name of y_col for plot title and saving
  y_col_name_full <- deparse(substitute(y_col))
  y_col_name <- if (grepl("\\$", y_col_name_full)) {
    sub(".*\\$(.*)", "\\1", y_col_name_full)
  } else {
    y_col_name_full
  }
  
  # Replace low numbers with NA (custom function)
  y_col <- lowNr_to_NA(y_col)
  
  # Create xts object
  xts_val <- xts::xts(y_col, order.by = x_col)
  
  # Identify timestamps with NA values (optional)
  na_times <- if (show_NA) {
    zoo::index(xts_val)[is.na(xts_val)]
  } else {
    NULL
  }
  
  # Create dygraph plot
  dygraph_plot <- dygraphs::dygraph(xts_val, main = y_col_name)
  dygraph_plot <- dygraphs::dyRangeSelector(dygraph_plot)
  
  # Add vertical lines for NA timestamps (if any)
  if (!is.null(na_times) && length(na_times) > 0) {
    for (t in na_times) {
      dygraph_plot <- dygraphs::dyEvent(
        dygraph_plot, t,
        label = "", color = "red", strokePattern = "dashed"
      )
    }
  }
  
  # Add max limit line (if provided)
  if (!is.null(max_val)) {
    dygraph_plot <- dygraphs::dyLimit(
      dygraph_plot, max_val,
      label = "max", color = "red", strokePattern = "dashed"
    )
  }
  
  # Add min limit line (if provided)
  if (!is.null(min_val)) {
    dygraph_plot <- dygraphs::dyLimit(
      dygraph_plot, min_val,
      label = "min", color = "red", strokePattern = "dashed"
    )
  }
  
  # Save widget as HTML file (if requested)
  if (save) {
    # Extract start and end date from x_col
    start_date <- format(min(x_col, na.rm = TRUE), "%d%m")
    end_date <- format(max(x_col, na.rm = TRUE), "%d%m")
    
    # Create filename with y_col_name and date range
    file_name <- paste0("Plots/", y_col_name, "_", start_date, "_", end_date, "_dygraph_simple.html")
    
    htmlwidgets::saveWidget(
      dygraph_plot,
      file = file_name,
      selfcontained = TRUE
    )
  }
  
  # Return dygraph object
  return(dygraph_plot)
}


interactive_dygraph_qc <- function(y_col, QC_col, x_col, max_val = NULL, min_val = NULL, save = FALSE, show_NA = FALSE) {
  # List of required packages
  required_pkgs <- c("dygraphs", "xts", "zoo", "htmlwidgets")
  
  # Check if each package is installed
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("The package '", pkg, "' is required but not installed. Please install it first."))
    }
  }
  
  # Extract column name for plot title and file naming
  y_col_name_full <- deparse(substitute(y_col))
  y_col_name <- if (grepl("\\$", y_col_name_full)) {
    sub(".*\\$(.*)", "\\1", y_col_name_full)
  } else {
    y_col_name_full
  }
  
  # Replace low values with NA
  y_col <- lowNr_to_NA(y_col)
  
  # Create base data frame
  df <- data.frame(
    time = as.POSIXct(x_col),
    value = y_col,
    qc = group_qc(QC_col)
  )
  
  # Extract NA timestamps if requested
  if (show_NA) {
    na_times <- df$time[is.na(df$value)]
  } else {
    na_times <- NULL
  }
  
  # Add a new factor column to split by QC group (0, 1, 2)
  df$qc_level <- factor(df$qc, levels = 0:2, labels = c("QC_0", "QC_1", "QC_2"))
  
  # Add row ID for reshaping
  df$row_id <- seq_len(nrow(df))
  
  # Reshape data to wide format (manually replacing pivot_wider)
  df_QC0 <- df$value
  df_QC0[df$qc_level != "QC_0"] <- NA
  df_QC1 <- df$value
  df_QC1[df$qc_level != "QC_1"] <- NA
  df_QC2 <- df$value
  df_QC2[df$qc_level != "QC_2"] <- NA
  
  # Add "All Data" column
  df$All_Data <- df$value
  
  # Combine final data frame for xts
  df_wide <- data.frame(
    time = df$time,
    All_Data = df$All_Data,
    QC_0 = df_QC0,
    QC_1 = df_QC1,
    QC_2 = df_QC2
  )
  
  # Create xts object
  xts_data <- xts::xts(df_wide[, c("All_Data", "QC_0", "QC_1", "QC_2")], order.by = df_wide$time)
  
  # Create dygraph plot
  dy <- dygraphs::dygraph(xts_data, main = y_col_name)
  dy <- dygraphs::dyOptions(dy, drawPoints = TRUE, pointSize = 2, strokeWidth = 1.2)
  dy <- dygraphs::dySeries(dy, "All_Data", label = "All Data", color = "lightgrey", drawPoints = FALSE)
  dy <- dygraphs::dySeries(dy, "QC_0", label = "QC Group 0 (1–3)", color = "#0072B2")
  dy <- dygraphs::dySeries(dy, "QC_1", label = "QC Group 1 (4–6)", color = "orange")
  dy <- dygraphs::dySeries(dy, "QC_2", label = "QC Group 2 (7–9)", color = "#D55E50")
  dy <- dygraphs::dyAxis(dy, "x", label = "Time")
  dy <- dygraphs::dyAxis(dy, "y", label = y_col_name)
  dy <- dygraphs::dyLegend(dy, show = "always", width = 500)
  dy <- dygraphs::dyHighlight(dy, highlightCircleSize = 4, highlightSeriesBackgroundAlpha = 0.2)
  dy <- dygraphs::dyRangeSelector(dy)
  
  # Add max/min lines if given
  if (!is.null(max_val)) {
    dy <- dygraphs::dyLimit(dy, max_val, color = "red", label = "Max")
  }
  if (!is.null(min_val)) {
    dy <- dygraphs::dyLimit(dy, min_val, color = "red", label = "Min")
  }
  
  # Mark NA values if requested
  if (!is.null(na_times) && length(na_times) > 0) {
    for (na_time in na_times) {
      dy <- dygraphs::dyEvent(dy, na_time, "Missing", labelLoc = "bottom", color = "red")
    }
  }
  
  # Save widget as HTML if requested
  if (save) {
    # Extract start and end date from x_col
    start_date <- format(min(x_col, na.rm = TRUE), "%d%m")
    end_date <- format(max(x_col, na.rm = TRUE), "%d%m")
    
    # Create filename with y_col_name and date range
    file_name <- paste0("Plots/", y_col_name, "_", start_date, "_", end_date, "_dygraph_qc.html")
    
    htmlwidgets::saveWidget(
      dy,
      file = file_name,
      selfcontained = TRUE
    )
  }
  
  return(dy)
}


noninteractive_ggplot_simple <- function(y_col, x_col, max_val = NULL, min_val = NULL, save = FALSE, show_NA = FALSE) {
  # Required package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is required but not installed. Please install it first.")
  }
  
  # Extract column name for labeling
  y_col_name_full <- deparse(substitute(y_col))
  y_col_name <- if (grepl("\\$", y_col_name_full)) {
    sub(".*\\$(.*)", "\\1", y_col_name_full)
  } else {
    y_col_name_full
  }
  
  # Apply custom NA conversion
  y_col <- lowNr_to_NA(y_col)
  
  # Create data frame for plotting
  df <- data.frame(x = x_col, y = y_col)
  
  # Create base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = y_col_name, x = "Time")
  
  # Show NA as vertical red lines
  if (show_NA) {
    na_times <- df$x[is.na(df$y)]
    if (length(na_times) > 0) {
      p <- p + ggplot2::geom_vline(xintercept = na_times, color = "red", linetype = "dashed")
    }
  }
  
  # Add horizontal max and min lines
  if (!is.null(max_val)) {
    p <- p +
      ggplot2::geom_hline(yintercept = max_val, color = "red", linetype = "dashed") +
      ggplot2::annotate("text", x = min(df$x), y = max_val, label = "max", vjust = -0.5, hjust = 0, color = "red")
  }
  if (!is.null(min_val)) {
    p <- p +
      ggplot2::geom_hline(yintercept = min_val, color = "red", linetype = "dashed") +
      ggplot2::annotate("text", x = min(df$x), y = min_val, label = "min", vjust = 1.5, hjust = 0, color = "red")
  }
  
  # Additional theme cleanup
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Save plot if requested
  if (save) {
    # Extract start and end date from x_col
    start_date <- format(min(x_col, na.rm = TRUE), "%d%m")
    end_date <- format(max(x_col, na.rm = TRUE), "%d%m")
    
    # Create filename with y_col_name and date range
    file_name <- paste0("Plots/", y_col_name, "_", start_date, "_", end_date, "_ggplot_simple.png")
    
    ggplot2::ggsave(
      filename = file_name,
      plot = p,
      width = 10,
      height = 6
    )
  }
  
  return(p)
}


noninteractive_ggplot_qc <- function(y_col, QC_col, x_col, max_val = NULL, min_val = NULL, save = FALSE, show_NA = FALSE) {
  # Required package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The package 'ggplot2' is required but not installed. Please install it first.")
  }
  
  # Extract column name for labeling
  y_col_name_full <- deparse(substitute(y_col))
  y_col_name <- if (grepl("\\$", y_col_name_full)) {
    sub(".*\\$(.*)", "\\1", y_col_name_full)
  } else {
    y_col_name_full
  }
  
  # Replace low values with NA
  y_col <- lowNr_to_NA(y_col)
  
  # Apply QC grouping
  qc_group <- group_qc(QC_col)
  
  # Create data frame
  df <- data.frame(x = x_col, y = y_col, qc = qc_group)
  df$qc <- factor(df$qc, levels = c("0", "1", "2"))
  
  # Define colorblind-friendly QC colors
  qc_colors <- c("0" = "#0072B2", "1" = "orange", "2" = "#D55E50")
  
  # Base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_point(ggplot2::aes(color = qc), shape = 4, size = 2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = qc_colors, name = "QC Group") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = y_col_name, x = "Time")
  
  # Add NA vertical lines if requested
  if (show_NA) {
    na_times <- df$x[is.na(df$y)]
    if (length(na_times) > 0) {
      p <- p + ggplot2::geom_vline(xintercept = na_times, color = "red", linetype = "dashed")
    }
  }
  
  # Add max/min limits
  if (!is.null(max_val)) {
    p <- p +
      ggplot2::geom_hline(yintercept = max_val, color = "red", linetype = "dashed") +
      ggplot2::annotate("text", x = min(df$x), y = max_val, label = "max", vjust = -0.5, hjust = 0, color = "red")
  }
  if (!is.null(min_val)) {
    p <- p +
      ggplot2::geom_hline(yintercept = min_val, color = "red", linetype = "dashed") +
      ggplot2::annotate("text", x = min(df$x), y = min_val, label = "min", vjust = 1.5, hjust = 0, color = "red")
  }
  
  # Additional styling
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  # Save if requested
  if (save) {
    # Extract start and end date from x_col
    start_date <- format(min(x_col, na.rm = TRUE), "%d%m")
    end_date <- format(max(x_col, na.rm = TRUE), "%d%m")
    
    # Create filename with y_col_name and date range
    file_name <- paste0("Plots/", y_col_name, "_", start_date, "_", end_date, "_ggplot_qc.png")
    
    ggplot2::ggsave(
      filename = file_name,
      plot = p,
      width = 10,
      height = 6
    )
  }
  
  return(p)
}


ggplot_radiative_balance <- function(data, 
                                     x_col, 
                                     y_col1, y_col2, y_col3, y_col4, 
                                     min_date = NULL, 
                                     max_date = NULL) {
  # Unquoted Spaltennamen erfassen
  x_col <- enquo(x_col)
  y_col1 <- enquo(y_col1)
  y_col2 <- enquo(y_col2)
  y_col3 <- enquo(y_col3)
  y_col4 <- enquo(y_col4)
  
  # Zeitstempel-Spalte extrahieren (z.B. date_time)
  time_values <- pull(data, !!x_col)
  
  # Falls keine min/max gesetzt, nimm gesamte Zeitspanne
  if (is.null(min_date)) min_date <- as.Date(min(time_values, na.rm = TRUE))
  else min_date <- as.Date(min_date)
  
  if (is.null(max_date)) max_date <- as.Date(max(time_values, na.rm = TRUE))
  else max_date <- as.Date(max_date)
  
  # Filtere nach Zeitraum
  data_filtered <- data %>%
    filter(as.Date(!!x_col) >= min_date & as.Date(!!x_col) <= max_date)
  
  # Plot erstellen
  p <- ggplot(data_filtered, aes(x = !!x_col)) +
    geom_line(aes(y = !!y_col1, colour = as_label(y_col1))) +
    geom_line(aes(y = !!y_col2, colour = as_label(y_col2))) +
    geom_line(aes(y = !!y_col3, colour = as_label(y_col3))) +
    geom_line(aes(y = !!y_col4, colour = as_label(y_col4))) +
    ylab(expression(Wm^-2)) +
    xlab("Time") +
    labs(
      colour = "Radiation",
      title = "Radiative balance",
      subtitle = paste0(
        format(min(pull(data_filtered, !!x_col)), "%d.%m.%Y"),
        " to ",
        format(max(pull(data_filtered, !!x_col)), "%d.%m.%Y")
      )
    ) +
    theme_minimal()
  
  return(p)
}

dygraph_radiative_balance <- function(data,
                                      x_col,
                                      y_col1, y_col2, y_col3, y_col4) {
  # Required packages
  required_pkgs <- c("dygraphs", "xts", "dplyr", "rlang")
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("The package '", pkg, "' is required but not installed. Please install it first."))
    }
  }
  
  # Unquoted column names to quosures
  x_col   <- rlang::enquo(x_col)
  y_col1  <- rlang::enquo(y_col1)
  y_col2  <- rlang::enquo(y_col2)
  y_col3  <- rlang::enquo(y_col3)
  y_col4  <- rlang::enquo(y_col4)
  
  # Select and arrange relevant columns
  data_selected <- dplyr::select(
    data,
    time  = !!x_col,
    SWin  = !!y_col1,
    SWout = !!y_col2,
    LWin  = !!y_col3,
    LWout = !!y_col4
  )
  data_selected <- dplyr::arrange(data_selected, time)
  
  # Convert to xts object
  xts_data <- xts::xts(data_selected[, -1], order.by = data_selected$time)
  
  # Create dygraph
  d <- dygraphs::dygraph(xts_data, main = "Radiative Balance")
  d <- dygraphs::dyOptions(d, colors = c("gold", "orange", "blue", "darkblue"))
  d <- dygraphs::dyRangeSelector(d)
  d <- dygraphs::dyAxis(d, "y", label = "W/m²")
  d <- dygraphs::dyLegend(d, show = "follow", hideOnMouseOut = TRUE)
  
  return(d)
}



plot_energybalance_ggplot <- function(ec_df, meteo_df, LE = "LE", H = "H", Rn = "Rn", SHF = "SHF_mean", 
                                      qc_H = "qc_H", qc_LE = "qc_LE", max_qc=7,
                                      align_metData_30minInt = TRUE, save = FALSE) {
  
  # List of required packages
  required_pkgs <- c("ggplot2", "stats", "base", "lubridate")
  
  # Check if each package is installed
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("The package '", pkg, "' is required but not installed. Please install it first."))
    }
  }
  
  ec_df <- filter_qc(df=ec_df, data_col=LE, qc_col=qc_LE, max_qc=max_qc)
  ec_df <- filter_qc(df=ec_df, data_col=H, qc_col=qc_H, max_qc=max_qc)
  
  if (align_metData_30minInt) {
    if (!inherits(meteo_df$date_time, "POSIXct")) {
      stop("'date_time' column in meteo_df must be POSIXct.")
    }
    # Round down all times to the nearest 30-minute interval
    meteo_df$time_30min <- lubridate::floor_date(meteo_df$date_time, unit = "30 minutes")
    
    # Keep only numeric columns for aggregation
    numeric_cols <- sapply(meteo_df, is.numeric)
    meteo_df <- meteo_df[, c("time_30min", names(numeric_cols[numeric_cols]))]
    
    # Aggregate over 30-minute intervals
    meteo_df <- stats::aggregate(. ~ time_30min, data = meteo_df, FUN = mean, na.rm = TRUE)
    
    # Rename for merging
    names(meteo_df)[names(meteo_df) == "time_30min"] <- "date_time"
  }
  
  # Double-check both dataframes have 'date_time'
  if (!"date_time" %in% names(ec_df) || !"date_time" %in% names(meteo_df)) {
    stop("Both data frames must contain a 'date_time' column for merging.")
  }
  
  # Find overlapping time range
  common_start <- max(min(ec_df$date_time), min(meteo_df$date_time))
  common_end   <- min(max(ec_df$date_time), max(meteo_df$date_time))
  
  # Subset to overlapping range
  ec_df <- ec_df[ec_df$date_time >= common_start & ec_df$date_time <= common_end, ]
  meteo_df <- meteo_df[meteo_df$date_time >= common_start & meteo_df$date_time <= common_end, ]
  
  # Merge by timestamp
  merged_df <- base::merge(ec_df, meteo_df, by = "date_time")
  
  # Calculate x and y
  merged_df$x_data <- merged_df[[Rn]] - merged_df[[SHF]]
  merged_df$x_data <- lowNr_to_NA(merged_df$x_data)
  merged_df$y_data <- merged_df[[LE]] + merged_df[[H]]
  merged_df$y_data <- lowNr_to_NA(merged_df$y_data)
  
  # Fit linear model
  fit <- stats::lm(y_data ~ x_data, data = merged_df)
  slope <- round(stats::coef(fit)[2], 2)
  
  # Determine plot limits for annotation placement
  x_range <- range(merged_df$x_data, na.rm = TRUE)
  y_range <- range(merged_df$y_data, na.rm = TRUE)
  x_pos <- x_range[2] - 0.05 * diff(x_range)
  y_pos <- y_range[1] + 0.05 * diff(y_range)
  
  # Format dates for the title
  formatted_start <- format(common_start, "%d. %b")
  formatted_end <- format(common_end, "%d. %b")
  plot_title <- paste("Energy balance between", formatted_start, "and", formatted_end)
  
  # Prepare plot object
  plot <- ggplot2::ggplot(merged_df, ggplot2::aes(x = x_data, y = y_data)) +
    ggplot2::geom_point(alpha = 0.6, color = "grey") +
    ggplot2::geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    ggplot2::annotate("text", x = x_pos, y = y_pos, 
                      label = paste0("Slope (m): ", slope), hjust = 1, size = 4) +
    ggplot2::labs(
      x = expression(R[n] - G ~ "["*W~m^{-2}*"]"),
      y = expression(LE + H ~ "["*W~m^{-2}*"]"),
      title = plot_title
    ) 
  
  # Additional styling
  plot <- plot + ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  
  if (save) {
    # Extract start and end date from merged_df$date_time
    start_date <- format(min(merged_df$date_time, na.rm = TRUE), "%d%m")
    end_date <- format(max(merged_df$date_time, na.rm = TRUE), "%d%m")
    
    file_name <- paste0("Plots/EnergyBalance_", start_date, "_", end_date, "_ggplot.png")
    
    ggplot2::ggsave(
      filename = file_name,
      plot = plot,
      width = 10,
      height = 6
    )
  }
  
  return(plot)
}


# Plot ustar ---------------------------------------------------------------
## diurnal range to compare with gap filling success 
# u_star is u. in the EC data 
# use x = "time" and y = "u."
# Falls jemand eine bessere Idee hat wie man die x-Achse mit nur den 
# vollen Stunden bennenen kann let me know :)

boxplot_diurnal_ustar <- function(ec_data, x, y) {
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(hms)
  
  x_label <- if (x == "time") "Time of day" else NULL 
  y_label <- if (y == "u." | y == "ustar" | y == "uStar" | y == "u_star") 
    expression(italic(u)*"*" ~ "[" * m ~ s^{-1} * "]") else NULL
  title_label <- if (y == "u." | y == "ustar" | y == "uStar" | y == "u_star") 
    expression("Diurnal range of " * italic(u)*"*") else NULL
  
  
  ec_filtered <- ec_data %>% 
    select(x_var = all_of(x), y_var = all_of(y)) %>% 
    filter(y_var >= 0) 
  
  time_labels <- c("00:00", " ", "01:00", " ", "02:00", " ", "03:00"," ", "04:00"," ",
                   "05:00"," ", "06:00"," ", "07:00"," ", "08:00"," ", "09:00"," ",
                   "10:00"," ", "11:00"," ", 
                   "12:00"," ", "13:00"," ", "14:00"," ", "15:00"," ", "16:00", 
                   " ","17:00"," ", "18:00", " ",
                   "19:00"," ", "20:00"," ", "21:00"," ", "22:00"," ", "23:00"," ")
  
  
  p <- ggplot(ec_filtered, aes(x = x_var, y = y_var, group = x_var)) +
    geom_boxplot(fill = "cornflowerblue") +
    scale_x_discrete(breaks = time_labels)+
    theme_minimal() +
    labs(title = title_label,
         x = x_label,
         y = y_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(p)
  
}




# only execute if not loaded in another script
#if (interactive()) {
#  #path='Data/eddypro_advanced_01_full_output_2025-05-26T212140_adv.csv'
#  path='Data/eddypro_advanced_01_full_output_2025-05-28T182940_adv.csv'
#  ec_df <- load_ec_data(path)
#  interactive_dygraph_simple(y_col=ec_df$LE, x_col=ec_df$date_time, show_NA=T, save=F)
#  interactive_dygraph_qc(y_col=ec_df$LE, QC_col=ec_df$qc_LE, x_col=ec_df$date_time, show_NA=T, save=F)
#  noninteractive_ggplot_simple(y_col=ec_df$LE, x_col=ec_df$date_time, show_NA=T, save=T)
#  noninteractive_ggplot_qc(y_col=ec_df$LE, QC_col=ec_df$qc_LE, x_col=ec_df$date_time, show_NA=T, save=T)

#  path_met <- 'Data/MeteoFull_PBioAtmo25_30May.csv'
#  met_df <- load_meteo_data(path_met)
#  plot_energybalance_ggplot(ec_df, met_df)
#}
