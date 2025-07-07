# 3) poster plots --------------------------------------------------------------

# adden: QC flags weg und schwarze Nulllinie. Außerdem keine Legende
# und vertikal stauchen

source("Code/HelperFunctions.R")
data <- read.csv("Data/ec_meteo_combined_gapfilled.csv")
data$date_time <- as.POSIXct(data$date_time)

banderole_plot <- function(y_col, QC_col, x_col, y_min, y_max, y_name, save = FALSE) {
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
    ggplot2::geom_point(ggplot2::aes(color = qc), shape = 16, size = 2, na.rm = TRUE) +  # shape = 4 ist x und = 16 ein Punkt
    ggplot2::scale_color_manual(values = qc_colors, name = "QC Group") +
    ggplot2::theme_minimal() +
    ggplot2::labs(y = y_name, x = "Time") +
    ggplot2::coord_cartesian(ylim = c(y_min, y_max))
  
  # # Additional styling
  # p <- p + ggplot2::theme_minimal() +
  #   ggplot2::theme(
  #     panel.background = ggplot2::element_rect(fill = "white", color = NA),
  #     plot.background = ggplot2::element_rect(fill = "white", color = NA)
  #   )
  
  # Save if requested
  if (save) {
    # Extract start and end date from x_col
    start_date <- format(min(x_col, na.rm = TRUE), "%d%m")
    end_date <- format(max(x_col, na.rm = TRUE), "%d%m")
    
    # Create filename with y_col_name and date range
    file_name <- paste0("Plots/", y_col_name, "_", start_date, "_", end_date, "_banderole_qc.png")
    
    ggplot2::ggsave(
      filename = file_name,
      plot = p,
      width = 42,    # in Zoll ≈ 1067 mm (knapp unter A0-Breite)
      height = 4,    # in Zoll ≈ 100 mm (dünne Banderole)
      units = "in",
      dpi = 500
    )
  }
  
  return(p)
}

banderole_plot(y_col=data$co2_flux_f, QC_col=data$qc_co2_flux, x_col=data$date_time,
               y_name=expression(CO[2]~"flux ["*mu*"mol"~m^-2~s^-1*"]"),
               y_min=-20, y_max=20, save=T)

banderole_plot(y_col=data$h2o_flux_f, QC_col=data$qc_h2o_flux, x_col=data$date_time,
               y_name=expression(H[2]*O~"flux ["*mu*"mol"~m^-2~s^-1*"]"),
               y_min=0, y_max=10, save=T)
