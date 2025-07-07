library(purrr)


sFillInit <- function(
    sDATA,
    ### Initializes data frame sTEMP for newly generated gap filled data and qualifiers.
  ##title<< sEddyProc$sFillInit - Initialize gap filling
  Var.s                   ##<< Variable to be filled
  , QFVar.s = 'none'      ##<< Quality flag of variable to be filled
  , QFValue.n = NA_real_  ##<< Value of quality flag for _good_ (original) data,
  ## other data is set to missing
  , FillAll.b = TRUE      ##<< Fill all values to estimate uncertainties
  #! , QF.V.b = TRUE      ##<< boolean vector of length nRow(sData), to allow
  # # specifying bad data directly (those entries that are set to FALSE)
) {
  ##author<< AMM
  'Initializes data frame sTEMP for newly generated gap filled data and qualifiers.'
  sTEMP <- sDATA
  # Check variable to fill and apply quality flag
  fCheckColNames(cbind(sDATA, sTEMP), c(Var.s, QFVar.s), 'sFillInit')
  Var.V.n <- fSetQF(cbind(sDATA, sTEMP), Var.s, QFVar.s, QFValue.n, 'sFillInit')
  #! Var.V.n[QF.V.b == FALSE]   <-   NA_real_
  
  # Abort if variable to be filled contains no data
  if (sum(!is.na(Var.V.n)) == 0) {
    warning('sFillInit::: Variable to be filled (', Var.s, ') contains no data at all!')
    return(-111)
  }
  ##details<<
  ## Description of newly generated variables with gap filled data and qualifiers: \cr
  ##details<<
  ## VAR\emph{_orig} - Original values used for gap filling \cr
  ## VAR\emph{_f   } - Original values and gaps filled with mean of selected
  ## datapoints (condition depending on gap filling method) \cr
  ## VAR\emph{_fqc} - Quality flag assigned depending on gap filling method and
  ## window length (0 = original data, 1 = most reliable, 2 = medium, 3 = least reliable) \cr
  ## VAR\emph{_fall} - All values considered as gaps (for uncertainty estimates) \cr
  ## VAR\emph{_fall_qc} - Quality flag assigned depending on gap filling method
  ## and window length (1 = most reliable, 2 = medium, 3 = least reliable) \cr
  ## VAR\emph{_fnum} - Number of datapoints used for gap-filling \cr
  ## VAR\emph{_fsd} - Standard deviation of datapoints used for gap
  ## filling (uncertainty) \cr
  ## VAR\emph{_fmeth} - Method used for gap filling (1 = similar meteo
  ## condition (sFillLUT with Rg, VPD, Tair), 2 = similar meteo
  ## (sFillLUT with Rg only), 3 = mean diurnal course (sFillMDC)) \cr
  ## VAR\emph{_fwin} - Full window length used for gap filling \cr
  
  lTEMP <- data.frame(
    VAR_orig = Var.V.n      # Original values of variable VAR used for gap filling
    , VAR_f = NA_real_       # Original values and filled gaps
    , VAR_fqc = NA_real_     # Quality flag assigned depending on gap filling
    ## method and window length
    , VAR_fall = NA_real_    # All values considered as gaps
    ## (for uncertainty estimates)
    , VAR_fall_qc = NA_real_ # Quality flag assigned depending on gap filling
    ## method and window length
    , VAR_fnum = NA_real_    # Number of datapoints used for gap-filling
    , VAR_fsd = NA_real_     # Standard deviation of data points used for filling
    , VAR_fmeth = NA_real_   # Method used for gap filling
    , VAR_fwin = NA_real_    # Full window length used for gap filling
  )
  
  # Set fqc to zero for original values
  lTEMP$VAR_f <- lTEMP$VAR_orig
  lTEMP$VAR_fqc <- ifelse(!is.na(lTEMP$VAR_orig), 0, NA_real_)
  
  # Set filling of only gaps
  if (FillAll.b == FALSE) lTEMP$VAR_fall <- lTEMP$VAR_orig #"prefill" with original data
  
  # Add units
  attr(lTEMP$VAR_f, 'units') <- attr(Var.V.n, 'units')
  attr(lTEMP$VAR_fall, 'units') <- attr(Var.V.n, 'units')
  attr(lTEMP$VAR_fsd, 'units') <- attr(Var.V.n, 'units')
  attr(lTEMP$VAR_fwin, 'units') <- 'days'
  
  ##details<<
  ## Long gaps (larger than 60 days) are not filled.
  #! Not congruent with PV-Wave, there the code is performed on single years
  #only with long gaps of 60 days in the beginning or end skipped.
  GapLength.V.n <- fCalcLengthOfGaps(lTEMP$VAR_orig)
  kMaxGap.n <- 48 * 60 #Halfhours in 60 days
  while (max(GapLength.V.n) > kMaxGap.n) {
    #Flag long gap with -9999.0
    End.i <- which(GapLength.V.n == max(GapLength.V.n))
    Start.i <- End.i - max(GapLength.V.n) + 1
    lTEMP$VAR_fall[Start.i:End.i] <- -9999.0 #Set to -9999.0 as a flag for long gaps
    GapLength.V.n[Start.i:End.i] <- -1 #Set to -1 since accounted for
    warning(
      'sMDSGapFill::: The long gap between position ', Start.i, ' and '
      , End.i, ' will not be filled!')
  }
  
  if (FillAll.b == T) {
    message(
      'Initialized variable \'', Var.s, '\' with ', sum(is.na(lTEMP$VAR_orig))
      , ' real gaps for gap filling of all ', sum(is.na(lTEMP$VAR_fall))
      , ' values (to estimate uncertainties).')
  } else {
    message('Initialized variable \'', Var.s, '\' with ', sum(is.na(lTEMP$VAR_orig)),
            ' real gaps for gap filling.')
  }
  
  # twutz: error prone if sTEMP already contains columns of lTEMP
  sTEMP <<- data.frame(c(sTEMP, lTEMP))
  return(invisible(NULL))
}


MDSGapFill <- function(
    ### MDS gap filling algorithm adapted after the PV-Wave code and paper by Markus Reichstein.
  
  sDATA, #<-- dataframe with data to be filled
  Var = Var.s                 ##<< Variable to be filled
  , QFVar = if (!missing(QFVar.s)) QFVar.s else 'none'       ##<< Quality flag
  ##  of variable to be filled
  , QFValue = if (!missing(QFValue.n)) QFValue.n else NA_real_   ##<< Value of
  ##  quality flag for _good_ (original) data, other data is set to missing
  , V1 = if (!missing(V1.s)) V1.s else 'Rg'            ##<< Condition variable 1
  ## (default: Global radiation 'Rg' in  W m-2)
  , T1 = if (!missing(T1.n)) T1.n else 50              ##<< Tolerance interval 1
  ## (default: 50 W m-2)
  , V2 = if (!missing(V2.s)) V2.s else 'VPD'           ##<< Condition variable 2
  ## (default: Vapour pressure deficit 'VPD' in hPa)
  , T2 = if (!missing(T2.n)) T2.n else 5               ##<< Tolerance interval 2
  ## (default: 5 hPa)
  , V3 = if (!missing(V3.s)) V3.s else 'Tair'          ##<< Condition variable 3
  ## (default: Air temperature 'Tair' in degC)
  , T3 = if (!missing(T3.n)) T3.n else 2.5             ##<< Tolerance interval 3
  ## (default: 2.5 degC)
  , FillAll = if (!missing(FillAll.b)) FillAll.b else TRUE       ##<< Fill
  ##  all values to estimate uncertainties
  , isVerbose = if (!missing(Verbose.b)) Verbose.b else TRUE       ##<< Print
  ##  status information to screen
  , suffix = if (!missing(Suffix.s)) Suffix.s else ''	      ##<< String
  ##  suffix needed for different processing setups on the same dataset
  ## (for explanations see below)
  , minNWarnRunLength = ##<< scalar integer:
    ## warn if number of subsequent
    ## numerically equal values exceeds this number.
    ## Set to Inf or NA for no warnings.
    ## defaults for "NEE" to records across 4 hours and no warning for others.
    if (Var == "NEE") 4 * 48/24 else NA_integer_
  , Var.s      ##<< deprecated
  , QFVar.s    ##<< deprecated
  , QFValue.n  ##<< deprecated
  , V1.s ##<< deprecated
  , T1.n ##<< deprecated
  , V2.s ##<< deprecated
  , T2.n ##<< deprecated
  , V3.s ##<< deprecated
  , T3.n ##<< deprecated
  , FillAll.b ##<< deprecated
  , Verbose.b ##<< deprecated
  , Suffix.s ##<< deprecated
  #! , QF.V.b = TRUE        ##<< boolean vector of length nRow(sData),
  ## to allow specifying bad data directly (those entries that are set to FALSE)
  , method = "Reichstein05" ##<< specify "Vekuri23" to use the skewness-bias
  ## reducing variant
) {
  varNamesDepr <- c(
    "Var.s","QFVar.s","QFValue.n","V1.s","T1.n"
    ,"V2.s","T2.n","V3.s","T3.n","FillAll.b","Verbose.b","Suffix.s")
  varNamesNew <- c(
    "Var","QFVar","QFValue","V1","T1"
    ,"V2","T2","V3","T3","FillAll","isVerbose","suffix")
  iDepr = which(!c(
    missing(Var.s),missing(QFVar.s),missing(QFValue.n),missing(V1.s)
    ,missing(T1.n),missing(V2.s),missing(T2.n),missing(V3.s),missing(T3.n)
    ,missing(FillAll.b),missing(Verbose.b),missing(Suffix.s)))
  if (length(iDepr)) warning(
    "Arguments names ",paste(varNamesDepr[iDepr], collapse = ",")
    ," have been deprecated."
    ," Please, use instead ", paste(varNamesNew[iDepr], collapse = ","))
  ##author<< AMM, TW
  ##references<<
  ## Reichstein, M. et al. (2005) On the separation of net ecosystem exchange
  ## into assimilation and ecosystem respiration: review and improved algorithm.
  ## Global Change Biology, 11, 1424-1439.
  TimeStart.p <- Sys.time()
  ##details<<
  ## Initialize temporal data frame sTEMP for newly generated gap filled data and
  ## qualifiers, see \code{\link{sEddyProc_sFillInit}} for explanations on suffixes.
  # sTEMP <<- sTEMP[, 1L, drop = FALSE]
  if (!is.null(sFillInit(sDATA, Var, QFVar, QFValue, FillAll)) ) #! , QF.V.b = QF.V.b)) )
    return(invisible(-111)) # Abort gap filling if initialization of sTEMP failed
  ##details<<
  ## Runs of numerically equal numbers hint to problems of the data and cause
  ## unreasonable estimates of uncertainty. This routine warns the user.
  if (is.finite(minNWarnRunLength) & (nrow(sTEMP) >= minNWarnRunLength)) {
    rl <- .runLength(as.vector(sTEMP$VAR_orig), minNRunLength = minNWarnRunLength)
    if (length(rl$index)) {
      rlSorted <- rl[rev(order(rl$nRep)),,drop = FALSE]
      warning(
        "Variable ", Var, " contains long runs of numerically equal numbers."
        , " Longest of ", rlSorted$nRep[1], " repeats of value "
        , sTEMP$VAR_orig[ rlSorted$index[1] ]
        , " starts at index ", rlSorted$index[1]
      )
    }
  }
  #+++ Handling of special cases of meteo condition variables V1, V2, V3
  # If variables are at default values but do not exist as columns, set to 'none'
  # (= disabled identifier).
  # This allows running MDS with less variables than prescribed in the default setting.
  # If meteo condition variable are same as variable to fill, also set to 'none'.
  # This prevents filling artificial gaps (for uncertainty estimates) with itself
  # as meteo condition variable.
  #! Attention: Non-congruent with MR PV-Wave. There artificial gaps in
  #Rg, VPD, Tair are filled with itself.
  if ( (V1 ==   'Rg' && !(V1 %in% c(colnames(sDATA)))) || (V1 == Var) )   V1 <- 'none'
  if ( (V2 ==  'VPD' && !(V2 %in% c(colnames(sDATA)))) || (V2 == Var) )   V2 <- 'none'
  if ( (V3 == 'Tair' && !(V3 %in% c(colnames(sDATA)))) || (V3 == Var) )   V3 <- 'none'
  # Check column names (with 'none' as dummy)
  # (Numeric type and plausibility have been checked on initialization of sEddyProc)
  fCheckColNames(cbind(sDATA, sTEMP), c(V1, V2, V3), 'sMDSGapFill')
  # Check tolerance entries (if condition variable is not 'none')
  NoneCols.b <- c(V1, V2, V3) %in% 'none'
  if (!fCheckValNum(T1) || !fCheckValNum(T2) || !fCheckValNum(T3) ) stop(
    'sMDSGapFill::: T1, T2, T3, T4.n, T5.n must be numeric '
    , '(if not specified, set to NA_real_)!')
  if (sum(is.na(c(T1, T2, T3)[!NoneCols.b])) ) stop(
    'sMDSGapFill::: If condition variable is specified (dummy name is \'none\'), '
    , 'the tolerance interval must be specified!')
  # Run gap filling scheme depending on auxiliary meteo data availability
  ##details<<
  ## MDS gap filling algorithm calls the subroutines Look Up Table
  ## \code{\link{sEddyProc_sFillLUT}}
  ## and Mean Diurnal Course \code{\link{sEddyProc_sFillMDC}} with different
  ## window sizes as described in the reference.
  ##details<<
  ## Vekari et al. 2023 proposed a modification that reduces bias
  ## that results from skewed distribution of radiation Rg.
  ## In order to use this modification, specify \code{method="Vekuri23"}.
  ##details<<
  calculate_gapstats <- if (method == "Reichstein05") {
    calculate_gapstats_Reichstein05
  } else if (method == "Vekuri23") {
    if (V1 != "Rg") warning(paste0(
      "Method Vekuri23 requires incoming shortware radiation (Rg), ",
      "but first covariate is called ", V1))
    calculate_gapstats_Vekuri23
  } else
    stop(paste0(
      "unknown MDS method '",method,"'. Specify 'Reichstein05' or 'Vekuri23'"))
  ## To run dataset only with MDC algorithm \code{\link{sEddyProc_sFillMDC}},
  ## set condition variable V1 to 'none'.
  # Check availablility of meteorological data for LUT
  Met.n <-
    if (
      V1 != 'none' && V2 != 'none' && V3 != 'none'
      && sum(!is.na(sDATA[, V1])) != 0 && sum(!is.na(sDATA[, V2])) != 0 &&
      sum(!is.na(sDATA[, V3])) != 0
    ) {
      #All three meteo conditions are available and valid to use:
      message(
        'Full MDS algorithm for gap filling of \''
        , attr(sTEMP$VAR_f, 'varnames'), '\' with LUT('
        , V1, ', ', V2, ', ', V3, ') and MDC.')
      3
    } else if (V1 != 'none' && sum(!is.na(sDATA[, V1])) != 0) {
      #Only one meteo condition available for LUT
      message(
        'Limited MDS algorithm for gap filling of \''
        , attr(sTEMP$VAR_f, 'varnames'), '\' with LUT(', V1, ' only) and MDC.')
      1
    } else {
      #No meteo condition available (use MDC only)
      message(
        'Restriced MDS algorithm for gap filling of \''
        , attr(sTEMP$VAR_f, 'varnames')
        , '\' with no meteo conditions and hence only MDC.')
      if (Var != 'Rg') warning(
        'sMDSGapFill::: No meteo available for MDS gap filling!')
      0
    }
  #+++ Full MDS algorithm
  # Step 1: Look-up table (method 1) with window size +-7 days
  if (Met.n == 3) sTEMP <- sFillLUT(sDATA,
    7, V1, T1, V2, T2, V3, T3, Verbose.b = isVerbose,
    calculate_gapstats = calculate_gapstats)
  # Step 2: Look-up table (method 1) with window size +-14 days
  if (Met.n == 3) sTEMP <- sFillLUT(sDATA,
    14, V1, T1, V2, T2, V3, T3, Verbose.b = isVerbose,
    calculate_gapstats = calculate_gapstats)
  # Step 3: Look-up table, Rg only (method 2) with window size +-7 days,
  if (Met.n == 3 || Met.n == 1) sTEMP <- sFillLUT(sDATA,
    7, V1, T1, Verbose.b = isVerbose,
    calculate_gapstats = calculate_gapstats)
  # Step 4: Mean diurnal course (method 3) with window size 0 (same day)
  sTEMP <- sFillMDC(sDATA,0, Verbose.b = isVerbose)
  # Step 5: Mean diurnal course (method 3) with window size +-1, +-2 days
  sTEMP <- sFillMDC(sDATA,1, Verbose.b = isVerbose)
  sTEMP <- sFillMDC(sDATA,2, Verbose.b = isVerbose)
  # Step 6: Look-up table (method 1) with window size +-21, +-28, ..., +-70
  if (Met.n == 3) for (WinDays.i in seq(21, 70, 7) ) sTEMP <- sFillLUT(sDATA,
    WinDays.i, V1, T1, V2, T2, V3, T3, Verbose.b = isVerbose,
    calculate_gapstats = calculate_gapstats)
  # Step 7: Look-up table (method 2) with window size +-14, +-21, ..., +-70
  if (Met.n == 3 || Met.n == 1) for (WinDays.i in seq(14, 70, 7) ) sTEMP <- sFillLUT(sDATA,
    WinDays.i, V1, T1, Verbose.b = isVerbose,
    calculate_gapstats = calculate_gapstats)
  # Step 8: Mean diurnal course (method 3) with window size +-7, +-14, ..., +-210 days
  for (WinDays.i in seq(7, 210, 7) ) sTEMP <- sFillMDC(sDATA,WinDays.i, Verbose.b = isVerbose)
  
  # Set long gaps again to NA
  sTEMP$VAR_fall <<- suppressMessages(fConvertGapsToNA(sTEMP$VAR_fall))
  
  # Message on gap filling
  TimeDiff.p <- as.numeric(Sys.time()) - as.numeric(TimeStart.p)
  message(
    'Finished gap filling of \'', Var, '\' in ', floor(TimeDiff.p)
    , ' seconds. Artificial gaps filled: '
    , length(sTEMP$VAR_fall) - sum(is.na(sTEMP$VAR_fall))
    , ', real gaps filled: ', sum(is.na(sTEMP$VAR_orig))
    , ', unfilled (long) gaps: ', sum(is.na(sTEMP$VAR_fall)), '.')
  
  ##details<< \describe{\item{Different processing setups on the same dataset}{
  ## Attention: When processing the same site data set with different setups for
  ## the gap filling or flux partitioning
  ## (e.g. due to different ustar filters),
  ## a string suffix is needed! This suffix is added to the result column names
  ## to distinguish the results of the different setups.
  ## }}
  # Rename new columns generated during gap filling
  colnames(sTEMP) <- gsub("VAR", Var, colnames(sTEMP))
  
  # Check for duplicate columns (to detect if different processing
  # setups were executed without different suffixes)
  if (length(names(which(table(colnames(sTEMP)) > 1))) )  warning(
    'sMDSGapFill::: Duplicated columns found! Please specify different suffix '
    , 'when processing different setups on the same dataset!')
  return(sTEMP)

}



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @export
sFillLUT <- function(
    ### Look-Up Table (LUT) algorithm of up to five conditions within prescribed window size
  sDATA,
  WinDays.i             ##<< Window size for filling in days
  , V1.s = 'none'          ##<< Condition variable 1
  , T1.n = NA_real_        ##<< Tolerance interval 1
  , V2.s = 'none'          ##<< Condition variable 2
  , T2.n = NA_real_        ##<< Tolerance interval 2
  , V3.s = 'none'          ##<< Condition variable 3
  , T3.n = NA_real_        ##<< Tolerance interval 3
  , V4.s = 'none'          ##<< Condition variable 4
  , T4.n = NA_real_        ##<< Tolerance interval 4
  , V5.s = 'none'          ##<< Condition variable 5
  , T5.n = NA_real_        ##<< Tolerance interval 5
  , Verbose.b = TRUE       ##<< Print status information to screen
  , calculate_gapstats = calculate_gapstats_Reichstein05 ##<< function computing gap-statistics
) {
  ##author<< AMM, TW
  #! Attention: For performance reasons, gap filled values and properties are
  #first written to single variables and local matrix lGF.M
  #! (rather than changing single values in sTEMP which copies the data frame each time!)
  #! Improved algorithm speed by more than a factor of 10 (maybe even 100...)
  lGF.M <- matrix(NA_real_, nrow = 0, ncol = 7, dimnames = list(
    NULL, c('index', 'mean', 'fnum', 'fsd', 'fmeth', 'fwin', 'fqc')))
  
  # Check if sTEMP has been initialized with new VAR_ columns
  if (!exists('VAR_f', sTEMP) ) stop(
    'sFillLUT::: Temporal data frame sTEMP for processing results has not '
    , 'been initalized with sFillInit!')
  
  # Determine gap positions
  ToBeFilled.V.i <- which(is.na(sTEMP$VAR_fall))
  if (length(ToBeFilled.V.i) > 0) {
    for (Pos.i in 1:length(ToBeFilled.V.i) ) {
      # Message on progress if wanted
      NoneCols.b <- c(V1.s, V2.s, V3.s, V4.s, V5.s) %in% 'none'
      if (Verbose.b && Pos.i == 1)  message(
        'Look up table with window size of ', WinDays.i, ' days with '
        , paste(c(V1.s, V2.s, V3.s, V4.s, V5.s)[!NoneCols.b], collapse = ' '))
      # Set window size
      Gap.i   <- ToBeFilled.V.i[Pos.i]
      if (T == T) {
        #! Non-congruent with MR PV-Wave
        Start.i <- Gap.i - (WinDays.i * 48)
        End.i   <- Gap.i + (WinDays.i * 48)
      } else {
        #! Code congruent with MR PV-Wave --> window size minus 1
        Start.i <- Gap.i - (WinDays.i * 48 - 1)
        End.i   <- Gap.i + (WinDays.i * 48 - 1)
      }
      
      if (Start.i <= 0) Start.i <- 1
      if (End.i > nrow(sTEMP) ) End.i <- nrow(sTEMP)
      
      # Special treatment of Rg to be congruent with MR PV-Wave,
      # not mentioned in paper
      isV1Rg <- grepl('Rg', V1.s)
      T1red.n <- if (isV1Rg) {
        # Reduce tolerance of radiation if variable name contains 'Rg' to
        # [20, 50] depending on measurement
        # at least 20, if Rg>20 then min(Rg, tol=50), larger than 20 for 20<Rg<50
        max(min(T1.n, sDATA[Gap.i, V1.s], na.rm = T), 20, na.rm = T)
      } else {
        T1.n
      }
      
      # For performance reasons, write sDATA subrange into vectors
      # (speed up about factor of 1.5)
      V1.V.n <- sDATA[Start.i:End.i, V1.s]
      V2.V.n <- sDATA[Start.i:End.i, V2.s]
      V3.V.n <- sDATA[Start.i:End.i, V3.s]
      V4.V.n <- sDATA[Start.i:End.i, V4.s]
      V5.V.n <- sDATA[Start.i:End.i, V5.s]
      SubGap.i <- Gap.i - (Start.i - 1)
      
      # Set LUT fill condition
      Rows.V.b <- !is.na(sTEMP$VAR_orig[Start.i:End.i])
      if (V1.s != 'none')
        Rows.V.b <- Rows.V.b & abs(V1.V.n - V1.V.n[SubGap.i]) < T1red.n  & !is.na(V1.V.n)
      if (V2.s != 'none')
        Rows.V.b <- Rows.V.b & abs(V2.V.n - V2.V.n[SubGap.i]) < T2.n  & !is.na(V2.V.n)
      if (V3.s != 'none')
        Rows.V.b <- Rows.V.b & abs(V3.V.n - V3.V.n[SubGap.i]) < T3.n  & !is.na(V3.V.n)
      if (V4.s != 'none')
        Rows.V.b <- Rows.V.b & abs(V4.V.n - V4.V.n[SubGap.i]) < T4.n  & !is.na(V4.V.n)
      if (V5.s != 'none')
        Rows.V.b <- Rows.V.b & abs(V5.V.n - V5.V.n[SubGap.i]) < T5.n  & !is.na(V5.V.n)
      lLUT.V.n <- subset(sTEMP$VAR_orig[Start.i:End.i], Rows.V.b)
      
      # Here using <=, alternative could use strict <
      #isV1belowT1 <- V1.V.n[Rows.V.b] < V1.V.n[SubGap.i]
      # usually fewer values with lower Rg -> use <=
      isV1belowT1 <- V1.V.n[Rows.V.b] <= V1.V.n[SubGap.i]
      
      is_nighttime <- if (isV1Rg) {
        # determine nighttime by radiation below 10W/m2 corresponding to uStar filter
        sDATA[Gap.i, V1.s] < 10
      } else {
        # without radiation: assume daytime-record
        FALSE
      }
      
      # If enough available data, fill gap
      if (length(lLUT.V.n) > 1) {
        lVAR_index.i <- Gap.i
        gapstats <- calculate_gapstats(lLUT.V.n, isV1belowT1, is_nighttime)
        lVAR_mean.n <- gapstats[1]
        lVAR_fnum.n <- gapstats[2]
        lVAR_fsd.n <- gapstats[3]
        
        #Set window size and quality flag
        ##details<< \describe{\item{Quality flags}{
        ## \itemize{
        ## \item 1: at least one variable and nDay <= 14
        ## \item 2: three variables and nDay in [14,56)
        ## or one variable and nDay in  [14,28)
        ## \item 3: three variables and nDay > 56
        ## or one variable and nDay > 28
        ## }
        ## }}
        #! Full window length, congruent with MR PV-Wave, in paper
        #single window sizes stated
        lVAR_fwin.n  <- 2 * WinDays.i
        lVAR_fmeth.n <- NA_real_; lVAR_fqc.n <- NA_real_;
        if (V1.s != 'none' && V2.s != 'none' && V3.s != 'none') {
          #Three conditions
          lVAR_fmeth.n <- 1
          #! Limit '14' congruent with MR PV-Wave, in paper different limit
          #of '28' (stated as single window size of 14 days)
          if (lVAR_fwin.n <= 14) lVAR_fqc.n <- 1
          if (lVAR_fwin.n >  14 & lVAR_fwin.n <= 56) lVAR_fqc.n <- 2
          if (lVAR_fwin.n >  56) lVAR_fqc.n <- 3
        }
        if (V1.s != 'none' && V2.s == 'none' && V3.s == 'none') {
          #One condition only
          lVAR_fmeth.n <- 2
          if (lVAR_fwin.n <= 14) lVAR_fqc.n <- 1
          if (lVAR_fwin.n >  14 & lVAR_fwin.n <= 28) lVAR_fqc.n <- 2
          if (lVAR_fwin.n >  28) lVAR_fqc.n <- 3
        }
        lGF.M <- rbind(lGF.M, c(
          lVAR_index.i, lVAR_mean.n, lVAR_fnum.n, lVAR_fsd.n
          , lVAR_fmeth.n, lVAR_fwin.n, lVAR_fqc.n))
      }
      if (Verbose.b && Pos.i %% 100 == 0)  message('.', appendLF = FALSE)
      if (Verbose.b && Pos.i %% 6000 == 0) message('\n.', appendLF = FALSE)
    }
    if (Verbose.b) message('', nrow(lGF.M))
  }
  # Copy gap filled values and properties to sTEMP
  if (nrow(lGF.M) > 0) {
    # Fill all rows in VAR_fall and co
    sTEMP[lGF.M[, 'index'], c(
      'VAR_fall', 'VAR_fnum', 'VAR_fsd', 'VAR_fmeth', 'VAR_fwin', 'VAR_fall_qc')] <<-
      lGF.M[, c('mean', 'fnum', 'fsd', 'fmeth', 'fwin', 'fqc')]
    # Only fill gaps in VAR_f and VAR_fqc
    Gaps.b <- is.na(sTEMP[lGF.M[, 'index'], 'VAR_f'])
    sTEMP[lGF.M[, 'index'], c('VAR_f', 'VAR_fqc')][Gaps.b, ] <<-
      as.data.frame(lGF.M[, c('mean', 'fqc') , drop = FALSE])[Gaps.b, ]
  }
  
  #Other columns are specific for full MR MDS algorithm
  return(sTEMP[, c(
    'VAR_orig', 'VAR_f', 'VAR_fall', 'VAR_fnum', 'VAR_fsd', 'VAR_fwin')])
  ##value<<
  ## LUT filling results in sTEMP data frame.
}


calculate_gapstats_Reichstein05 <- function(lLUT.V.n, ...) {
  lVAR_mean.n <- mean(lLUT.V.n)
  lVAR_fnum.n <- length(lLUT.V.n)
  lVAR_fsd.n <- sd(lLUT.V.n)
  c(lVAR_mean.n, lVAR_fnum.n, lVAR_fsd.n)
}

calculate_gapstats_Vekuri23 <- function(
    ### Compute the gap statistic by separately aggregating reocrds lower or higher than the gap-covariate
  lLUT.V.n, ##<< numeric vector of observations
  islower,  ##<< logical vector indicating if the observation was lower than gaps covariate
  is_nighttime, ##<< logical scalar indicating whether record to fill is nighttime
  n_min=1   ##<< minimum number of observations in each group.
  ## If there are fewer observations in on the groups
  ## return the (skewness-biased) Reichstein05 estimate.
) {
  ##details<<
  ## In order to reduce bias due to skewed distributions of covariate.
  if (isTRUE(is_nighttime))
    return(calculate_gapstats_Reichstein05(lLUT.V.n))
  vals_lower <- lLUT.V.n[islower]
  vals_higher <- lLUT.V.n[!islower]
  if ((length(vals_lower) < n_min) || (length(vals_higher) < n_min))
    return(calculate_gapstats_Reichstein05(lLUT.V.n))
  mean_lower <- mean(vals_lower)
  mean_higher <- mean(vals_higher)
  lVAR_mean.n <- (mean_lower  + mean_higher)/2
  lVAR_fnum.n <- length(lLUT.V.n)
  lVAR_fsd.n <- sd(lLUT.V.n) # keep uncertainty of population instead of estimate
  # ##details<<
  # ## Error propagation assumes upper and lower aggregate to be uncorrelated.
  # ## If uncertainty of one group could not be estimated, e.g. if it has only
  # ## one observation, then the uncertainty of
  # ## the mean across the two values is returned.
  # var_lower <- var(vals_lower)
  # var_higher <- var(vals_higher)
  # lVAR_fsd.n <- if (!is.finite(var_lower) || !is.finite(var_higher)) {
  #   sd(c(mean_lower, mean_higher))
  # } else {
  #   sqrt(var_lower + var_higher) / 2
  # }
  c(lVAR_mean.n, lVAR_fnum.n, lVAR_fsd.n)
}


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#' @export
sFillMDC <- function(
    ### Mean Diurnal Course (MDC) algorithm based on average values within +/- one hour of adjacent days
  sDATA,
  WinDays.i            ##<< Window size for filling in days
  , Verbose.b = TRUE   ##<< Print status information to screen
) {
  ##author<< AMM
  #! Attention: For performance reasons, gap filled values and properties are first written to single
  #! variables and local matrix lGF.M
  #! (rather than changing single values in sTEMP which copies the data frame each time!)
  #! Improved algorithm speed by more than a factor of 10 (maybe even 100...)
  lGF.M <- matrix(NA_real_, nrow = 0, ncol = 7, dimnames = list(NULL, c(
    'index', 'mean', 'fnum', 'fsd', 'fmeth', 'fwin', 'fqc')))
  
  # Determine gap positions
  ToBeFilled.V.i <- which(is.na(sTEMP$VAR_fall))
  if (length(ToBeFilled.V.i) > 0) {
    for (Pos.i in 1:length(ToBeFilled.V.i)) {
      # Message on progress if wanted
      if (Verbose.b && Pos.i == 1) message(
        'Mean diurnal course with window size of ', WinDays.i, ' days: .', sep = '')
      
      # Set index within window size
      Gap.i   <- ToBeFilled.V.i[Pos.i]
      Index.V.i <- numeric(0)
      for (Day.i in (0:WinDays.i))
        if (Day.i == 0) {
          Index.V.i <- c(Index.V.i, Gap.i + (-2:2))
        } else {
          Index.V.i <- c(
            Index.V.i, Gap.i + c(-Day.i * 48 + (-2:2))
            , Gap.i + c(Day.i * 48 + (-2:2)))
        }
      Index.V.i <- Index.V.i[Index.V.i > 0 & Index.V.i <= nrow(sTEMP)]
      
      # If enough available data, fill gap
      #iRecsSimilar <- Index.V.i[!is.na(sTEMP$VAR_orig[Index.V.i]) ]
      lMDC.V.n <- subset(sTEMP$VAR_orig[Index.V.i], !is.na(sTEMP$VAR_orig[Index.V.i]))
      
      if (length(lMDC.V.n) > 1) {
        #if (Gap.i == 15167  ) recover()
        
        lVAR_index.i <- Gap.i
        lVAR_mean.n <- mean(lMDC.V.n)
        lVAR_fnum.n  <- length(lMDC.V.n)
        lVAR_fsd.n  <- sd(lMDC.V.n)
        lVAR_fmeth.n  <- 3
        
        #Set window size and quality flag
        if (T == T) {
          #! Non-congruent with MR PV-Wave
          ##! Full window length, not congruent with MR PV-Wave (see below),
          ##in paper single window sizes stated
          lVAR_fwin.n <- 2 * WinDays.i + 1
        } else {
          #! Code if required to be congruent with MR PV-Wave -->
          #window calculation changes depending on size
          lVAR_fwin.n <- if (WinDays.i < 7) {
            2 * WinDays.i + 1
          } else {
            WinDays.i + 1
          }
        }
        
        ##details<< \describe{\item{Quality flag}{
        ## \itemize{
        ## \item 1: nDay <= 1
        ## \item 2: nDay [2,5)
        ## \item 3: nDay > 5
        ## }
        ## }}
        if (lVAR_fwin.n <= 1) lVAR_fqc.n <- 1
        if (lVAR_fwin.n >  1 & lVAR_fwin.n <= 5) lVAR_fqc.n <- 2
        if (lVAR_fwin.n >  5) lVAR_fqc.n <- 3
        
        lGF.M <- rbind(lGF.M, c(
          lVAR_index.i, lVAR_mean.n, lVAR_fnum.n, lVAR_fsd.n
          , lVAR_fmeth.n, lVAR_fwin.n, lVAR_fqc.n))
      }
      if (Verbose.b && Pos.i %% 100 == 0)  message('.', appendLF = FALSE)
      if (Verbose.b && Pos.i %% 6000 == 0) message('\n.', appendLF = FALSE)
    }
    if (Verbose.b) message('', nrow(lGF.M))
  }
  # Copy gap filled values and properties to sTEMP
  if (nrow(lGF.M) > 0) {
    # Fill all rows in VAR_fall and co
    sTEMP[lGF.M[, 'index'], c(
      'VAR_fall', 'VAR_fnum', 'VAR_fsd', 'VAR_fmeth', 'VAR_fwin', 'VAR_fall_qc')] <<-
      lGF.M[, c('mean', 'fnum', 'fsd', 'fmeth', 'fwin', 'fqc')]
    # Only fill gaps in VAR_f and VAR_fqc
    Gaps.b <- is.na(sTEMP[lGF.M[, 'index'], 'VAR_f'])
    # twutz: inserted drop = FALSE, otherwise one-row matrix was not
    # converted to data.frame correctly
    sTEMP[lGF.M[, 'index'], c('VAR_f', 'VAR_fqc')][Gaps.b, ] <<-
      as.data.frame(lGF.M[, c('mean', 'fqc') , drop = FALSE])[Gaps.b, ]
  }
  
  #Other columns are specific for full MR MDS algorithm
  return(sTEMP[, c(
    'VAR_orig', 'VAR_f', 'VAR_fall', 'VAR_fnum', 'VAR_fsd', 'VAR_fwin')])
  ##value<<
  ## MDC filling results in sTEMP data frame.
}


fCheckColNames <- function(
    ##description<<
  ## Check if specified columns exists in data frame
  Data.F                ##<< Data frame
  , ColNames.V.s         ##<< Vector of variables to be checked
  , CallFunction.s = ''    ##<< Name of function called from
)
  ##author<<
  ## AMM
  # TEST: Data.F <- Date.F.x; ColNames.V.s <- c('Year.n', 'none', 'Month.n', 'test'); CallFunction.s <- 'Dummy'
{
  #Exclude dummy 'none'
  NoneCols.b <- ColNames.V.s %in% 'none'
  #Check if specified columns exist in data frame
  NameCols.b <- ColNames.V.s[!NoneCols.b] %in% colnames(Data.F)
  if (!all(NameCols.b) ) {
    ColNames.s <- paste(ColNames.V.s[!NoneCols.b][!NameCols.b], collapse = ', ', sep = '')
    stop(CallFunction.s, ':::fCheckColNames::: Missing specified columns in dataset: ', ColNames.s, '!')
  }
  
  ##value<<
  ## Function stops on errors.
}

fSetQF <- function(
    ##title<<
  ## Check and set quality flag
  ##description<<
  ## Generate new vector from data and quality flag column.
  Data.F                ##<< Data frame
  , Var.s                ##<< Variable to be filtered
  , QFVar.s              ##<< Quality flag of variable to be filtered
  , QFValue.n            ##<< Numeric value of quality flag for _good_ data, other data is set to missing
  , CallFunction.s = ''    ##<< Name of function called from
)
  ##author<<
  ## AMM
  # TEST: Data.F <- EddyData.F; Var.s <- 'NEE'; QFVar.s <- 'QF'; QFValue.n <- 0; CallFunction.s = 'dummy'
{
  # Check if specified columns exist and are numeric
  SubCallFunc.s <- paste(CallFunction.s, 'fSetQF', sep = ':::')
  fCheckColNames(Data.F, c(Var.s, QFVar.s), SubCallFunc.s)
  fCheckColNum(Data.F, c(Var.s, QFVar.s), SubCallFunc.s)
  
  ##details<<
  ## Quality flag will be applied to variable - unless quality flag variables is called 'none' (dummy).
  if (QFVar.s != 'none') {
    # Check quality flag value
    if (fCheckValNum(QFValue.n) != T)
      stop(CallFunction.s, ':::fSetQF::: Quality flag \'', QFVar.s, '\' has a non-numeric value: ', QFValue.n, '!')
    # Only use data values when good data (quality flag is equal to flag value)
    Var.V.n <- ifelse(Data.F[, QFVar.s] == QFValue.n, Data.F[, Var.s], NA_real_)
    if (sum(!is.na(Var.V.n)) == 0)
      warning(CallFunction.s, ':::fSetQF::: Variable \'', Var.s, '\' contains no data after applying quality flag \'', QFVar.s, '\' with value ', QFValue.n, '!')
  } else {
    # Use all data
    Var.V.n <- Data.F[, Var.s]
    if (sum(!is.na(Var.V.n)) == 0)
      warning(CallFunction.s, ':::fSetQF::: Variable \'', Var.s, '\' contains no data!')
  }
  # Add units
  attr(Var.V.n, 'units') <- attr(Data.F[[Var.s]], 'units')
  attr(Var.V.n, 'varnames') <- if (QFVar.s == 'none') { paste(Var.s, sep = '')
  } else { paste(Var.s, '.', QFVar.s, '_', round(QFValue.n, digits = 3), sep = '') }
  
  Var.V.n
  ##value<<
  ## Numeric vector with _good_ data.
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++ Variable check functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fCheckValString <- function(
    ##description<<
  ## Check if variable is a non-empty character string
  Value.s               ##<< Value to be checked if string
)
  ##author<<
  ## AMM
  ##details<<
  ## See test_CheckValue.R for more details.
{
  if ( length(Value.s) == 0) {
    FALSE
  } else if (!is.na(Value.s) && (!is.character(Value.s) || !nzchar(Value.s)) ) {
    FALSE
  } else {
    TRUE
  }
  ##value<<
  ## Boolean value if true of false.
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fCheckValNum <- function(
    ##description<<
  ## Check if variable is a numeric
  Value.n               ##<< Value to be checked if numeric (but can also be NA of any type)
)
  ##author<<
  ## AMM
  ##details<<
  ## See test_CheckValue.R for more details.
{
  if ( length(Value.n) == 0) {
    FALSE
  } else if (!is.na(Value.n) && !is.numeric(Value.n) ) {
    FALSE
  } else {
    TRUE
  }
  ##value<<
  ## Boolean value if true of false.
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fCheckColNames <- function(
    ##description<<
  ## Check if specified columns exists in data frame
  Data.F                ##<< Data frame
  , ColNames.V.s         ##<< Vector of variables to be checked
  , CallFunction.s = ''    ##<< Name of function called from
)
  ##author<<
  ## AMM
  # TEST: Data.F <- Date.F.x; ColNames.V.s <- c('Year.n', 'none', 'Month.n', 'test'); CallFunction.s <- 'Dummy'
{
  #Exclude dummy 'none'
  NoneCols.b <- ColNames.V.s %in% 'none'
  #Check if specified columns exist in data frame
  NameCols.b <- ColNames.V.s[!NoneCols.b] %in% colnames(Data.F)
  if (!all(NameCols.b) ) {
    ColNames.s <- paste(ColNames.V.s[!NoneCols.b][!NameCols.b], collapse = ', ', sep = '')
    stop(CallFunction.s, ':::fCheckColNames::: Missing specified columns in dataset: ', ColNames.s, '!')
  }
  
  ##value<<
  ## Function stops on errors.
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fCheckColNum <- function(
    ##description<<
  ## Check if specified columns are numeric
  Data.F                ##<< Data frame
  , ColNames.V.s         ##<< Vector of variables to be checked, with 'none' as dummy
  , CallFunction.s = ''    ##<< Name of function called from
  , isWarnMissing = TRUE  ##<< set to FALSE to avoid warning when columns to
  ##  to check are missing
) {
  ##author<< TW
  #Exclude dummy 'none' from checking
  ColNames.V.s <- ColNames.V.s[ColNames.V.s != "none"]
  #Exclude columns not present in dataset from checking
  iMissing <- which(!(ColNames.V.s %in% names(Data.F)))
  colNamesCheck <- if (length(iMissing)) {
    if (isTRUE(isWarnMissing)) warning(
      "missing columns ", paste(ColNames.V.s[iMissing], collapse = ","))
    ColNames.V.s[-iMissing]
  } else ColNames.V.s
  isNotNumeric <- map_lgl(colNamesCheck, ~!is.numeric(Data.F[[.]]))
  if (sum(isNotNumeric)) {
    # report the first occurence of a nonnumeric
    colNameFirst <- colNamesCheck[isNotNumeric][1]
    x <- Data.F[[colNameFirst]]
    xn <- suppressWarnings(as.numeric(x))
    indexFirst <- which(!is.na(x) & is.na(xn))[1]
    stop(CallFunction.s, ":::fCheckColNum::: Detected following columns in ",
         "dataset to be non numeric: ",
         paste(colNamesCheck[isNotNumeric], collapse = ","), '! ',
         "First occurence of non-numeric value at column '",colNameFirst,
         "' at row ", indexFirst, " is '",x[indexFirst],"'.")
  }
  ##value<<
  ## Function stops on errors.
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fCheckOutsideRange <- function(
    ##description<<
  ## Check if specified variable is outside of provided boundaries
  Data.F                ##<< Data frame
  , VarName.s            ##<< Variable (column) name
  , Condition.V.s        ##<< Logical condition for outside values
  , CallFunction.s = ''   ##<< Name of function called from
)
  ##author<<
  ## AMM
  ##details<<
  ## Example of condition structure: c(' <= ', 0) or c(' <= ', 0, '|', '>', 20)
  ## Allowed relational operators: < <= == >= > !=
  ## Allowed logical operators: & |
  # TEST: Data.F <- Date.F.x; VarName.s <- 'Rg';  CallFunction.s <- 'test'; Condition.V.s <- c(' <= ', 0, '|', '>', 20); Condition.V.s <- c(' <= ', 0)
{
  fCheckColNames(Data.F, VarName.s, paste(CallFunction.s, 'fCheckOutsideRange', sep = ':::'))
  fCheckColNum(Data.F, VarName.s, paste(CallFunction.s, 'fCheckOutsideRange',  sep = ':::'))
  Var.V.n <- Data.F[, VarName.s]
  
  # Check condition
  CondText.s <- if (length(Condition.V.s) == 2 && Condition.V.s[1]  %in% c('<', ' <= ', ' == ', ' >= ', '>', ' != ') && nzchar(Condition.V.s[2]) ) {
    # One condition
    paste('Var.V.n ', Condition.V.s[1], ' ', Condition.V.s[2], ' & !is.na(Var.V.n)', sep = '')
  } else if (length(Condition.V.s) == 5 && all(Condition.V.s[c(1, 4)]  %in% c('<', ' <= ', ' == ', ' >= ', '>', ' != '))
             && all(nzchar(Condition.V.s[2]), nzchar(Condition.V.s[5])) && (Condition.V.s[3] %in% c('|', '&')) ) {
    # Two conditions
    paste('(Var.V.n ', Condition.V.s[1], ' ', Condition.V.s[2], '  ', Condition.V.s[3], ' Var.V.n ', Condition.V.s[4], ' ', Condition.V.s[5], ') & !is.na(Var.V.n)', sep = '')
  } else {
    stop(CallFunction.s, ':::fCheckOutsideRange::: Incorrect condition definition: ', paste(Condition.V.s, collapse = ' '), '!')
  }
  
  # Warning message
  Outside.b <- eval(parse(text = CondText.s))
  Outside.n <- sum(Outside.b)
  if (Outside.n > 0)
    warning(CallFunction.s, ':::fCheckOutsideRange::: Variable outside (plausible) range in ', Outside.n,
            ' cases! Invalid values with \'', VarName.s, ' ',
            paste(Condition.V.s, collapse = ' '), '\': ', paste(format(Var.V.n[Outside.b][1:(min(Outside.n, 50))], digits = 2), collapse = ', '), ' ...')
  
  return(invisible(NULL))
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fCheckColPlausibility <- function(
    ##description<<
  ## Check plausibility of common (eddy) variables
  Data.F                ##<< Data frame
  , VarName.V.s          ##<< Variable (column) names
  , CallFunction.s = ''   ##<< Name of function called from
)
  ##author<<
  ## AMM
  # TEST: VarName.V.s <- c('Rg_s'); v.i <- 1
{
  # Check column names
  SubCallFunc.s <- paste(CallFunction.s, 'fCheckColPlausibility', sep = ':::')
  fCheckColNames(Data.F, VarName.V.s, SubCallFunc.s)
  # Strip variable name to before dot '.' (because quality flag setting after dot)
  VarName.V.s <- sub('[.]. * ', '', VarName.V.s)
  
  ##details<<
  ## Variables CONTAINing the following abbreviations are checked for plausibility
  # Separated checks for upper and lower limit to have separate warnings
  for (v.i in 1:length(VarName.V.s)) {
    ## 'Rg' - global radiation, W m-2
    if (grepl('Rg', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 1200), SubCallFunc.s)
    }
    ## 'PotRad' - potential global radiation, W m-2
    if (grepl('PotRad', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 3000), SubCallFunc.s)	#TODO plausible upper bound
    }
    ## 'PPFD' or 'ppfd' - photosynthetic active radiation, umol m-2 s-1
    if (grepl('PPFD', VarName.V.s[v.i], ignore.case = TRUE) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 2500), SubCallFunc.s)
    }
    ## 'PAR' or 'par' - photosynthetic active radiation, umol m-2 s-1
    if (grepl('PAR', VarName.V.s[v.i], ignore.case = TRUE) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 2500), SubCallFunc.s)
    }
    ## 'Ta' - air temperature in degree Celsius
    if (grepl('Ta', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', -70), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 60), SubCallFunc.s)
    }
    ## 'Ts' - soil temperature in degree Celsius
    if (grepl('Ts', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', -20), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 80), SubCallFunc.s)
    }
    ## 'VPD' - vapour pressure deficit in hPa
    if (grepl('VPD', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0), SubCallFunc.s)
      #twutz 2301: at dry sites (ES-Abr) regularly high VPD
      #fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 50), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 65), SubCallFunc.s)
    }
    ## 'Rh' - relative humidity in %
    if (grepl('Rh', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', -10), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 110), SubCallFunc.s)
    }
    ## 'NEE' - in umol CO2 m-2 s-1 oder g C m-2 day-1
    if (grepl('NEE', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', -50), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 100), SubCallFunc.s)
    }
    ## 'ustar' - in m s-1
    if (grepl('ustar', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', -1), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 50), SubCallFunc.s)
    }
    ## 'E_0' - in degK
    if (grepl('E_0', VarName.V.s[v.i]) )
    {
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0), SubCallFunc.s)
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('>', 600), SubCallFunc.s)
    }
    ## FLUXNET _fqc, 0: data are original, 1: gapfilled high quality, 2: gapfilled medium quality, 3: gapfilled low quality
    if (grepl('_fqc', VarName.V.s[v.i]) && !grepl('_fqcOK', VarName.V.s[v.i], ignore.case = TRUE) ) # 0 is best
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0, '|', '>', 3), SubCallFunc.s)
    ## FLUXNET _fqcOK: 1 = if data are orginal or high quality gapfilled (_fqc was 0 or 1), O = otherwise
    if (grepl('_fqcOK', VarName.V.s[v.i], ignore.case = TRUE) ) # 1 (= 100%) is best
      fCheckOutsideRange(Data.F, VarName.V.s[v.i], c('<', 0, '|', '>', 1), SubCallFunc.s)
  }
  
  ##value<<
  ## Function produces warnings if variable outside range.
}


fCalcLengthOfGaps <- function(
    ##description<<
  ## Calculate length of gaps (with flag NA) in specified vector
  Data.V.n              ##<< Numeric vector with gaps (missing values, NAs)
)
  ##author<<
  ## AMM
  # TEST: Data.V.n <- sDATA$NEE
{
  Data.V.b <- is.na(Data.V.n)
  #Determine cumulative sums of gaps
  CumSum.V.n <- cumsum(Data.V.b)
  #Calculate sum from start of gap
  LenGaps.V.n <- CumSum.V.n-cummax((!Data.V.b) * CumSum.V.n)
  
  LenGaps.V.n
  ##value<<
  ## An integer vector with length of gap since start of gap.
}


fConvertGapsToNA <- function(
    ##description<<
  ## Converts all gap flags to NA
  Data.F                ##<< Data frame to be converted
  , GapFlag.n =-9999.0    ##<< Flag value used for gaps, defaults to -9999.0
  #TEST: Data.F <- ResultData.D
)
  ##author<<
  ## AMM
{
  Gaps.i <- sum(Data.F == GapFlag.n, na.rm = T)
  Data.F[Data.F == GapFlag.n] <- NA_real_
  message('Number of \'', GapFlag.n, '\' convertered to NA: ', Gaps.i)
  
  Data.F
  ##value<<
  ## Data frame with NAs instead of gaps.
}

