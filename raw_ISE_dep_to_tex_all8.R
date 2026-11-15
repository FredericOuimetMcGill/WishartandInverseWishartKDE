## Build one CSV and LaTeX summary table per kappa directory.

# Parameters -------------------------------------------------------------
kappa_values <- 4:9                    # Cases to process; kappa = 3 is excluded.
base_directory <- NULL                 # NULL = directory containing this script.
input_filename <- "raw_ISE_dep.csv"     # Raw ISE stores 10^6 * true ISE.
output_csv_filename <- "table_summary_dep.csv"
output_tex_filename <- "table_summary_dep.tex"
digits_RISE <- 0                       # Decimal places in the LaTeX table.
n_levels <- c(100, 200, 300)            # Sample sizes in the table.
table_label <- "tab:summary"            # Existing reference in each case's manuscript.

# Place this script in the parent of kappa_4, ..., kappa_9 and use Source
# in RStudio or run it with Rscript. Only base R is required.
# Each selected directory receives table_summary_dep.csv and
# table_summary_dep.tex, computed ONLY from that directory's raw_ISE_dep.csv.
# CSV statistics are unrounded and use the same 10^5 * true RISE units as
# the LaTeX table. The caption contains the corresponding kappa value.
# Example manuscript input: \input{kappa_4/table_summary_dep.tex}
# By default the tables retain tab:summary for separate manuscripts.
# If including multiple cases in one manuscript, set table_label to
# "tab:summary.kappa_{kappa}"; {kappa} is replaced by the current value.

local({
  script_directory <- function() {
    frames <- sys.frames()
    for (i in rev(seq_along(frames))) {
      path <- frames[[i]]$ofile
      if (is.character(path) && length(path) == 1L && file.exists(path)) {
        return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
      }
    }
    args <- commandArgs(trailingOnly = FALSE)
    arg <- args[startsWith(args, "--file=")]
    if (length(arg)) {
      return(dirname(normalizePath(sub("^--file=", "", arg[1L]),
                                   winslash = "/", mustWork = TRUE)))
    }
    getwd()  # Fallback when pasting the script into the console.
  }

  if (!is.numeric(kappa_values) || !length(kappa_values) ||
      any(!is.finite(kappa_values)) ||
      any(kappa_values <= 0 | kappa_values != floor(kappa_values)) ||
      anyDuplicated(kappa_values)) {
    stop("kappa_values must contain distinct positive integers.")
  }
  if (is.null(base_directory)) base_directory <- script_directory()
  base_directory <- normalizePath(base_directory, winslash = "/", mustWork = TRUE)
  if (!dir.exists(base_directory)) stop("base_directory is not a directory.")

  directories <- file.path(base_directory, paste0("kappa_", kappa_values))
  input_files <- file.path(directories, input_filename)
  if (any(!file.exists(input_files))) {
    stop("Missing input file(s):\n",
         paste(input_files[!file.exists(input_files)], collapse = "\n"))
  }
  if (input_filename %in% c(output_csv_filename, output_tex_filename) ||
      identical(output_csv_filename, output_tex_filename)) {
    stop("Input and output filenames must be different.")
  }

  write_case_tables <- function(kappa, case_directory) {
    input_file <- file.path(case_directory, input_filename)
    output_csv <- file.path(case_directory, output_csv_filename)
    output_tex <- file.path(case_directory, output_tex_filename)
    table_label <- gsub("{kappa}", as.character(kappa), table_label, fixed = TRUE)

    # Read raw data ----------------------------------------------------------

    df_raw <- read.csv(input_file, stringsAsFactors = FALSE)
    required_columns <- c("nobs", "Mmod", "Smod", "kernel", "criterion", "ISE")
    missing_columns <- setdiff(required_columns, names(df_raw))
    if (length(missing_columns)) {
      stop("Missing columns in ", input_file, ": ", paste(missing_columns, collapse = ", "))
    }
    if (!is.numeric(df_raw$ISE)) stop("ISE must be numeric in ", input_file)

    # Match the kernel names used by the simulation and timing scripts.
    kernels <- c("Wishart", "inverseWishart", "smlnorm", "smnorm")
    criteria <- c("lscv", "lcv")
    if (anyNA(df_raw$kernel) || any(!df_raw$kernel %in% kernels)) {
      stop("Unknown kernel in raw_ISE_dep.csv. Expected Wishart, inverseWishart, smlnorm, or smnorm.")
    }
    if (anyNA(df_raw$criterion) || any(!df_raw$criterion %in% criteria)) {
      stop("Unknown criterion in raw_ISE_dep.csv. Expected lscv or lcv.")
    }
    df_raw <- df_raw[is.finite(df_raw$ISE) & df_raw$ISE >= 0, ]
    if (nrow(df_raw) == 0L) stop("No finite nonnegative ISE results in raw_ISE_dep.csv.")

    expected_groups <- expand.grid(
      nobs = n_levels, Mmod = 1:3, Smod = 1:3, kernel = kernels, criterion = criteria,
      stringsAsFactors = FALSE
    )
    group_columns <- names(expected_groups)
    observed_groups <- unique(df_raw[group_columns])
    observed_groups$available <- TRUE
    coverage <- merge(expected_groups, observed_groups, by = group_columns, all.x = TRUE)
    if (anyNA(coverage$available)) {
      warning("Missing ISE results for ", sum(is.na(coverage$available)),
              " model/sample-size/method groups. Missing entries are shown as ***; incomplete comparisons are not shaded.")
    }

    # Compute RISE metric: 10^5 * sqrt(ISE_true) = 100 * sqrt(raw_ISE) because raw_ISE = ISE_true * 10^6
    # raw ISE column is assumed to be named 'ISE'
    df_raw$RISE <- 100 * sqrt(df_raw$ISE)

    # Aggregate to get mean, median, standard deviation, and IQR of RISE per (nobs, Mmod, Smod, kernel, criterion)

    summ_mean <- aggregate(
      RISE ~ nobs + Mmod + Smod + kernel + criterion,
      data = df_raw,
      FUN  = mean,
      na.rm = TRUE
    )
    names(summ_mean)[names(summ_mean) == "RISE"] <- "mean_RISE"

    summ_median <- aggregate(
      RISE ~ nobs + Mmod + Smod + kernel + criterion,
      data = df_raw,
      FUN  = median,
      na.rm = TRUE
    )
    names(summ_median)[names(summ_median) == "RISE"] <- "median_RISE"

    # Sample standard deviation of RISE across the available replications.
    summ_SD <- aggregate(
      RISE ~ nobs + Mmod + Smod + kernel + criterion,
      data = df_raw,
      FUN  = sd,
      na.rm = TRUE
    )
    names(summ_SD)[names(summ_SD) == "RISE"] <- "SD_RISE"

    summ_IQR <- aggregate(
      RISE ~ nobs + Mmod + Smod + kernel + criterion,
      data = df_raw,
      FUN  = IQR,
      na.rm = TRUE
    )
    names(summ_IQR)[names(summ_IQR) == "RISE"] <- "IQR_RISE"

    # Merge summaries into one data.frame (this will be our 'df' below)
    df <- merge(
      summ_mean,
      summ_median,
      by = c("nobs", "Mmod", "Smod", "kernel", "criterion"),
      all = TRUE
    )
    df <- merge(
      df,
      summ_SD,
      by = c("nobs", "Mmod", "Smod", "kernel", "criterion"),
      all = TRUE
    )
    df <- merge(
      df,
      summ_IQR,
      by = c("nobs", "Mmod", "Smod", "kernel", "criterion"),
      all = TRUE
    )

    # Helper: get formatted statistic or *** if missing; shade first/second smallest ----

    get_stat <- function(df, n, M, S, kernel, crit, stat, digits = 1) {
      row <- df[
        df$nobs == n &
          df$Mmod == M &
          df$Smod == S &
          df$kernel == kernel &
          df$criterion == crit,
      ]
      
      # Missing entry -> ***
      if (nrow(row) == 0L || !is.finite(row[[stat]][1])) {
        return("***")
      }
      
      val <- row[[stat]][1]
      # Rounded to nearest integer if digits = 0
      fmt <- formatC(val, format = "f", digits = digits)
      
      vals <- vapply(
        methods,
        function(m_info) {
          rank_row <- df[
            df$nobs == n &
              df$Mmod == M &
              df$Smod == S &
              df$kernel == m_info$kernel &
              df$criterion == m_info$criterion,
          ]
          if (nrow(rank_row) == 0L || !is.finite(rank_row[[stat]][1])) {
            return(NA_real_)
          }
          rank_row[[stat]][1]
        },
        FUN.VALUE = numeric(1)
      )
      # Rank only complete eight-method comparisons.
      if (any(!is.finite(vals))) return(fmt)
      vals <- sort(unique(vals))
      
      if (length(vals) >= 1L && abs(val - vals[1]) < 1e-8) {
        return(paste0("\\cellcolor{black!30}", fmt))
      }
      if (length(vals) >= 2L && abs(val - vals[2]) < 1e-8) {
        return(paste0("\\cellcolor{black!10}", fmt))
      }
      
      fmt
    }

    # Methods mapping (kernel + criterion -> LaTeX label + CSV values) ------

    methods <- list(
      list(tex = "\\mathrm{W}_{\\mathrm{lscv}}", kernel = "Wishart", criterion = "lscv"),
      list(tex = "\\mathrm{W}_{\\mathrm{lcv}}",  kernel = "Wishart", criterion = "lcv"),
      list(tex = "\\mathrm{IW}_{\\mathrm{lscv}}", kernel = "inverseWishart", criterion = "lscv"),
      list(tex = "\\mathrm{IW}_{\\mathrm{lcv}}",  kernel = "inverseWishart", criterion = "lcv"),
      list(tex = "\\mathrm{LG}_{\\mathrm{lscv}}", kernel = "smlnorm", criterion = "lscv"),
      list(tex = "\\mathrm{LG}_{\\mathrm{lcv}}",  kernel = "smlnorm", criterion = "lcv"),
      list(tex = "\\mathrm{G}_{\\mathrm{lscv}}",  kernel = "smnorm", criterion = "lscv"),
      list(tex = "\\mathrm{G}_{\\mathrm{lcv}}",   kernel = "smnorm", criterion = "lcv")
    )

    # Write unrounded CSV summaries in the same units as the LaTeX table.
    # Keep all expected groups; unavailable statistics are written as NA.
    csv_summary <- merge(expected_groups, df, by = group_columns, all.x = TRUE)
    method_keys <- vapply(methods, function(m) paste(m$kernel, m$criterion), character(1))
    method_order <- match(paste(csv_summary$kernel, csv_summary$criterion), method_keys)
    csv_summary <- csv_summary[order(csv_summary$Mmod, csv_summary$Smod,
                                     csv_summary$nobs, method_order), ]
    csv_summary <- data.frame(kappa = kappa, csv_summary, row.names = NULL)
    write.csv(csv_summary, output_csv, row.names = FALSE, na = "NA")

    # Open .tex file ---------------------------------------------------------

    con <- file(output_tex, open = "wt")
    on.exit(close(con), add = TRUE)

    cat("%%% TABLE 1\n", file = con)
    cat("\\begin{table}[!ht]\n", file = con)
    cat(paste0("\\caption{Comparison of the Wishart (W), inverse-Wishart (IW), log-Gaussian (LG), and Gaussian (G) methods in terms of the mean, median, standard deviation (SD), and interquartile range (IQR) of $1024$ RISE values ($\\times 10^5$) with $\\kappa = ", kappa, "$ for different configurations of autoregressive coefficient matrices ($M_i, \\, i\\in \\{1,2,3\\}$), innovation covariance matrices ($\\Sigma_j, \\, j\\in \\{1,2,3\\}$), and sample sizes $n \\in \\{100, 200, 300\\}$. Dark gray cells indicate the smallest value among the 8 competitors for the corresponding model, sample size, and statistic; light gray cells indicate the second smallest.}\n"), file = con)
    cat(paste0("\\label{", table_label, "}\n\n"), file = con)
    cat("\\fontsize{7.9}{7}\\selectfont", file = con)
    cat("\\renewcommand{\\arraystretch}{0}\n", file = con)
    cat("\\setlength{\\tabcolsep}{1.2pt}\n", file = con)
    cat("\\setlength{\\aboverulesep}{0.35ex}\n", file = con)
    cat("\\setlength{\\belowrulesep}{0.35ex}\n", file = con)
    cat("\\centering\n", file = con)
    cat("{\n", file = con)
    cat("\\begin{tabular}{@{}cc|ccc|ccc|ccc|ccc@{}}\n", file = con)
    cat("\\toprule\n", file = con)
    cat("\\multicolumn{2}{c|}{} & $n = 100$ & $n = 200$ & $n = 300$ & $n = 100$ & $n = 200$ & $n = 300$ & $n = 100$ & $n = 200$ & $n = 300$ & $n = 100$ & $n = 200$ & $n = 300$ \\\\\n", file = con)
    cat("\\cmidrule(r){1-2}\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}\\cmidrule(l){12-14}\n", file = con)
    cat("Model & Method & \\multicolumn{3}{c|}{Mean RISE} & \\multicolumn{3}{c|}{Median RISE} & \\multicolumn{3}{c|}{SD RISE} & \\multicolumn{3}{c}{IQR RISE} \\\\\n", file = con)
    cat("\\midrule\n", file = con)

    # Write body: loop over (M_i, Sigma_j) blocks ----------------------------

    for (M in 1:3) {
      for (S in 1:3) {
        # Eight methods for this (M,S) block
        for (k in seq_along(methods)) {
          m_info <- methods[[k]]
          model_str  <- if (k == 1L) sprintf("$(M_%d,\\Sigma_%d)$", M, S) else ""
          method_tex <- sprintf("$%s$", m_info$tex)
          
          # Mean RISE for n = 100, 200, 300
          mean_vals <- vapply(
            n_levels,
            function(n) get_stat(
              df,
              n = n,
              M = M,
              S = S,
              kernel = m_info$kernel,
              crit = m_info$criterion,
              stat = "mean_RISE",
              digits = digits_RISE
            ),
            FUN.VALUE = character(1)
          )
          
          # Median RISE for n = 100, 200, 300
          med_vals <- vapply(
            n_levels,
            function(n) get_stat(
              df,
              n = n,
              M = M,
              S = S,
              kernel = m_info$kernel,
              crit = m_info$criterion,
              stat = "median_RISE",
              digits = digits_RISE
            ),
            FUN.VALUE = character(1)
          )
          
          # Standard deviation of RISE for n = 100, 200, 300
          sd_vals <- vapply(
            n_levels,
            function(n) get_stat(
              df,
              n = n,
              M = M,
              S = S,
              kernel = m_info$kernel,
              crit = m_info$criterion,
              stat = "SD_RISE",
              digits = digits_RISE
            ),
            FUN.VALUE = character(1)
          )
          
          # IQR RISE for n = 100, 200, 300
          iqr_vals <- vapply(
            n_levels,
            function(n) get_stat(
              df,
              n = n,
              M = M,
              S = S,
              kernel = m_info$kernel,
              crit = m_info$criterion,
              stat = "IQR_RISE",
              digits = digits_RISE
            ),
            FUN.VALUE = character(1)
          )
          
          row_tex <- sprintf(
            "%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
            model_str,
            method_tex,
            mean_vals[1], mean_vals[2], mean_vals[3],
            med_vals[1], med_vals[2], med_vals[3],
            sd_vals[1], sd_vals[2], sd_vals[3],
            iqr_vals[1], iqr_vals[2], iqr_vals[3]
          )
          cat(row_tex, "\n", file = con)
        }
        
        # Add \midrule between blocks, but not after the last one
        if (!(M == 3 && S == 3)) {
          cat("\\midrule\n", file = con)
        }
      }
    }

    # Close table -----------------------------------------------------------

    cat("\\bottomrule\n", file = con)
    cat("\\end{tabular}}\n", file = con)
    cat("\\end{table}\n", file = con)

    message("kappa = ", kappa, ": wrote ", output_csv, " and ", output_tex)
    invisible(csv_summary)
  }

  for (i in seq_along(kappa_values)) {
    message("Processing kappa = ", kappa_values[i], " ...")
    tryCatch(
      write_case_tables(kappa_values[i], directories[i]),
      error = function(e) {
        stop("kappa = ", kappa_values[i], ": ", conditionMessage(e), call. = FALSE)
      }
    )
  }
})
