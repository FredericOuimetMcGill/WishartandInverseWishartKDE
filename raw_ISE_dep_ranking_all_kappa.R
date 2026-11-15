# Parameters -------------------------------------------------------------
kappa_values <- 4:9           # Cases to process; kappa = 3 is excluded.
tolerance_percentages <- c(2, 5)  # Two positive percentages, in increasing order.
base_directory <- NULL        # NULL = directory containing this R script.
input_filename <- "raw_ISE_dep.csv"
global_csv_filename <- "ranking_counts_all_kappa.csv"
global_tex_filename <- "ranking_counts_all_kappa.tex"
total_csv_filename <- "ranking_counts_total.csv"
total_tex_filename <- "ranking_counts_total.tex"
case_csv_filename <- "ranking_counts.csv"
table_label <- "tab:ranking.counts"
total_table_label <- "tab:ranking.counts.total"

# Put this script in the parent directory containing kappa_4, ..., kappa_9.
# Run it with RStudio's Source button or:
#   Rscript raw_ISE_dep_ranking_all_kappa.R
# Only base R is needed. It writes an eight-method CSV inside each kappa_N,
# a combined CSV/LaTeX table, and a total-count CSV/LaTeX table in the parent
# directory. Minimum, within-2%, and within-5% counts are cumulative; the
# percentages are controlled above. In the manuscript, use:
#   \input{ranking_counts_all_kappa.tex}
#   \input{ranking_counts_total.tex}

local({
  script_directory <- function() {
    frames <- sys.frames()
    for (i in rev(seq_along(frames))) {
      path <- frames[[i]]$ofile
      if (is.character(path) && length(path) == 1L && file.exists(path)) {
        return(dirname(normalizePath(path, winslash = "/")))
      }
    }
    args <- commandArgs(trailingOnly = FALSE)
    arg <- args[startsWith(args, "--file=")]
    if (length(arg)) {
      return(dirname(normalizePath(sub("^--file=", "", arg[1L]),
                                   winslash = "/", mustWork = TRUE)))
    }
    getwd()  # Fallback when the code is pasted directly into the console.
  }

  if (!is.numeric(kappa_values) || !length(kappa_values) ||
      any(!is.finite(kappa_values)) ||
      any(kappa_values <= 0 | kappa_values != floor(kappa_values)) ||
      anyDuplicated(kappa_values)) {
    stop("kappa_values must contain distinct positive integers.")
  }
  if (!is.numeric(tolerance_percentages) || length(tolerance_percentages) != 2L ||
      any(!is.finite(tolerance_percentages)) || any(tolerance_percentages <= 0) ||
      is.unsorted(tolerance_percentages, strictly = TRUE)) {
    stop("tolerance_percentages must contain two distinct positive percentages ",
         "in increasing order, for example c(2, 5).")
  }
  if (is.null(base_directory)) base_directory <- script_directory()
  base_directory <- normalizePath(base_directory, winslash = "/", mustWork = TRUE)
  if (!dir.exists(base_directory)) stop("base_directory is not a directory.")

  # Zero tolerance is the minimum group; the next groups use the percentages.
  decimal_text <- function(x) format(x, scientific = FALSE, trim = TRUE,
                                     digits = 12L, decimal.mark = ".")
  percentages <- vapply(tolerance_percentages, decimal_text, character(1L))
  percentage_tags <- gsub(".", "p", percentages, fixed = TRUE)
  tolerances <- c(0, tolerance_percentages / 100)
  group_tags <- c("minimum", paste0("within_", percentage_tags, "pct"))
  directories <- file.path(base_directory, paste0("kappa_", kappa_values))
  input_files <- file.path(directories, input_filename)
  if (any(!file.exists(input_files))) {
    stop("Missing input file(s):\n", paste(input_files[!file.exists(input_files)],
                                          collapse = "\n"))
  }

  methods <- data.frame(
    Method = c("W_lscv", "W_lcv", "IW_lscv", "IW_lcv",
               "LG_lscv", "LG_lcv", "G_lscv", "G_lcv"),
    kernel = rep(c("Wishart", "inverseWishart", "smlnorm", "smnorm"), each = 2L),
    criterion = rep(c("lscv", "lcv"), times = 4L),
    symbol = rep(c("W", "IW", "LG", "G"), each = 2L),
    stringsAsFactors = FALSE
  )
  displayed <- match(c("W_lscv", "IW_lscv", "W_lcv", "IW_lcv",
                       "LG_lcv", "LG_lscv", "G_lscv", "G_lcv"), methods$Method)
  gaussian_methods <- which(methods$kernel == "smnorm")
  dispersion_methods <- which(methods$kernel != "smnorm")
  method_key <- paste(methods$kernel, methods$criterion, sep = ":")
  statistic_names <- c("Mean", "Median", "SD", "IQR")
  statistic_columns <- unlist(lapply(group_tags, function(tag) {
    paste0(statistic_names, "_", tag)
  }), use.names = FALSE)
  dispersion_columns <- which(rep(statistic_names, length(tolerances)) %in%
                                c("SD", "IQR"))

  # Process each file separately: never pool replications across kappas.
  summarize_case <- function(input_file, kappa) {
    fail <- function(...) stop("kappa = ", kappa, ": ", ..., call. = FALSE)
    raw <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
    required <- c("nobs", "Mmod", "Smod", "kernel", "criterion", "ISE")
    missing_columns <- setdiff(required, names(raw))
    if (length(missing_columns)) {
      fail("Missing columns: ", paste(missing_columns, collapse = ", "))
    }
    if (!nrow(raw)) fail("The input file is empty.")
    for (column in c("nobs", "Mmod", "Smod")) {
      x <- raw[[column]]
      if (!is.numeric(x) || any(!is.finite(x)) ||
          any(x <= 0 | x != floor(x))) {
        fail(column, " must contain finite positive integers.")
      }
    }
    if (!is.numeric(raw$ISE) || any(!is.finite(raw$ISE)) || any(raw$ISE < 0)) {
      fail("ISE must contain only finite, nonnegative numbers.")
    }
    raw$method_id <- match(paste(raw$kernel, raw$criterion, sep = ":"), method_key)
    if (anyNA(raw$method_id)) fail("Unknown or missing kernel/criterion.")

    # Same scaling as the original summary script: 10^5 * true RISE,
    # because the input stores 10^6 * true ISE. Scaling does not affect counts.
    raw$RISE <- 100 * sqrt(raw$ISE)
    configurations <- split(
      seq_len(nrow(raw)),
      interaction(raw$nobs, raw$Mmod, raw$Smod, drop = TRUE, lex.order = TRUE)
    )
    design <- lapply(raw[c("Mmod", "Smod", "nobs")], function(x) {
      sort(unique(as.numeric(x)))
    })
    if (length(configurations) != prod(lengths(design))) {
      fail("Incomplete grid of Mmod, Smod, and nobs combinations.")
    }
    counts <- matrix(0L, nrow = nrow(methods), ncol = length(statistic_columns),
                     dimnames = list(methods$Method, statistic_columns))

    for (indices in configurations) {
      data <- raw[indices, ]
      replications <- tabulate(data$method_id, nbins = nrow(methods))
      config_label <- paste0("(M_", data$Mmod[1L], ", Sigma_", data$Smod[1L],
                             ", n = ", data$nobs[1L], ")")
      if (any(replications < 2L)) {
        fail("At least two valid replications per method are required in ",
             config_label, ". Missing/insufficient methods: ",
             paste(methods$Method[replications < 2L], collapse = ", "))
      }
      if (length(unique(replications)) > 1L) {
        warning("kappa = ", kappa, ": unequal replication counts in ", config_label)
      }

      # Summarize replicate-level RISE, without rounding.
      # sd() is sample SD; IQR(type = 7) uses R's default quantiles.
      statistics <- t(vapply(seq_len(nrow(methods)), function(j) {
        x <- data$RISE[data$method_id == j]
        c(Mean = mean(x), Median = median(x), SD = sd(x), IQR = IQR(x, type = 7))
      }, FUN.VALUE = setNames(numeric(4L), statistic_names)))
      if (any(!is.finite(statistics))) fail("Nonfinite summary in ", config_label)

      # Mean and median: minimum among ALL EIGHT methods.
      # SD and IQR: minimum among ONLY THE SIX non-Gaussian methods.
      minima <- vapply(seq_along(statistic_names), function(j) {
        eligible <- if (j <= 2L) seq_len(nrow(methods)) else dispersion_methods
        min(statistics[eligible, j])
      }, FUN.VALUE = numeric(1L))
      # Inclusive threshold; a zero minimum allows only zero values.
      # At tolerance zero, all exact ties for the minimum are counted.
      for (g in seq_along(tolerances)) {
        columns <- (4L * (g - 1L) + 1L):(4L * g)
        qualifies <- sweep(statistics, 2L, (1 + tolerances[g]) * minima, FUN = "<=")
        qualifies[gaussian_methods, 3:4] <- FALSE
        counts[, columns] <- counts[, columns] + qualifies
      }
    }

    # Excluded dispersion counts are NA internally, never zero.
    # Both CSV and LaTeX exports render them as exactly three hyphens.
    counts[gaussian_methods, dispersion_columns] <- NA_integer_
    result <- data.frame(kappa = kappa, Method = methods$Method[displayed],
                         counts[displayed, , drop = FALSE],
                         row.names = NULL, check.names = FALSE)
    names(result)[-(1:2)] <- statistic_columns
    list(table = result, design = design, n_configurations = length(configurations))
  }

  tables <- vector("list", length(kappa_values))
  reference_design <- NULL
  for (i in seq_along(kappa_values)) {
    message("Processing kappa = ", kappa_values[i], " ...")
    case <- summarize_case(input_files[i], kappa_values[i])
    if (is.null(reference_design)) {
      reference_design <- case$design
      n_configurations <- case$n_configurations
    } else if (!identical(case$design, reference_design)) {
      stop("Configuration grids differ across kappas; a common table caption ",
           "would be misleading. Check ", input_files[i])
    }
    tables[[i]] <- case$table
  }
  combined <- do.call(rbind, tables)
  rownames(combined) <- NULL
  total_counts <- Reduce("+", lapply(tables, function(x) {
    as.matrix(x[statistic_columns])
  }))
  totals <- data.frame(Method = methods$Method[displayed], total_counts,
                       row.names = NULL, check.names = FALSE)

  # Construct complete LaTeX tables, including grouped headers and notes.
  tex_set <- function(x) paste0("\\{", paste(x, collapse = ","), "\\}")
  M_levels <- reference_design$Mmod
  S_levels <- reference_design$Smod
  n_levels <- reference_design$nobs
  pair_indices <- if (identical(M_levels, S_levels)) {
    paste0("$i,j\\in", tex_set(M_levels), "$")
  } else {
    paste0("$i\\in", tex_set(M_levels), "$ and $j\\in", tex_set(S_levels), "$")
  }
  method_tex <- paste0("$\\mathrm{", methods$symbol[displayed],
                       "}_{\\mathrm{", methods$criterion[displayed], "}}$")
  names(method_tex) <- methods$Method[displayed]
  row_end <- " \\\\"
  group_headers <- c(
    "\\multicolumn{4}{c}{Minimum}",
    vapply(percentages, function(p) {
      paste0("\\multicolumn{4}{c}{Within $", p,
             "\\%$ of the minimum}")
    }, character(1L))
  )
  note <- paste0(
    "\\vspace{-4mm}$^{*}$SD and IQR comparisons exclude $\\mathrm{G}_{\\mathrm{lscv}}$ and $\\mathrm{G}_{\\mathrm{lcv}}$ since they always have the worst mean and median.\\vspace{-7mm}"
  )
  design_description <- paste0(
    "the $", n_configurations, "$ configurations combine the $",
    length(M_levels) * length(S_levels), "$ matrix pairs $(M_i,\\Sigma_j)$, ",
    pair_indices, ", with the $", length(n_levels), "$ sample sizes $n\\in",
    tex_set(n_levels), "$"
  )
  main_caption <- paste0(
    "Numbers of configurations in which the Wishart (W), inverse-Wishart (IW), ",
    "log-Gaussian (LG), and Gaussian (G) methods attain the minimum or fall ",
    "within the indicated percentage of the minimum for each RISE summary statistic, for $\\kappa\\in",
    tex_set(kappa_values), "$. The summaries are the mean, median, standard ",
    "deviation (SD), and interquartile range (IQR) of $1024$ RISE values for ",
    "each method, model, and sample size. ",
    "For each $\\kappa$, ", design_description, ". Here $M_i$ and $\\Sigma_j$ ",
    "are the autoregressive coefficient and innovation covariance matrices, ",
    "respectively. Dark and light gray cells indicate the highest and ",
    "second-highest distinct counts, respectively, within each $\\kappa$ ",
    "and column."
  )
  total_caption <- paste0(
    "Totals of the counts in Table~\\ref{", table_label,
    "} over $\\kappa\\in", tex_set(kappa_values), "$."
  )

  # Match Table~\ref{tab:summary}: black!30 and black!10.
  # Here larger counts rank first. All ties at a rank share its shade.
  shade_counts <- function(values) {
    shades <- rep("", length(values))
    ranked_values <- sort(unique(values[is.finite(values)]), decreasing = TRUE)
    if (length(ranked_values) >= 1L) {
      shades[which(values == ranked_values[1L])] <- "\\cellcolor{black!30}"
    }
    if (length(ranked_values) >= 2L) {
      shades[which(values == ranked_values[2L])] <- "\\cellcolor{black!10}"
    }
    shades
  }

  render_table <- function(data, by_kappa, caption, label) {
    n_leading <- if (by_kappa) 2L else 1L
    stub <- if (by_kappa) c("$\\kappa$", "Method") else "Method"
    rule <- paste(vapply(seq_along(group_headers), function(g) {
      left <- n_leading + 4L * (g - 1L) + 1L
      paste0("\\cmidrule(lr){", left, "-", left + 3L, "}")
    }, character(1L)), collapse = "")
    second_header <- paste(c(rep("", n_leading),
                             rep(c("Mean", "Median", "SD$^{*}$", "IQR$^{*}$"),
                                 length(group_headers))), collapse = " & ")
    # Rank each count column independently within its kappa block.
    # The total-count table has one block containing all methods.
    shading <- matrix("", nrow = nrow(data), ncol = length(statistic_columns))
    blocks <- if (by_kappa) {
      split(seq_len(nrow(data)), data$kappa)
    } else {
      list(seq_len(nrow(data)))
    }
    for (rows in blocks) {
      for (j in seq_along(statistic_columns)) {
        shading[rows, j] <- shade_counts(data[[statistic_columns[j]]][rows])
      }
    }
    body <- character()
    for (i in seq_len(nrow(data))) {
      new_kappa <- by_kappa && (i == 1L || data$kappa[i] != data$kappa[i - 1L])
      if (new_kappa && i > 1L) body <- c(body, "\\midrule")
      values <- as.numeric(unlist(data[i, statistic_columns], use.names = FALSE))
      cells <- vapply(seq_along(values), function(j) {
        if (is.na(values[j])) "\\texttt{---}" else
          paste0(shading[i, j], "$", values[j], "$")
      }, character(1L))
      leading <- unname(method_tex[data$Method[i]])
      if (by_kappa) {
        kappa_cell <- if (new_kappa) paste0("$", data$kappa[i], "$") else ""
        leading <- c(kappa_cell, leading)
      }
      body <- c(body, paste0(paste(c(leading, cells), collapse = " & "), row_end))
    }
    title <- c(paste0("\\caption{", caption, "}"), paste0("\\label{", label, "}"))
    c(
      "\\begin{table}[!htp]",
      title,
      "\\centering",
      "\\fontsize{7.9}{7}\\selectfont",
      "\\setlength{\\tabcolsep}{4pt}",
      "\\renewcommand{\\arraystretch}{0.15}",
      paste0("\\begin{tabular}{", if (by_kappa) "cl" else "l",
             paste(rep("c", length(statistic_columns)), collapse = ""), "}"),
      "\\toprule",
      paste0(paste(c(stub, group_headers), collapse = " & "), row_end),
      rule,
      paste0(second_header, row_end),
      "\\midrule",
      body,
      "\\bottomrule",
      "\\end{tabular}",
      if (by_kappa) c(
        "",
        "\\smallskip",
        "\\begin{minipage}{0.98\\textwidth}",
        "\\fontsize{7.9}{7}\\selectfont",
        note,
        "\\end{minipage}"
      ) else character(),
      "\\end{table}"
    )
  }
  tex <- render_table(combined, TRUE, main_caption, table_label)
  total_tex_lines <- render_table(totals, FALSE, total_caption, total_table_label)

  # Write outputs after all selected cases have been computed and checked.
  for (i in seq_along(kappa_values)) {
    write.csv(tables[[i]], file.path(directories[i], case_csv_filename),
              row.names = FALSE, na = "---", fileEncoding = "UTF-8")
  }
  global_csv <- file.path(base_directory, global_csv_filename)
  global_tex <- file.path(base_directory, global_tex_filename)
  total_csv <- file.path(base_directory, total_csv_filename)
  total_tex <- file.path(base_directory, total_tex_filename)
  write.csv(combined, global_csv, row.names = FALSE, na = "---", fileEncoding = "UTF-8")
  writeLines(tex, global_tex, useBytes = TRUE)
  write.csv(totals, total_csv, row.names = FALSE, na = "---", fileEncoding = "UTF-8")
  writeLines(total_tex_lines, total_tex, useBytes = TRUE)
  message("Created per-case CSV files for kappa = ", paste(kappa_values, collapse = ", "))
  message("Created: ", global_csv, "\nCreated: ", global_tex)
  message("Created: ", total_csv, "\nCreated: ", total_tex)
})
