# Parameters -------------------------------------------------------------
kappa_values <- 4:9                    # Cases to process; kappa = 3 is excluded.
base_directory <- NULL                 # NULL = directory containing this script.
input_filename <- "raw_ISE_dep.csv"
output_filename <- "RISE_boxplots.png"  # Written inside each selected kappa directory.
text_scale <- 1.5                      # Scale all existing text sizes.

# Place this script in the parent of kappa_4, ..., kappa_9 and use Source
# in RStudio or run it with Rscript. Each selected directory must contain
# raw_ISE_dep.csv and receives its own RISE_boxplots.png (12 x 16 in, 300 dpi).
# Each figure uses only that directory's data. Existing figures are replaced.
# If needed, install the plotting packages once with:
# install.packages(c("ggplot2", "dplyr", "ggh4x"))

library(ggplot2)
library(dplyr)
library(ggh4x)

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
  if (identical(input_filename, output_filename)) {
    stop("Input and output filenames must be different.")
  }

  plot_case <- function(kappa, case_directory) {
    raw_file <- file.path(case_directory, input_filename)
    output_file <- file.path(case_directory, output_filename)

    dat <- read.csv(raw_file, stringsAsFactors = FALSE)
    required_columns <- c("ISE", "kernel", "criterion", "Mmod", "Smod", "nobs")
    missing_columns <- setdiff(required_columns, names(dat))
    if (length(missing_columns)) {
      stop("Missing columns in ", raw_file, ": ", paste(missing_columns, collapse = ", "))
    }
    if (!is.numeric(dat$ISE)) stop("ISE must be numeric in ", raw_file)

    # The CSV must contain both bandwidth criteria for all four kernels.
    kernel_codes <- c(Wishart = "W", inverseWishart = "IW", smlnorm = "LG", smnorm = "G")
    method_order <- c("W_lscv", "W_lcv", "IW_lscv", "IW_lcv", "LG_lscv", "LG_lcv", "G_lscv", "G_lcv")
    if (anyNA(dat$kernel) || any(!dat$kernel %in% names(kernel_codes))) {
      stop("Unknown kernel in raw_ISE_dep.csv. Expected Wishart, inverseWishart, smlnorm, or smnorm.")
    }
    if (anyNA(dat$criterion) || any(!dat$criterion %in% c("lscv", "lcv"))) {
      stop("Unknown criterion in raw_ISE_dep.csv. Expected lscv or lcv.")
    }

    ## --- Transform ISE -> RISE = sqrt(ISE / 1e6) --------------------------
    dat <- dat %>%
      filter(is.finite(ISE), ISE >= 0) %>%
      mutate(
        RISE = sqrt(ISE / 1e6),
        kernel_short = unname(kernel_codes[kernel]),
        method = factor(
          paste(kernel_short, criterion, sep = "_"),
          levels = method_order
        ),
        Mmod_f = factor(Mmod, levels = 1:3, labels = paste0("M", 1:3)),
        Smod_f = factor(Smod, levels = 1:3, labels = paste0("Sigma[", 1:3, "]")),
        nobs_f = factor(nobs, levels = sort(unique(nobs)))
      )

    # Do not overwrite the figure with an incomplete comparison.
    expected_groups <- expand.grid(
      Mmod = 1:3, Smod = 1:3, nobs = c(100, 200, 300), method = method_order,
      stringsAsFactors = FALSE
    )
    observed_groups <- unique(as.data.frame(dat[c("Mmod", "Smod", "nobs", "method")]))
    missing_groups <- dplyr::anti_join(expected_groups, observed_groups, by = c("Mmod", "Smod", "nobs", "method"))
    if (nrow(missing_groups) > 0L) {
      stop("Missing finite ISE results for ", nrow(missing_groups),
           " model/sample-size/method groups. Regenerate raw_ISE_dep.csv including inverseWishart with lscv and lcv before plotting.")
    }

    ## --- Grouped colour palette for the 8 methods ------------------------
    method_cols <- c(
      "W_lscv"  = "#0072B2", # blue
      "W_lcv"   = "#009E73", # bluish green
      "IW_lscv" = "#56B4E9", # sky blue
      "IW_lcv"  = "#D55E00", # vermillion
      "LG_lscv" = "#6A3D9A", # purple
      "LG_lcv"  = "#CC6677", # muted red
      "G_lscv"  = "#E69F00", # orange
      "G_lcv"   = "#B79F00"  # golden yellow
    )

    ## (1) labels with lscv / lcv as subscripts -----------------------------
    method_labels <- c(
      expression(W[lscv]),
      expression(W[lcv]),
      expression(IW[lscv]),
      expression(IW[lcv]),
      expression(LG[lscv]),
      expression(LG[lcv]),
      expression(G[lscv]),
      expression(G[lcv])
    )

    ## --- Alpha for sample sizes -------------------------------------------
    n_levels <- levels(dat$nobs_f)
    n_sizes  <- length(n_levels)

    # smallest n = opaque (alpha = 1), largest n = pretty pale
    alpha_vals <- seq(1.0, 0.15, length.out = n_sizes)
    names(alpha_vals) <- n_levels

    # greys for sample-size legend outlines (100 = black, larger n = lighter grey)
    if (n_sizes > 1) {
      grey_vals <- c(
        "#000000",
        gray.colors(n_sizes - 1, start = 0.3, end = 0.5)
      )
    } else {
      grey_vals <- "#000000"
    }
    names(grey_vals) <- n_levels

    ## (3) row-wise upper bounds based on whiskers (not outliers) ----------

    whisker_by_row <- dat %>%
      group_by(Mmod_f, Smod_f, method, nobs_f) %>%
      summarise(
        Q1 = quantile(RISE, 0.25, na.rm = TRUE),
        Q3 = quantile(RISE, 0.75, na.rm = TRUE),
        IQR = Q3 - Q1,
        upper_whisker = max(RISE[RISE <= Q3 + 1.5 * IQR], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      group_by(Mmod_f) %>%
      summarise(
        row_upper = max(upper_whisker, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(row_upper = row_upper * 1.02)  # tiny headroom above top whisker

    ylim_M1 <- c(0, whisker_by_row$row_upper[whisker_by_row$Mmod_f == "M1"])
    ylim_M2 <- c(0, whisker_by_row$row_upper[whisker_by_row$Mmod_f == "M2"])
    ylim_M3 <- c(0, whisker_by_row$row_upper[whisker_by_row$Mmod_f == "M3"])

    ## --- Plot -------------------------------------------------------------
    p_RISE <- ggplot(
      dat,
      aes(
        x     = method,
        y     = RISE,
        fill  = method,                # 8 methods -> 8 colours inside
        alpha = nobs_f,                # sample size -> transparency (fill + outline)
        group = interaction(method, nobs_f)
      )
    ) +
      geom_boxplot(
        position      = position_dodge(width = 0.7),
        width         = 0.6,
        colour        = "black",        # always black outline (fades via alpha)
        size          = 0.3,
        outlier.shape = NA              # (3) remove outliers: only box + whiskers
      ) +
      facet_grid(
        Mmod_f ~ Smod_f,
        labeller = label_parsed,
        scales   = "free_y"
      ) +
      # x-axis labels with subscripts
      scale_x_discrete(limits = method_order, labels = method_labels, drop = FALSE) +
      # coloured method legend, with subscripts
      scale_fill_manual(
        name   = "Method",
        values = method_cols,
        breaks = method_order,
        labels = method_labels,
        drop   = FALSE
      ) +
      scale_alpha_manual(
        name   = "Sample size",
        values = alpha_vals
      ) +
      # we don't map colour to a scale (constant black in geom), so no colour legend
      labs(
        x     = NULL,
        y     = "RISE",
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0.02, 0.06))   # small margin above/below data
      ) +
      theme_bw(base_size = 11) +
      theme(
        aspect.ratio     = 5/4,   
        strip.background = element_rect(colour = "black", fill = "grey95"),
        strip.text       = element_text(
          face = "bold",
          size = 11 * text_scale
        ),
        axis.text.x      = element_text(
          angle  = 0,
          hjust  = 0.5,
          vjust  = 0.5,
          colour = "black",
          size   = 9 * text_scale
        ),
        axis.text.y      = element_text(
          colour = "black",
          size   = 11 * text_scale
        ),
        axis.title.y     = element_blank(),
        legend.position  = "bottom",
        legend.box       = "vertical",
        legend.title     = element_text(
          face = "bold",
          size = 11 * text_scale
        ),
        legend.text      = element_text(
          size = 11 * text_scale
        ),
        ## (2) bring the Sample size legend closer to the Method legend
        legend.spacing.y  = grid::unit(0.2, "lines"),
        plot.margin       = margin(t = 3, r = 0.1, b = 3, l = 0.1, unit = "mm"),
        legend.margin     = margin(t = 0, r = 0, b = 2, l = 0, unit = "mm"),
        legend.box.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "mm")
      ) +
      guides(
        # coloured method legend: solid colours with no transparency
        fill = guide_legend(
          order = 1,
          nrow = 1,
          byrow = TRUE,
          override.aes = list(alpha = 1)
        ),
        # sample-size legend: white boxes with black/grey outlines, alpha still mapped
        alpha = guide_legend(
          order = 2,
          override.aes = list(
            fill   = "white",
            colour = grey_vals
          )
        )
      ) +
      ## (3) per-row y-limits just above the max whisker on that row
      ggh4x::facetted_pos_scales(
        y = list(
          Mmod_f == "M1" ~ scale_y_continuous(
            limits = ylim_M1,
            expand = expansion(mult = c(0.02, 0.06))
          ),
          Mmod_f == "M2" ~ scale_y_continuous(
            limits = ylim_M2,
            expand = expansion(mult = c(0.02, 0.06))
          ),
          Mmod_f == "M3" ~ scale_y_continuous(
            limits = ylim_M3,
            expand = expansion(mult = c(0.02, 0.06))
          )
        )
      )

    ## --- Save -------------------------------------------------------------
    ggsave(
      filename = output_file,
      plot     = p_RISE,
      width    = 12,            # any width you like
      height   = 12 * 4/3,      # for panel ratio 4/3
      dpi      = 300
    )

    message("kappa = ", kappa, ": wrote ", output_file)
    invisible(output_file)
  }

  for (i in seq_along(kappa_values)) {
    message("Processing kappa = ", kappa_values[i], " ...")
    tryCatch(
      plot_case(kappa_values[i], directories[i]),
      error = function(e) {
        stop("kappa = ", kappa_values[i], ": ", conditionMessage(e), call. = FALSE)
      }
    )
  }
})
