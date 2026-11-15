########################################
## Main hyperparameter                ##
########################################

kappa <- 4L  # WAR(1) degrees of freedom; change this value for another case.

if (!is.numeric(kappa) || length(kappa) != 1L || !is.finite(kappa) ||
    kappa < 2 || kappa != floor(kappa) || kappa > .Machine$integer.max) {
  stop("kappa must be a finite integer between 2 and .Machine$integer.max.")
}

########################################
## Libraries                          ##
########################################

require("ksm")        # SPD kernel density tools (rWAR, bandwidth_optim, etc.)
require("parallel")   # detectCores
require("doFuture")   # parallel backend for foreach
require("foreach")    # %dopar%
require("future")     # parallel plans

########################################
## Inverse Wishart KDE and selectors  ##
########################################

# These helpers implement the manuscript's variable-interchanged kernel:
# f_IW(S) = mean_i K_{1/b + d + 1, b X_i}(S).
# ksm currently supports only Wishart, smlnorm and smnorm kernels.
# To reuse just these functions in another simulation script, evaluate this
# section without running the experiment below.

iw_logsumexp <- function(z) {
  zmax <- max(z)
  if (!is.finite(zmax)) return(zmax)
  zmax + log(sum(exp(z - zmax)))
}

iw_lmgamma <- function(a, d) {
  d * (d - 1) * log(pi) / 4 + sum(lgamma(a - (seq_len(d) - 1) / 2))
}

iw_spd_array <- function(x) {
  if (is.matrix(x)) x <- array(x, dim = c(nrow(x), ncol(x), 1L))
  if (length(dim(x)) != 3L || dim(x)[1] != dim(x)[2] || dim(x)[3] < 1L ||
      any(!is.finite(x))) stop("Expected a finite d x d x n array of SPD matrices.")
  for (i in seq_len(dim(x)[3])) {
    Xi <- x[, , i, drop = FALSE]
    dim(Xi) <- dim(x)[1:2]
    if (!isSymmetric(Xi, tol = 1e-10)) stop("Matrices must be symmetric.")
    chol(Xi) # also validates positive definiteness
  }
  x
}

iw_prepare <- function(x, overlap = FALSE) {
  x <- iw_spd_array(x)
  d <- dim(x)[1]
  n <- dim(x)[3]
  xvec <- matrix(x, nrow = d * d, ncol = n)
  inverses <- matrix(NA_real_, nrow = d * d, ncol = n)
  logdet <- numeric(n)
  for (i in seq_len(n)) {
    Xi <- matrix(x[, , i], nrow = d, ncol = d)
    Ri <- chol(Xi)
    logdet[i] <- 2 * sum(log(diag(Ri)))
    inverses[, i] <- as.vector(chol2inv(Ri))
  }
  out <- list(x = x, d = d, n = n, xvec = xvec,
              inverses = inverses, logdet = logdet)
  if (overlap) {
    logdet_sum <- matrix(NA_real_, nrow = n, ncol = n)
    for (i in seq_len(n)) {
      for (j in seq_len(i)) {
        Zij <- matrix(x[, , i] + x[, , j], nrow = d, ncol = d)
        logdet_sum[i, j] <- logdet_sum[j, i] <- 2 * sum(log(diag(chol(Zij))))
      }
    }
    out$logdet_sum <- logdet_sum
  }
  out
}

# Construct once per sample and bandwidth before numerical integration.
# The returned function accepts evaluation matrices x and an optional log flag.
make_kdens_inverseWishart <- function(xs, b) {
  if (length(b) != 1L || !is.finite(b) || b <= 0) stop("b must be positive.")
  sample <- iw_prepare(xs)
  d <- sample$d
  alpha <- 1 / (2 * b)
  a <- alpha + (d + 1) / 2
  log_weights <- -a * sample$logdet -
    d * a * log(2 * b) - iw_lmgamma(a, d)
  log_n <- log(sample$n)
  
  function(x, log = FALSE) {
    x <- iw_spd_array(x)
    if (dim(x)[1] != d) stop("x and xs must have the same matrix dimension.")
    neval <- dim(x)[3]
    xvec <- matrix(x, nrow = d * d, ncol = neval)
    logdet <- vapply(seq_len(neval), function(i) {
      2 * sum(base::log(diag(chol(matrix(x[, , i], nrow = d, ncol = d)))))
    }, numeric(1))
    # Rows are evaluation points; columns are observations in the KDE.
    log_kernel <- outer(alpha * logdet, log_weights, "+") -
      crossprod(xvec, sample$inverses) / (2 * b)
    log_density <- apply(log_kernel, 1L, iw_logsumexp) - log_n
    if (log) log_density else exp(log_density)
  }
}

# x: evaluation matrices; xs: observations. All densities use Lebesgue measure.
kdens_inverseWishart <- function(x, xs, b, log = FALSE) {
  make_kdens_inverseWishart(xs, b)(x, log = log)
}

bandwidth_optim_inverseWishart <- function(
    x, criterion = c("lscv", "lcv"), tol = 1e-4,
    h = 1L, bounds = c(1e-4, 10)
) {
  criterion <- match.arg(criterion)
  if (length(bounds) != 2L || any(!is.finite(bounds))) stop("Invalid bandwidth bounds.")
  bounds <- sort(bounds)
  if (bounds[1] <= 0 || bounds[2] <= bounds[1]) stop("Invalid bandwidth bounds.")
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0) stop("tol must be positive.")
  if (length(h) != 1L || !is.finite(h) || h < 1 || h != as.integer(h)) {
    stop("h must be a positive integer.")
  }
  cache <- iw_prepare(x, overlap = criterion == "lscv")
  d <- cache$d
  n <- cache$n
  keep <- abs(outer(seq_len(n), seq_len(n), "-")) >= h
  nkeep <- rowSums(keep)
  if (any(nkeep == 0L)) stop("h leaves at least one validation point without observations.")
  trace_mat <- crossprod(cache$xvec, cache$inverses)
  logdet_pair <- outer(cache$logdet, cache$logdet, "+")
  cdim <- (d + 1) / 2
  rdim <- d * cdim
  
  objective <- function(b) {
    alpha <- 1 / (2 * b)
    a <- alpha + cdim
    log_kernel <- outer(alpha * cache$logdet, -a * cache$logdet, "+") -
      trace_mat / (2 * b) - d * a * log(2 * b) - iw_lmgamma(a, d)
    log_kernel[!keep] <- -Inf
    log_cv <- apply(log_kernel, 1L, iw_logsumexp) - log(nkeep)
    if (criterion == "lcv") return(-mean(log_cv))
    
    # This squared-density integral equals the overlap for the Wishart KDE.
    log_overlap <- iw_lmgamma(2 * alpha + cdim, d) -
      2 * iw_lmgamma(a, d) - rdim * log(2 * b) +
      alpha * logdet_pair - (2 * alpha + cdim) * cache$logdet_sum
    log_l2 <- iw_logsumexp(as.vector(log_overlap)) - 2 * log(n)
    # Match ksm's LSCV convention: pool all retained validation pairs.
    log_cross <- log(2) + iw_logsumexp(log_kernel[keep]) - log(sum(nkeep))
    
    # A strictly increasing signed-log transform preserves the LSCV minimizer
    # and avoids overflow in the difference of two positive sums.
    if (log_l2 == log_cross) return(0)
    sign_q <- if (log_l2 > log_cross) 1 else -1
    log_abs_q <- max(log_l2, log_cross) +
      log(-expm1(-abs(log_l2 - log_cross)))
    sign_q * (max(log_abs_q, 0) + log1p(exp(-abs(log_abs_q))))
  }
  
  # Match ksm::bandwidth_optim, including its endpoint-bound expansions.
  opt <- stats::optimize(objective, interval = bounds, tol = tol)$minimum
  if ((opt - bounds[1]) < 1e-3) {
    opt <- stats::optimize(objective, interval = c(1e-8, bounds[1] * 1.1),
                           tol = tol)$minimum
  }
  if ((bounds[2] - opt) < 1e-3) {
    opt <- stats::optimize(objective, interval = c(bounds[2] * 0.9, 1e4),
                           tol = tol)$minimum
  }
  opt
}

########################################
## Simulation helper: WAR(1)          ##
########################################

# Simulation study serial dependence from WAR(1)
# Returns a 2 x 2 x n array of SPD matrices
simu_rWAR <- function(
    n,
    Mmod = c("M1", "M2", "M3"),
    Smod = c("S1", "S2", "S3"),
    K = kappa
) {
  Mmod <- match.arg(Mmod)
  Smod <- match.arg(Smod)
  n    <- as.integer(n)
  K    <- as.integer(K)
  stopifnot(n >= 1, K >= 2) # otherwise matrix not invertible / degenerate
  
  # Autoregressive coefficient matrix M
  if (Mmod == "M1") {
    M <- matrix(
      c(0.9, 0, 1, 0),
      byrow = TRUE, nrow = 2, ncol = 2
    )
  } else if (Mmod == "M2") {
    M <- matrix(
      c(0.3, -0.3, -0.3, 0.3),
      byrow = TRUE, nrow = 2, ncol = 2
    )
  } else { # "M3"
    M <- diag(rep(0.5, 2))
  }
  
  # Innovation covariance matrix S (as specified in the manuscript)
  if (Smod == "S1") {
    S <- matrix(0.9, 2, 2) + diag(rep(0.1, 2))
  } else if (Smod == "S2") {
    S <- matrix(0.95, 2, 2) + diag(rep(0.05, 2))
  } else { # "S3"
    S <- matrix(
      c(1, 0.99, 0.99, 1),
      byrow = TRUE, nrow = 2, ncol = 2
    )
  }
  
  ksm::rWAR(
    n = n,
    M = M,
    S = S,
    K = K
  )
}

########################################
## Hyper-parameters                   ##
########################################

# Replications
RR   <- 1:100

# Sample sizes
nobs <- c(100L, 200L, 300L)

# WAR(1) models: 3 × 3 = 9 models
Mmod <- c("M1", "M2", "M3")
Smod <- c("S1", "S2", "S3")

# Kernels and criteria / methods
# Wishart        + {lscv, lcv} -> W_lscv, W_lcv
# inverseWishart + {lscv, lcv} -> IW_lscv, IW_lcv
# smlnorm + {lscv, lcv}  -> LG_lscv, LG_lcv
# smnorm  + {lscv, lcv}  -> G_lscv, G_lcv
kernels  <- c("Wishart", "inverseWishart", "smlnorm", "smnorm")
criteria <- c("lscv", "lcv")

# Output path
path <- getwd()

########################################
## Parallel setup on laptop           ##
########################################

# Use all available cores minus 1 (at least 1 worker)
n_cores_total <- parallel::detectCores()
n_workers     <- max(1L, n_cores_total - 1L)

doFuture::registerDoFuture()
future::plan(future::multisession, workers = n_workers)

cat("Total cores:", n_cores_total, " -> using", n_workers, "workers\n")

########################################
## Storage for raw timing results     ##
########################################

raw_results <- data.frame(
  replication = integer(),
  nobs        = integer(),
  Mmod        = character(),
  Smod        = character(),
  kernel      = character(),
  criterion   = character(),
  method      = character(),
  bandwidth   = numeric(),
  time_sec    = numeric(),
  stringsAsFactors = FALSE
)

########################################
## Main parallel loop over RR         ##
########################################

start_time <- Sys.time()

res_list <- foreach::foreach(
  r = RR,
  .combine  = "rbind",
  .packages = c("ksm")
) %dopar% {
  
  set.seed(r)
  
  out_rows_all <- vector(
    "list",
    length(nobs) * length(Mmod) * length(Smod) * length(kernels) * length(criteria)
  )
  idx_out <- 0L
  
  for (i in seq_along(nobs)) {
    n_i <- nobs[i]
    h_i <- ceiling(n_i^0.25)  # h-lag adjustment for serial dependence
    
    for (j in seq_along(Mmod)) {
      for (k in seq_along(Smod)) {
        
        # Simulate WAR(1) sample
        xs <- simu_rWAR(
          n    = n_i,
          Mmod = Mmod[j],
          Smod = Smod[k],
          K    = kappa
        )
        
        for (l in seq_along(kernels)) {
          ker <- kernels[l]
          
          for (c_idx in seq_along(criteria)) {
            crit <- criteria[c_idx]
            
            # Time ONLY the bandwidth selection
            t0 <- Sys.time()
            band_obj <- try(
              if (ker == "inverseWishart") {
                bandwidth_optim_inverseWishart(
                  x = xs, criterion = crit, h = h_i
                )
              } else {
                ksm:::bandwidth_optim(
                  x         = xs,
                  criterion = crit,
                  h         = h_i,
                  kernel    = ker
                )
              },
              silent = TRUE
            )
            t1 <- Sys.time()
            
            if (inherits(band_obj, "try-error")) {
              band_val <- NA_real_
              time_val <- NA_real_
            } else {
              band_val <- as.numeric(band_obj)
              time_val <- as.numeric(difftime(t1, t0, units = "secs"))
            }
            
            method_name <- paste(
              switch(ker, Wishart = "W", inverseWishart = "IW", smlnorm = "LG", smnorm = "G"),
              crit,
              sep = "_"
            )
            
            idx_out <- idx_out + 1L
            out_rows_all[[idx_out]] <- data.frame(
              replication = r,
              nobs        = n_i,
              Mmod        = Mmod[j],
              Smod        = Smod[k],
              kernel      = ker,
              criterion   = crit,
              method      = method_name,
              bandwidth   = band_val,
              time_sec    = time_val,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  do.call(rbind, out_rows_all)
}

raw_results <- res_list
if (any(!is.finite(raw_results$bandwidth) | !is.finite(raw_results$time_sec))) {
  warning("Some bandwidth selections failed; their timings are missing. Inspect raw_results before reporting the table.")
}

# Reset plan to sequential
future::plan(sequential)

elapsed_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("Total elapsed time:", round(elapsed_time_minutes, 2), "minutes\n")

# ########################################
# ## Save raw timings                   ##
# ########################################
# 
# raw_output_file <- file.path(path, "bandwidth_timing_raw_local.csv")
# write.csv(raw_results, raw_output_file, row.names = FALSE)
# cat("Raw timing results saved to:", raw_output_file, "\n")

########################################
## Summarise computation times        ##
########################################

# Ensure character columns
raw_results$Mmod      <- as.character(raw_results$Mmod)
raw_results$Smod      <- as.character(raw_results$Smod)
raw_results$kernel    <- as.character(raw_results$kernel)
raw_results$criterion <- as.character(raw_results$criterion)
raw_results$method    <- as.character(raw_results$method)

# Mean time by (nobs, model, method)
summary_times <- aggregate(
  time_sec ~ nobs + Mmod + Smod + method,
  data = raw_results,
  FUN  = mean,
  na.rm = TRUE
)

# Also compute sd, median, IQR if needed
summary_sd <- aggregate(
  time_sec ~ nobs + Mmod + Smod + method,
  data = raw_results,
  FUN  = sd,
  na.rm = TRUE
)

summary_median <- aggregate(
  time_sec ~ nobs + Mmod + Smod + method,
  data = raw_results,
  FUN  = median,
  na.rm = TRUE
)

summary_iqr <- aggregate(
  time_sec ~ nobs + Mmod + Smod + method,
  data = raw_results,
  FUN  = IQR,
  na.rm = TRUE
)

# Merge into one table
summary_results <- merge(
  summary_times,
  summary_sd,
  by = c("nobs", "Mmod", "Smod", "method"),
  suffixes = c("_mean", "_sd")
)

summary_results <- merge(
  summary_results,
  summary_median,
  by = c("nobs", "Mmod", "Smod", "method")
)
names(summary_results)[names(summary_results) == "time_sec"] <- "time_sec_median"

summary_results <- merge(
  summary_results,
  summary_iqr,
  by = c("nobs", "Mmod", "Smod", "method")
)
names(summary_results)[names(summary_results) == "time_sec"] <- "time_sec_IQR"

# # Save summary
# summary_output_file <- file.path(path, "bandwidth_timing_summary_local.csv")
# write.csv(summary_results, summary_output_file, row.names = FALSE)
# cat("Summary timing results saved to:", summary_output_file, "\n")

########################################
## Generate LaTeX table (.tex)        ##
########################################

# We use the mean times (in seconds) from `summary_results$time_sec_mean`
sample_sizes  <- c(100L, 200L, 300L)
methods_order <- c("W_lscv", "IW_lscv", "LG_lscv", "G_lscv", "W_lcv", "IW_lcv", "LG_lcv", "G_lcv")
Mmod_vec      <- c("M1", "M2", "M3")
Smod_vec      <- c("S1", "S2", "S3")

# Format mean time in milliseconds
fmt_ms <- function(x) {
  if (is.na(x)) return("--")
  as.character(round(1000 * x))  # seconds -> ms, rounded to integer
}

latex_lines <- c(
  "\\begin{table}[H]",
  paste0(
    "\\caption{Mean bandwidth computation time (in milliseconds) for each ",
    "method, model and sample size, with $\\kappa = ", kappa,
    "$. Bold indicates the smallest mean among the eight competitors for ",
    "the corresponding model and sample size.}"
  ),
  "\\label{tab:mean.time}",
  "\\footnotesize",
  "\\setlength{\\aboverulesep}{1.2pt}",
  "\\setlength{\\belowrulesep}{1.2pt}",
  "\\renewcommand{\\arraystretch}{0.4}",
  "\\setlength{\\tabcolsep}{1.5pt}",
  "\\centering",
  "\\resizebox{\\textwidth}{!}{%",
  "\\begin{tabular}{c|cccccccc|cccccccc|cccccccc}",
  "\\toprule",
  " & \\multicolumn{8}{c|}{$n = 100$} & \\multicolumn{8}{c|}{$n = 200$} & \\multicolumn{8}{c}{$n = 300$} \\\\",
  "\\cmidrule(lr){2-9}\\cmidrule(lr){10-17}\\cmidrule(l){18-25}",
  " & \\multicolumn{4}{c}{lscv} & \\multicolumn{4}{c|}{lcv} & \\multicolumn{4}{c}{lscv} & \\multicolumn{4}{c|}{lcv} & \\multicolumn{4}{c}{lscv} & \\multicolumn{4}{c}{lcv} \\\\",
  "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-13}\\cmidrule(lr){14-17}\\cmidrule(lr){18-21}\\cmidrule(l){22-25}",
  "Model & $W$ & $\\mathrm{IW}$ & $\\mathrm{LG}$ & $\\mathrm{G}$ & $W$ & $\\mathrm{IW}$ & $\\mathrm{LG}$ & $\\mathrm{G}$ & $W$ & $\\mathrm{IW}$ & $\\mathrm{LG}$ & $\\mathrm{G}$ & $W$ & $\\mathrm{IW}$ & $\\mathrm{LG}$ & $\\mathrm{G}$ & $W$ & $\\mathrm{IW}$ & $\\mathrm{LG}$ & $\\mathrm{G}$ & $W$ & $\\mathrm{IW}$ & $\\mathrm{LG}$ & $\\mathrm{G}$ \\\\",
  "\\midrule"
)

# 9 model rows: (M1,S1), (M1,S2), ..., (M3,S3)
for (Mj in Mmod_vec) {
  for (Sj in Smod_vec) {
    # Label like (M_{1}, \Sigma_{1})
    M_idx <- substring(Mj, 2)  # "M1" -> "1"
    S_idx <- substring(Sj, 2)  # "S2" -> "2"
    model_label <- sprintf("$(M_{%s}, \\Sigma_{%s})$", M_idx, S_idx)
    
    row_cells <- c(model_label)
    
    # For each sample size, get 8 means and bold the smallest
    for (n_i in sample_sizes) {
      means_sec <- rep(NA_real_, length(methods_order))
      
      # Collect means for the 8 methods at this (n_i, Mj, Sj)
      for (m_idx in seq_along(methods_order)) {
        meth <- methods_order[m_idx]
        idx <- which(
          summary_results$nobs   == n_i &
            summary_results$Mmod   == Mj  &
            summary_results$Smod   == Sj  &
            summary_results$method == meth
        )
        if (length(idx) > 0L) {
          means_sec[m_idx] <- summary_results$time_sec_mean[idx[1]]
        }
      }
      
      # Rank only complete eight-method comparisons
      if (any(!is.finite(means_sec))) {
        min_val <- NA_real_
      } else {
        min_val <- min(means_sec, na.rm = TRUE)
      }
      
      # Format and bold the smallest mean for this sample size
      for (m_idx in seq_along(methods_order)) {
        t_sec <- means_sec[m_idx]
        cell  <- fmt_ms(t_sec)
        if (!is.na(min_val) && !is.na(t_sec) && t_sec == min_val) {
          cell <- paste0("\\textbf{", cell, "}")
        }
        row_cells <- c(row_cells, cell)
      }
    }
    
    latex_lines <- c(
      latex_lines,
      paste(paste(row_cells, collapse = " & "), "\\\\")
    )
  }
}

latex_lines <- c(
  latex_lines,
  "\\bottomrule",
  "\\end{tabular}}",
  "\\end{table}"
)

tex_file <- file.path(path, "bandwidth_timing_table.tex")
writeLines(latex_lines, con = tex_file)
cat("LaTeX table written to:", tex_file, "\n")
