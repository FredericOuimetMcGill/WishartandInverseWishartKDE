setwd(this.path::here())
library(ggplot2)
library(patchwork)
library(ksm)

## inverse Wishart KDE: f_IW(S) = mean_i K_{1/b + d + 1, b X_i}(S)
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

# For this variable-interchanged kernel, the integrated squared density
# and the pooled h-step cross-validation sum equal their Wishart counterparts.
# Thus the LSCV objective and its minimizing bandwidth are identical.
bandwidth_optim_inverseWishart <- function(
    x, h = 1L, tol = 1e-4, bounds = c(1e-4, 10)
) {
  ksm::bandwidth_optim(
    x = x, criterion = "lscv", kernel = "Wishart",
    h = h, tol = tol, bounds = bounds
  )
}

printPlot <- TRUE
if (printPlot) {
  dir.create("fig", showWarnings = FALSE)
}
ggplot2::set_theme(theme_classic())
update_theme(
  panel.background = element_rect(fill = "white", color = NA),
  panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
  panel.grid.minor = element_line(
    color = "gray90",
    linetype = "dotted",
    linewidth = 0.2
  ),
  legend.position = "bottom",
  legend.text = element_text(size = 11),
  legend.key.size = unit(2, "lines"),
  # legend.background = element_rect(fill = "white", color = NA),
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 12, color = "black"),
  axis.ticks = element_blank()
)

## DATA ANALYSIS
# Load data on realized variance
data(realvar, package = "ksm")
rvar <- data.frame(
  date = as.Date(dimnames(realvar)$date),
  rvarAMZN = realvar[1, 1, ],
  rvarSPY = realvar[2, 2, ],
  rcov = realvar[1, 2, ]
)


# Exploratory data analysis plots
g_vcov <- ggplot(
  rvar |>
    tidyr::pivot_longer(
      cols = c("rvarAMZN", "rcov", "rvarSPY"),
      names_to = "component",
      values_to = "value"
    ) |>
    dplyr::filter(component != "rcov") |>
    dplyr::mutate(
      component = relevel(
        factor(
          component,
          levels = c("rvarAMZN", "rvarSPY"),
          labels = c("var(AMZN)", "var(SPY)")
        ),
        ref = "var(AMZN)"
      )
    ),
  aes(x = date, y = value, color = component)
) +
  geom_line() +
  scale_color_viridis_d(option = "D") +
  scale_y_continuous(
    transform = "log10",
    labels = scales::label_log()
  ) +
  labs(
    x = "date",
    y = "",
    color = NULL
  ) #+
#theme(legend.position = "inside",
#      legend.position.inside = c(0.1,0.9),legend.box.background = NULL)

rvar$rcor <- apply(realvar, 3, function(x) {
  x[1, 2] / sqrt(x[1, 1] * x[2, 2])
})
rvar$determinant <- apply(realvar, 3, function(x) {
  x[1, 1] * x[2, 2] - x[1, 2]^2
})
eigenvalues <- t(apply(realvar, 3, function(x) {
  eigen(x)$values
}))
rvar$eigen_1 <- eigenvalues[, 1]
rvar$eigen_2 <- eigenvalues[, 2]

g_det <- ggplot(
  data = rvar,
  mapping = aes(x = date, y = log(determinant))
) +
  geom_line() +
  scale_y_continuous(
    labels = scales::label_number(
      drop0trailing = TRUE,
      style_negative = "minus"
    )
  )

g_cor <- ggplot(
  data = rvar,
  mapping = aes(x = date, y = rcor)
) +
  geom_line() +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::label_number(
      style_negative = "minus",
      drop0trailing = TRUE
    )
  ) +
  labs(y = "realized correlation")
g_eigen <- ggplot(
  data = rvar |>
    tidyr::pivot_longer(
      cols = tidyr::starts_with("eigen"),
      names_to = "eigenvalues",
      values_to = "values",
      names_prefix = "eigen_"
    ),
  mapping = aes(x = date, y = values, col = eigenvalues)
) +
  geom_line() +
  scale_color_viridis_d() +
  scale_y_continuous(
    transform = "log10",
    labels = scales::trans_format('log10', scales::math_format(10^.x))
  ) +
  labs(y = "eigenvalues (log scale)", color = "") +
  theme(legend.position = "none")

g_eigen + g_det
if (printPlot) {
  ggsave(
    filename = "fig/rvar_eigen_det.pdf",
    width = 10,
    height = 5
  )
}
# Collect the variance legend above the left-hand realized-variance plot.
g_vcov_cor <- (g_vcov + g_cor + plot_layout(guides = "collect")) &
  scale_x_date(labels = function(x) paste0(
    month.abb[as.integer(format(x, "%m"))], ". ", format(x, "%Y"))) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "left"
  )
print(g_vcov_cor)
if (printPlot) {
  ggsave(
    filename = "fig/rvar_vcov_cor.pdf",
    plot = g_vcov_cor,
    width = 10,
    height = 5
  )
}
# Find the optimal bandwidth with h-step cross validation
h <- ceiling(nrow(rvar)^0.25)

bopt_Wishart <- ksm::bandwidth_optim(
  x = realvar,
  criterion = "lscv",
  kernel = "Wishart",
  h = h
)
bopt_logGauss <- ksm::bandwidth_optim(
  x = realvar,
  criterion = "lscv",
  kernel = "smlnorm",
  h = h
)

# inverse Wishart fit with the same h-step LSCV convention.
bopt_inverseWishart <- bandwidth_optim_inverseWishart(x = realvar, h = h)
kdens_IW_fit <- make_kdens_inverseWishart(xs = realvar, b = bopt_inverseWishart)
bandwidth_comparison <- data.frame(
  estimator = c("Wishart", "inverse Wishart"),
  criterion = "lscv",
  h = h,
  bandwidth = c(bopt_Wishart, bopt_inverseWishart)
)
print(bandwidth_comparison)

# Plot sequence
b_seq <- seq(0.005, to = 0.05, length.out = 101)
lscv_crit <- sapply(b_seq, function(b) {
  lscv_kern_Wishart(realvar, b, h = h)
})

# Current ksm returns the scalar LSCV value; older versions return two
# log terms. Retain the representation expected by the original plot below.
if (is.null(dim(lscv_crit))) {
  lscv_crit <- rbind(log(pmax(lscv_crit, 0)), log(pmax(-lscv_crit, 0)))
}

# The plotted ordinate is log(-LSCV(b)), defined where LSCV(b) < 0.
# Maximizing this curve is equivalent to minimizing the LSCV criterion.
g_band <- ggplot(
  data = data.frame(
    x = b_seq,
    y = log(-exp(lscv_crit[1, ]) + exp(lscv_crit[2, ]))
  ),
  mapping = aes(x = x, y = y)
) +
  geom_line() +
  geom_vline(xintercept = bopt_Wishart, linetype = "dashed") +
  labs(
    y = expression(log(-LSCV(b))),
    x = "bandwidth (lag h=4)",
    title = "Wishart"
  )


# W and IW have the same pooled h-step LSCV curve.
lscv_crit_IW <- lscv_crit
g_band_IW <- ggplot(
  data = data.frame(
    x = b_seq,
    y = log(-exp(lscv_crit_IW[1, ]) + exp(lscv_crit_IW[2, ]))
  ),
  mapping = aes(x = x, y = y)
) +
  geom_line() +
  geom_vline(xintercept = bopt_inverseWishart, linetype = "dashed") +
  labs(
    y = expression(log(-LSCV(b))),
    x = paste0("bandwidth (lag h=", h, ")"),
    title = "inverse Wishart"
  )


## Application
ran1 <- range(realvar[1, 1, ])
ran2 <- range(realvar[2, 2, ])
pseq <- seq(0, 1, by = 0.01)
q1 <- seq(from = ran1[1], to = ran1[2], length.out = 100L)
q2 <- seq(from = ran2[1], to = ran2[2], length.out = 100L)
q1 <- quantile(realvar[1, 1, ], probs = pseq)
q2 <- quantile(realvar[2, 2, ], probs = pseq)
# q1 <- exp(seq(from = log(ran1[1]), log(ran1[2]), length.out = 100L))
# q2 <- exp(seq(from = log(ran2[1]), log(ran2[2]), length.out = 100L))
# Create an array with all the matrices
qc <- expand.grid(q1, q2)
# Define the set of s12 values we want to loop over
corv <- seq(0, 0.75, by = 0.25)

S <- array(dim = c(2, 2, nrow(qc) * length(corv)))
k = 0
for (j in seq_along(corv)) {
  for (i in seq_len(nrow(qc))) {
    k <- k + 1
    S[1, 1, k] <- qc[i, 1]
    S[2, 2, k] <- qc[i, 2]
    S[1, 2, k] <- S[2, 1, k] <- sqrt(qc[i, 1] * qc[i, 2]) * corv[j]
  }
}
kdens_v <- ksm::kdens_Wishart(
  x = S,
  xs = realvar,
  b = bopt_Wishart,
  log = FALSE
)
res <- array(
  kdens_v,
  dim = c(length(q1), length(q2), length(corv))
)

S_density <- S
log_kdens_v_IW <- kdens_IW_fit(x = S_density, log = TRUE)
kdens_v_IW <- exp(log_kdens_v_IW)
res_IW <- array(kdens_v_IW, dim = dim(res))
log_res_IW <- array(log_kdens_v_IW, dim = dim(res))

if (printPlot) {
  for (k in seq_along(corv)) {
    pdf(
      file = paste0("fig/rvar_dens_W_", k, ".pdf"),
      width = 7,
      height = 5
    )
    par(mar = c(4, 4, 2.8, 0)) # Reduce the margins
    ext <- log10(range(res[,, k]))
    filled.contour(
      x = pseq,
      y = pseq,
      z = log10(res[,, k]) - ext[2],
      levels = c(
        seq(-50, -20, length.out = 5),
        seq(-19, 0, by = 1)
      ),
      color.palette = function(n) {
        viridis::magma(n)
      },
      xlab = "quantile level of realized variance of AMZN",
      ylab = "quantile level of realized variance of SPY",
      # xlab = expression(s[11]),
      # ylab = expression(s[22]),
      cex.lab = 1.25,
      cex.axis = 1.15,
      mgp = c(2.5, 0.3, 0),
      key.title = title(cex.main = 2), # Legend title size
      key.axes = axis(4, cex.axis = 1.5) # Legend tick labels size
    )
    mtext(
      bquote("W, corr =" ~ .(corv[k])),
      side = 3,
      line = 0.5,
      cex = 1.25,
      adj = 0.36
    )
    
    dev.off()
  }
}

# inverse Wishart contours on the same grid and with the same levels.
if (printPlot) {
  for (k in seq_along(corv)) {
    pdf(
      file = paste0("fig/rvar_dens_IW_", k, ".pdf"),
      width = 7,
      height = 5
    )
    par(mar = c(4, 4, 2.8, 0)) # Reduce the margins
    ext <- range(log_res_IW[,, k]) / log(10)
    filled.contour(
      x = pseq,
      y = pseq,
      z = log_res_IW[,, k] / log(10) - ext[2],
      levels = c(
        seq(-50, -20, length.out = 5),
        seq(-19, 0, by = 1)
      ),
      color.palette = function(n) {
        viridis::magma(n)
      },
      xlab = "quantile level of realized variance of AMZN",
      ylab = "quantile level of realized variance of SPY",
      # xlab = expression(s[11]),
      # ylab = expression(s[22]),
      cex.lab = 1.25,
      cex.axis = 1.15,
      mgp = c(2.5, 0.3, 0),
      key.title = title(cex.main = 2), # Legend title size
      key.axes = axis(4, cex.axis = 1.5) # Legend tick labels size
    )
    mtext(
      bquote("IW, corr =" ~ .(corv[k])),
      side = 3,
      line = 0.5,
      cex = 1.25,
      adj = 0.36
    )
    
    dev.off()
  }
}

# Different plot idea - get 2d profile likelihood
corv <- seq(-0.2, 0.99, by = 0.01)
prof_mat <- matrix(NA, nrow = length(q1), ncol = length(q2))
prof_mat_IW <- matrix(NA, nrow = length(q1), ncol = length(q2))
S <- array(dim = c(2, 2, length(corv)))
for (i in seq_along(q1)) {
  for (j in seq_along(q2)) {
    for (k in seq_along(corv)) {
      S[1, 1, k] <- q1[i]
      S[2, 2, k] <- q2[j]
      S[1, 2, k] <- S[2, 1, k] <- sqrt(q1[i] * q2[j]) * corv[k]
    }
    
    kdens_v <- ksm::kdens_Wishart(
      x = S,
      xs = realvar,
      b = bopt_Wishart
    )
    prof_mat[i, j] <- corv[which.max(kdens_v)]
    
    log_kdens_profile_IW <- kdens_IW_fit(x = S, log = TRUE)
    prof_mat_IW[i, j] <- corv[which.max(log_kdens_profile_IW)]
  }
}

prof <- data.frame(cbind(expand.grid(x = pseq, y = pseq), c(prof_mat)))
colnames(prof) <- c("x", "y", "correlation")
# Use the same colour scale for the W and IW correlation legends.
correlation_limits <- range(prof_mat, prof_mat_IW, finite = TRUE)
g_prof <- ggplot(data = prof, mapping = aes(x = x, y = y, fill = correlation)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(
    option = "A",
    limits = correlation_limits,
    labels = scales::label_number(drop0trailing = TRUE),
    guide = guide_colourbar(
      direction = "vertical",
      theme = theme(
        legend.key.width = grid::unit(3, "mm"),
        legend.key.height = grid::unit(45, "mm"),
        legend.title.position = "right",
        legend.title = element_text(angle = 90, margin = margin(l = 3)),
        legend.text.position = "right",
        legend.margin = margin(0, 0, 0, 0)
      )
    )
  ) +
  scale_x_continuous(
    expand = expansion(),
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  scale_y_continuous(
    expand = expansion(),
    labels = scales::label_number(drop0trailing = TRUE)
  ) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    expand = FALSE
  ) +
  labs(
    x = "quantile level of realized variance of AMZN",
    y = "quantile level of realized variance of SPY",
    title = "Wishart"
  )
# Match the original correlation-profile plot and combined PDF.
prof_IW <- data.frame(
  expand.grid(x = pseq, y = pseq),
  correlation = c(prof_mat_IW)
)
g_prof_IW <- (g_prof + prof_IW) +
  labs(title = "inverse Wishart") +
  theme(legend.position = "none")

# Build both plot rows without legends, with identical panel and title spacing.
profile_plot_width <- 10
profile_plot_height <- 5
g_band_profile_W_body <- patchwork::patchworkGrob(
  g_band + (g_prof + theme(legend.position = "none"))
)
g_band_profile_IW_body <- patchwork::patchworkGrob(g_band_IW + g_prof_IW)
shared_widths <- grid::unit.pmax(
  g_band_profile_W_body$widths, g_band_profile_IW_body$widths
)
shared_heights <- grid::unit.pmax(
  g_band_profile_W_body$heights, g_band_profile_IW_body$heights
)
g_band_profile_W_body$widths <- g_band_profile_IW_body$widths <- shared_widths
g_band_profile_W_body$heights <- g_band_profile_IW_body$heights <- shared_heights

# Add a narrow legend column beside each profile panel. The PDF width grows
# by exactly this column's width, preserving every original panel dimension.
legend_plot <- ggplotGrob(
  g_prof + theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.justification = "center",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0)
  )
)
correlation_legend <- legend_plot$grobs[[
  which(legend_plot$layout$name == "guide-box-right")
]]
legend_gap <- grid::unit(1.5, "mm")
profile_legend_width <- grid::convertWidth(
  grid::grobWidth(correlation_legend) + legend_gap, "inches", valueOnly = TRUE
)
# Fix the added column in physical units to match the extra PDF width exactly.
legend_column_width <- grid::unit(profile_legend_width, "inches")
profile_pdf_width <- profile_plot_width + profile_legend_width

add_profile_legend <- function(plot_body) {
  right_panel <- plot_body$layout[plot_body$layout$name == "panel-2", ]
  stopifnot(nrow(right_panel) == 1L)
  legend_aligned <- grid::grobTree(
    correlation_legend,
    vp = grid::viewport(
      x = legend_gap, y = 0.5,
      width = grid::grobWidth(correlation_legend),
      height = grid::grobHeight(correlation_legend),
      just = c("left", "center")
    )
  )
  plot_with_legend <- gtable::gtable_add_cols(
    plot_body, widths = legend_column_width, pos = right_panel$r
  )
  gtable::gtable_add_grob(
    plot_with_legend, legend_aligned,
    t = right_panel$t, b = right_panel$b, l = right_panel$r + 1L,
    clip = "off", name = "correlation-legend"
  )
}
g_band_profile_W <- add_profile_legend(g_band_profile_W_body)
g_band_profile_IW <- add_profile_legend(g_band_profile_IW_body)

grid::grid.newpage()
grid::grid.draw(g_band_profile_W)
grid::grid.newpage()
grid::grid.draw(g_band_profile_IW)
if (printPlot) {
  ggsave(
    filename = "fig/rvar_bandwidth_W.pdf",
    plot = g_band_profile_W,
    width = profile_pdf_width,
    height = profile_plot_height
  )
  ggsave(
    filename = "fig/rvar_bandwidth_IW.pdf",
    plot = g_band_profile_IW,
    width = profile_pdf_width,
    height = profile_plot_height
  )
  # Preserve the previous plot scale when including the wider PDFs in LaTeX.
  cat("To keep the same panel sizes in LaTeX, multiply the previous includegraphics width by ",
      format(profile_pdf_width / profile_plot_width, digits = 8), ".\n", sep = "")
}

# pdf(
#   file = paste0("fig/rvar_profile_W.pdf"),
#   width = 8,
#   height = 6
# )
# par(mar = c(4, 4, 0.5, 0))
# filled.contour(
#   x = pseq,
#   y = pseq,
#   z = prof_mat,
#   color.palette = function(n) {
#     viridis::magma(n)
#   },
#   levels = seq(0, 1, length.out = 21L),
#   xlab = "quantile level of realized variance of AMZN",
#   ylab = "quantile level of realized variance of SPY",
#   cex.lab = 1.25,
#   cex.axis = 1.15,
#   mgp = c(2.5, 0.3, 0),
#   key.title = title(cex.main = 2), # Legend title size
#   key.axes = axis(4, cex.axis = 1.5) # Legend tick labels size
# )
# dev.off()

kdens_v <- ksm::kdens_Wishart(
  x = S,
  xs = realvar,
  b = bopt_Wishart,
  log = FALSE
)

# Evaluate IW on the same final profile grid as the preceding Wishart call.
kdens_v_IW <- kdens_IW_fit(x = S, log = FALSE)

# Restore the full grid for the original 10,000-point timing benchmark.
S <- S_density

# Check timing of different kernel estimators
# Very small differences in d=2

timing <- bench::mark(
  Wishart = ksm::kdens_Wishart(
    x = S[,, 1:10000],
    xs = realvar,
    b = bopt_Wishart,
    log = FALSE
  ),
  "inverse Wishart" = kdens_inverseWishart(
    x = S[,, 1:10000],
    xs = realvar,
    b = bopt_inverseWishart,
    log = FALSE
  ),
  logGauss = ksm::kdens_smlnorm(
    x = S[,, 1:10000],
    xs = realvar,
    b = bopt_logGauss,
    log = FALSE
  ),
  check = FALSE,
  filter_gc = FALSE,
  iterations = 100,
  time_unit = "ms"
)
mean_timing <- unlist(lapply(timing$time, mean))
se_timing <- unlist(lapply(timing$time, sd)) / 10


# Timing summaries in milliseconds for all three benchmarked estimators.
timing_comparison <- data.frame(
  estimator = as.character(timing$expression),
  mean_ms = 1000 * as.numeric(mean_timing),
  se_ms = 1000 * as.numeric(se_timing)
)
print(timing_comparison)


