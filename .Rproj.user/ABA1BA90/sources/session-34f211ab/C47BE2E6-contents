###############################################################################
# KM2IPD single-trial data generator
# - Three-layer modular design:
#   Layer 1: Pure generation (no plotting, no I/O)
#   Layer 2: Derived outputs (axis, risk table, summaries)
#   Layer 3: I/O + plotting wrapper (optional)
###############################################################################

## ---- 0) Libraries & setup ----
suppressPackageStartupMessages({
  library(survival)
  library(survRM2)
  library(IPDfromKM)
  library(SurvdigitizeR)
  library(data.table)
  library(ggplot2)
  library(here)
  library(KMtoIPD)
  library(kmdata)
})



`%||%` <- function(x, y) if (!is.null(x)) x else y

#-----------------------------#
# Scenario library
#-----------------------------#
make_scenario_library <- function(
    mean_time = 10,
    weibull_shapes = c(0.6, 1.0, 2.0),
    lognormal_sigmas = c(1, 2),
    hr_set = c(0.67, 0.85, 1.00)
){
  weibull_scale_from_mean <- function(mean, shape){
    mean / gamma(1 + 1/shape)
  }

  lib <- list()

  # Weibull scenarios (PH via HR mapping)
  for (sh in weibull_shapes) {
    sc <- weibull_scale_from_mean(mean_time, sh)
    shape_tag <- gsub("\\.", "p", sprintf("%.1f", sh))

    for (hr in hr_set) {
      hr_tag <- sprintf("HR%03d", round(hr * 100))
      id <- paste0("W_shape", shape_tag, "_", hr_tag)

      scale_ctrl <- sc
      scale_trt  <- sc / (hr^(1/sh))  # PH: H_trt(t)=HR*H_ctrl(t)

      lib[[id]] <- list(
        dist = "weibull",
        label = sprintf("Weibull(shape=%.1f, mean=%.1f) with PH HR=%.2f", sh, mean_time, hr),
        control = list(shape = sh, scale = scale_ctrl),
        treat   = list(shape = sh, scale = scale_trt),
        effect  = list(type = "PH", HR = hr),
        mean_time = mean_time
      )
    }
  }

  # Lognormal scenarios (MR shift; mean fixed)
  for (sg in lognormal_sigmas) {
    mu_ctrl <- log(mean_time) - 0.5 * sg^2

    for (hr in hr_set) {
      MR <- 1 / hr
      mr_tag <- sprintf("MR%03d", round(MR * 100))
      sg_tag <- sprintf("S%d", as.integer(sg))
      id <- paste0("LN_", sg_tag, "_", mr_tag)

      mu_trt <- mu_ctrl + log(MR)

      lib[[id]] <- list(
        dist = "lognormal",
        label = sprintf("Lognormal(sigma=%.0f, mean=%.1f) with median ratio=%.2f", sg, mean_time, MR),
        control = list(meanlog = mu_ctrl, sdlog = sg),
        treat   = list(meanlog = mu_trt,  sdlog = sg),
        effect  = list(type = "MedianRatio", MR = MR),
        mean_time = mean_time
      )
    }
  }

  lib
}

#==============================================================#
# Layer 1: Pure generation
#==============================================================#

# Build arm-specific event-time RNG functions from a scenario definition
make_event_generators <- function(scn){
  stopifnot(is.list(scn), !is.null(scn$dist))
  dist <- scn$dist

  rT_control <- if (dist == "weibull") {
    function(n) stats::rweibull(n, shape = scn$control$shape, scale = scn$control$scale)
  } else if (dist == "lognormal") {
    function(n) stats::rlnorm(n, meanlog = scn$control$meanlog, sdlog = scn$control$sdlog)
  } else {
    stop("Unsupported dist in scenario: ", dist)
  }

  rT_treat <- if (dist == "weibull") {
    function(n) stats::rweibull(n, shape = scn$treat$shape, scale = scn$treat$scale)
  } else {
    function(n) stats::rlnorm(n, meanlog = scn$treat$meanlog, sdlog = scn$treat$sdlog)
  }

  list(dist = dist, rT_control = rT_control, rT_treat = rT_treat)
}

# Generate latent event times for both arms
generate_latent_event_times <- function(scn, n_control, n_treatment){
  stopifnot(n_control > 0, n_treatment > 0)
  gen <- make_event_generators(scn)

  T_control <- gen$rT_control(as.integer(n_control))
  T_treat   <- gen$rT_treat(as.integer(n_treatment))

  list(
    dist = gen$dist,
    T_control = T_control,
    T_treat = T_treat
  )
}

# Create deterministic random streams for censoring, so calibration is stable/reproducible
make_censor_uniforms <- function(n_control, n_treatment){
  n_control <- as.integer(n_control)
  n_treatment <- as.integer(n_treatment)

  list(
    U_C_control = stats::runif(n_control),
    U_C_treat   = stats::runif(n_treatment),

    # informative censoring: selection + within-(0, Ti) time
    U_inf_sel_control  = stats::runif(n_control),
    U_inf_sel_treat    = stats::runif(n_treatment),
    U_inf_time_control = stats::runif(n_control),
    U_inf_time_treat   = stats::runif(n_treatment)
  )
}

# Generate censoring times C given tau and fixed random streams
make_censor_times <- function(
    tau,
    censoring = c("random","informative","exp","front","back"),
    target_censoring = 0.25,
    T_control,
    T_treat,
    U
){
  censoring <- match.arg(censoring)
  stopifnot(is.finite(tau), tau > 0)
  stopifnot(target_censoring >= 0, target_censoring <= 0.95)

  # No staggered entry: everyone has the same administrative follow-up cap
  tauA <- rep(tau, length(T_control))
  tauB <- rep(tau, length(T_treat))

  if (censoring == "random") {
    # Uniform(0, tau) LTFU clock
    C_control <- U$U_C_control * tauA
    C_treat   <- U$U_C_treat   * tauB

  } else if (censoring == "exp") {
    # Exponential LTFU clock with P(C <= tau) approximately target_censoring
    p0 <- max(0, min(0.95, target_censoring))
    lambda <- -log(1 - p0) / max(1e-9, tau)
    C_control <- -log(1 - U$U_C_control) / lambda
    C_treat   <- -log(1 - U$U_C_treat)   / lambda

  } else if (censoring == "front") {
    # Early censoring: Beta(a<1, b>1) concentrated near 0
    C_control <- stats::qbeta(U$U_C_control, 0.7, 2.0) * tauA
    C_treat   <- stats::qbeta(U$U_C_treat,   0.7, 2.0) * tauB

  } else if (censoring == "back") {
    # Late censoring: Beta(a>1, b<1) concentrated near 1
    C_control <- stats::qbeta(U$U_C_control, 2.0, 0.7) * tauA
    C_treat   <- stats::qbeta(U$U_C_treat,   2.0, 0.7) * tauB

  } else {
    # Informative censoring:
    # With probability p_inf, censor occurs uniformly BEFORE the (possibly administratively truncated) event time.
    p_inf <- max(0, min(0.95, target_censoring))
    TiA <- pmin(T_control, tauA)
    TiB <- pmin(T_treat,   tauB)

    C_control <- ifelse(U$U_inf_sel_control < p_inf, U$U_inf_time_control * TiA, Inf)
    C_treat   <- ifelse(U$U_inf_sel_treat   < p_inf, U$U_inf_time_treat   * TiB, Inf)
  }

  list(C_control = C_control, C_treat = C_treat, tauA = tauA, tauB = tauB)
}

# Realize observed IPD given latent T, censor times C, and admin tau
realize_ipd <- function(
    tau,
    censoring,
    target_censoring,
    T_control,
    T_treat,
    U
){
  CC <- make_censor_times(
    tau = tau,
    censoring = censoring,
    target_censoring = target_censoring,
    T_control = T_control,
    T_treat = T_treat,
    U = U
  )

  C_control <- CC$C_control
  C_treat   <- CC$C_treat
  tauA <- CC$tauA
  tauB <- CC$tauB

  # Observed time and status
  X_control <- pmin(T_control, C_control, tauA)
  D_control <- as.integer(T_control <= C_control & T_control <= tauA)

  X_treat <- pmin(T_treat, C_treat, tauB)
  D_treat <- as.integer(T_treat <= C_treat & T_treat <= tauB)

  ipd <- data.frame(
    time   = c(X_control, X_treat),
    status = c(D_control, D_treat),
    arm    = c(rep("Control", length(X_control)), rep("Treatment", length(X_treat))),
    stringsAsFactors = FALSE
  )
  ipd <- ipd[order(ipd$arm, ipd$time), , drop = FALSE]

  list(
    ipd = ipd,
    censor_prop = mean(ipd$status == 0),
    n_events = sum(ipd$status == 1),
    C_control = C_control,
    C_treat = C_treat,
    tauA = tauA,
    tauB = tauB
  )
}

# Calibrate tau to achieve target censor proportion (deterministic calibration)
calibrate_tau <- function(
    censoring,
    target_censoring,
    scn,
    T_control,
    T_treat,
    U,
    tol = 1e-4,
    tau_cap = 1e5
){
  censoring <- match.arg(censoring, c("random","informative","exp","front","back"))
  stopifnot(target_censoring >= 0, target_censoring <= 0.95)

  # Informative: treat target_censoring as p_inf, and choose tau long enough so admin censor is negligible.
  if (censoring == "informative") {
    base <- (scn$mean_time %||% 10)
    tau_star <- max(5, 3 * base)

    # Inflate tau until overall censor proportion is close to p_inf from above
    tol_inf <- 0.01
    for (k in 1:12) {
      cp <- realize_ipd(tau_star, censoring, target_censoring, T_control, T_treat, U)$censor_prop
      if (cp - target_censoring <= tol_inf) break
      tau_star <- min(tau_star * 1.5, tau_cap)
    }
    return(tau_star)
  }

  # Otherwise: root solve for tau so censor_prop(tau) - target = 0
  f_root <- function(tau){
    realize_ipd(tau, censoring, target_censoring, T_control, T_treat, U)$censor_prop - target_censoring
  }

  tau_lo <- 1e-6
  tau_hi <- max(5, (scn$mean_time %||% 10))

  f_lo <- f_root(tau_lo)
  f_hi <- f_root(tau_hi)

  iter <- 0L
  while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && iter < 25L) {
    tau_hi <- min(tau_hi * 1.8, tau_cap)
    f_hi <- f_root(tau_hi)
    iter <- iter + 1L
    if (tau_hi >= tau_cap) break
  }

  if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi <= 0) {
    stats::uniroot(f_root, interval = c(tau_lo, tau_hi), tol = tol)$root
  } else {
    warning(sprintf(
      "Could not bracket tau for target_censoring=%.3f (censoring=%s); using tau=%.3f.",
      target_censoring, censoring, tau_hi
    ))
    min(tau_hi, tau_cap)
  }
}

#==============================================================#
# Layer 2: Derived outputs (axis, risk table, summaries)
#==============================================================#

# Make a "pretty" x-axis end for plotting/digitization
compute_pretty_axis <- function(ipd, multiple = 8L){
  stopifnot(is.data.frame(ipd), all(c("time","status","arm") %in% names(ipd)))

  x_end_raw <- max(ipd$time, na.rm = TRUE)

  pretty_x_end <- function(x, multiple = 8L){
    x1 <- ceiling(as.numeric(x)[1])
    as.numeric(ceiling(x1 / multiple) * multiple)
  }

  x_end_plot <- pretty_x_end(x_end_raw, multiple = multiple)
  x_increment <- x_end_plot / multiple
  ticks <- seq(0, x_end_plot, by = x_increment)  # exactly multiple+1 points

  list(
    x_end_raw = x_end_raw,
    x_end_plot = x_end_plot,
    x_increment = x_increment,
    ticks = ticks
  )
}

# Risk table at specific times (aligned to axis ticks)
make_risk_table <- function(ipd, times){
  stopifnot(is.data.frame(ipd), all(c("time","status","arm") %in% names(ipd)))
  stopifnot(is.numeric(times), length(times) >= 2)

  nrisk <- function(tt, tvec) sum(tvec >= tt)

  time_ctrl <- ipd$time[ipd$arm == "Control"]
  time_trt  <- ipd$time[ipd$arm == "Treatment"]

  rt <- data.frame(time = times)
  rt$nrisk_control <- vapply(times, nrisk, numeric(1), tvec = time_ctrl)
  rt$nrisk_treat   <- vapply(times, nrisk, numeric(1), tvec = time_trt)
  rt$nrisk_all     <- rt$nrisk_control + rt$nrisk_treat
  rt
}

# Build a compact totals list for downstream benchmarking/digitization
summarize_trial <- function(scn, scenario_id, ipd, tau_star, axis, risk_table){
  stopifnot(is.data.frame(ipd), all(c("time","status","arm") %in% names(ipd)))

  status_ctrl <- ipd$status[ipd$arm == "Control"]
  status_trt  <- ipd$status[ipd$arm == "Treatment"]

  list(
    scenario_id = scenario_id,
    scenario_label = scn$label,
    dist = scn$dist,
    effect = scn$effect,

    tau_study_end = tau_star,

    x_end_raw  = axis$x_end_raw,
    x_end_plot = axis$x_end_plot,
    x_increment = axis$x_increment,

    censor_prop = mean(ipd$status == 0),
    total_events_ctrl = sum(status_ctrl == 1),
    total_events_trt  = sum(status_trt  == 1),

    risk_table_times = risk_table$time
  )
}

#==============================================================#
# Layer 3: Optional plotting + I/O wrappers
#==============================================================#

# Plot KM curves with a digitization-friendly axis; optionally save as PNG.
plot_km <- function(
    ipd,
    axis,
    add_censor_marks = TRUE,
    file = NULL,
    width = 1400,
    height = 900,
    res = 150
){
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required for plot_km(). Please install it.")
  }

  if (!is.null(file)) {
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  graphics::plot(
    0, 0, type = "n",
    xlim = c(0, axis$x_end_plot),
    ylim = c(0, 1),
    xlab = "Time",
    ylab = "Survival",
    xaxt = "n", yaxt = "n",
    xaxs = "i",  # <-- KEY: no x-axis expansion beyond xlim
    bty = "l"
  )

  # Integer ticks
  x_ticks <- axis$ticks
  y_ticks <- seq(0, 1, by = 0.2)
  graphics::axis(1, at = x_ticks, labels = as.integer(x_ticks))
  graphics::axis(2, at = y_ticks, labels = y_ticks)

  draw_arm <- function(dat, col){
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
    ss  <- summary(fit)
    graphics::lines(
      stats::stepfun(ss$time, c(1, ss$surv)),
      col = col, lwd = 2.5, do.points = FALSE
    )

    if (isTRUE(add_censor_marks)) {
      cen_times <- dat$time[dat$status == 0]
      if (length(cen_times)) {
        xx <- c(0, ss$time); yy <- c(1, ss$surv)
        ux <- !duplicated(xx); xx <- xx[ux]; yy <- yy[ux]
        cen_surv <- stats::approx(
          x = xx, y = yy, xout = cen_times,
          method = "constant", f = 0, rule = 2
        )$y
        graphics::points(cen_times, cen_surv, pch = 3, col = "black", cex = 0.8)
      }
    }
  }

  draw_arm(ipd[ipd$arm == "Control",   c("time","status")], col = "#1f77b4")
  draw_arm(ipd[ipd$arm == "Treatment", c("time","status")], col = "#d62728")

  invisible(TRUE)
}

# Write outputs to disk (CSV + optional KM plot)
write_trial_outputs <- function(
    ipd,
    risk_table,
    totals,
    out_dir,
    prefix,
    save_ipd = FALSE,
    save_risk_table = FALSE,
    save_totals_rds = FALSE
){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (isTRUE(save_ipd)) {
    ipd_path <- file.path(out_dir, paste0(prefix, "_IPD_true.csv"))
    utils::write.csv(ipd, ipd_path, row.names = FALSE)
  }
  if (isTRUE(save_risk_table)) {
    rt_path <- file.path(out_dir, paste0(prefix, "_risk_table.csv"))
    utils::write.csv(risk_table, rt_path, row.names = FALSE)
  }
  if (isTRUE(save_totals_rds)) {
    rds_path <- file.path(out_dir, paste0(prefix, "_totals.rds"))
    saveRDS(totals, rds_path)
  }

  invisible(TRUE)
}

# High-level wrapper: generate + derived outputs + optional plot/I/O
simulate_trial <- function(
    scenario_id,
    scenario_lib = NULL,
    n_control = 250,
    n_treatment = 250,
    censoring = c("random","informative","exp","front","back"),
    target_censoring = 0.30,
    seed = NULL,

    # Optional outputs
    out_dir = NULL,   # if NULL, no disk writing
    prefix = "sim1",
    plot = c("none","km"),
    add_censor_marks = TRUE
){
  censoring <- match.arg(censoring)
  plot <- match.arg(plot)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(scenario_lib)) scenario_lib <- make_scenario_library(mean_time = 10)
  scn <- scenario_lib[[scenario_id]]
  if (is.null(scn)) stop("scenario_id not found in scenario_lib: ", scenario_id)

  # ---- Layer 1: pure generation ----
  latent <- generate_latent_event_times(scn, n_control, n_treatment)
  U <- make_censor_uniforms(n_control, n_treatment)

  tau_star <- calibrate_tau(
    censoring = censoring,
    target_censoring = target_censoring,
    scn = scn,
    T_control = latent$T_control,
    T_treat = latent$T_treat,
    U = U
  )

  R <- realize_ipd(
    tau = tau_star,
    censoring = censoring,
    target_censoring = target_censoring,
    T_control = latent$T_control,
    T_treat = latent$T_treat,
    U = U
  )
  ipd <- R$ipd

  # ---- Layer 2: derived outputs ----
  axis <- compute_pretty_axis(ipd, multiple = 8L)
  risk_table <- make_risk_table(ipd, axis$ticks)
  totals <- summarize_trial(scn, scenario_id, ipd, tau_star, axis, risk_table)

  # ---- Layer 3: optional plot/I/O ----
  km_png <- NULL
  if (plot == "km") {
    if (is.null(out_dir)) {
      # plot to current device
      plot_km(ipd, axis, add_censor_marks = add_censor_marks, file = NULL)
    } else {
      km_png <- file.path(out_dir, paste0(prefix, "_KM_with_censor.png"))
      plot_km(ipd, axis, add_censor_marks = add_censor_marks, file = km_png)
    }
  }

  if (!is.null(out_dir)) {
    write_trial_outputs(
      ipd = ipd,
      risk_table = risk_table,
      totals = totals,
      out_dir = out_dir,
      prefix = prefix
    )
  }

  list(
    ipd = ipd,
    axis = axis,
    risk_table = risk_table,
    totals = totals,
    km_png = km_png
  )
}

###############################################################################
# TEST EXAMPLES (copy-paste runnable)
###############################################################################

# Example 1: Build library and inspect scenario IDs
lib <- make_scenario_library(mean_time = 10)
names(lib)[1:5]  # show a few scenario IDs

# # Example 2: Run a single trial (Weibull, HR=0.85) with random censoring
# res1 <- simulate_trial(
#   scenario_id = "W_shape1p0_HR085",
#   scenario_lib = lib,
#   n_control = 250,
#   n_treatment = 250,
#   censoring = "random",
#   target_censoring = 0.30,
#   seed = 123,
#   out_dir = NULL,
#   plot = "none"
# )
# str(res1$totals)
# res1$risk_table
#
# # Example 3: Same scenario, compare censoring types (front vs back) quickly
# res_front <- simulate_trial(
#   scenario_id = "W_shape1p0_HR085",
#   scenario_lib = lib,
#   n_control = 250,
#   n_treatment = 250,
#   censoring = "front",
#   target_censoring = 0.30,
#   seed = 123,
#   out_dir = NULL,
#   plot = "none"
# )
#
# res_back <- simulate_trial(
#   scenario_id = "W_shape1p0_HR085",
#   scenario_lib = lib,
#   n_control = 250,
#   n_treatment = 250,
#   censoring = "back",
#   target_censoring = 0.30,
#   seed = 123,
#   out_dir = NULL,
#   plot = "none"
# )
#
# c(front = res_front$totals$censor_prop, back = res_back$totals$censor_prop)


# Example 4 (optional): Save outputs + KM plot to a folder
out_dir_demo <- here::here("test")

res3 <- simulate_trial(
  scenario_id = "W_shape0p6_HR085",
  scenario_lib = lib,
  n_control = 250,
  n_treatment = 250,
  censoring = "back",
  target_censoring = 0.30,
  seed = 252766,
  out_dir = out_dir_demo,
  prefix = "W_shape0p6_HR085",
  plot = "km",
  add_censor_marks = FALSE
)
res3$km_png
