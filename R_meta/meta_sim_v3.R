###############################################################################
# meta_sim_v3.R
# KM2IPD meta-analysis simulation (reduced design):
#   Fixed DGM (target HR=0.75) + small grid over (missing% x heterogeneity)
#   Compare 4 analysis routes: AD_only, Hybrid, RM, Pooled_IPD
#
# Core pipeline per replication:
#   1) For k=1..K trials: generate true IPD with trial-specific HR_k
#   2) Estimate per-trial Cox logHR from true IPD ("reported HR")
#   3) Numeric digitization of KM -> digitized points
#   4) IPD reconstruction via IPDfromKM -> reconstructed IPD
#   5) Estimate per-trial Cox logHR from reconstructed IPD
#   6) Impose random missingness on reported HR with probability p_miss
#   7) Meta-analyze under 4 routes + summarize bias/RMSE trends
#
# Notes:
# - Heterogeneity is on logHR scale: logHR_k ~ Normal(log(mu_HR), tau_logHR^2)
# - Baseline event times are Weibull PH with fixed shape + mean time.
# - Censoring is "random" (Uniform(0, tau)), with tau calibrated to hit target
#   censoring proportion (deterministically using fixed uniforms per trial).
###############################################################################

## ---- 0) Package checks ----
.require_pkgs <- function(pkgs){
  ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok)) {
    missing <- pkgs[!ok]
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "), "\n",
      "Please install them, e.g.: install.packages(c(",
      paste(sprintf("\"%s\"", missing), collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.require_pkgs(c("survival", "data.table", "metafor", "ggplot2", "IPDfromKM", "dplyr", "tidyr", "scales"))

suppressPackageStartupMessages({
  library(survival)
  library(data.table)
  library(metafor)
  library(ggplot2)
  library(IPDfromKM)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y

## ---- 1) DGM helpers (Weibull PH + calibrated random censoring) ----

# Convert Weibull mean -> scale
.weibull_scale_from_mean <- function(mean_time, shape){
  mean_time / gamma(1 + 1/shape)
}

# Build a Weibull PH scenario list (control + treatment) from HR
make_scn_weibull_ph <- function(hr, shape = 1.0, mean_time = 10){
  stopifnot(is.finite(hr), hr > 0)
  stopifnot(is.finite(shape), shape > 0)
  stopifnot(is.finite(mean_time), mean_time > 0)

  scale_ctrl <- .weibull_scale_from_mean(mean_time, shape)
  scale_trt  <- scale_ctrl / (hr^(1/shape))  # PH mapping: H_trt(t) = HR * H_ctrl(t)

  list(
    dist = "weibull",
    control = list(shape = shape, scale = scale_ctrl),
    treat   = list(shape = shape, scale = scale_trt),
    effect  = list(type = "PH", HR = hr),
    mean_time = mean_time
  )
}

# Deterministic censoring uniforms (so tau calibration is stable)
make_censor_uniforms_simple <- function(n_control, n_treatment){
  list(
    U_C_control = stats::runif(as.integer(n_control)),
    U_C_treat   = stats::runif(as.integer(n_treatment))
  )
}

# Realize IPD under random censoring:
#   T: latent event times
#   C: Uniform(0, tau) LTFU
#   Admin cap: tau
realize_ipd_random <- function(tau, T_control, T_treat, U){
  stopifnot(is.finite(tau), tau > 0)
  C_control <- U$U_C_control * tau
  C_treat   <- U$U_C_treat   * tau

  X_control <- pmin(T_control, C_control, tau)
  D_control <- as.integer(T_control <= C_control & T_control <= tau)

  X_treat <- pmin(T_treat, C_treat, tau)
  D_treat <- as.integer(T_treat <= C_treat & T_treat <= tau)

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
    total_events_ctrl = sum(D_control == 1),
    total_events_trt  = sum(D_treat  == 1)
  )
}

# Calibrate tau so censor_prop(tau) ~= target_censoring (random censoring)
calibrate_tau_random <- function(
    target_censoring,
    T_control,
    T_treat,
    U,
    tol = 1e-4,
    tau_hi_init = NULL,
    tau_cap = 1e5
){
  stopifnot(target_censoring >= 0, target_censoring <= 0.95)

  f_root <- function(tau){
    realize_ipd_random(tau, T_control, T_treat, U)$censor_prop - target_censoring
  }

  tau_lo <- 1e-6
  tau_hi <- tau_hi_init %||% max(5, stats::median(c(T_control, T_treat), na.rm = TRUE))
  tau_hi <- min(max(tau_hi, 1e-3), tau_cap)

  f_lo <- f_root(tau_lo)
  f_hi <- f_root(tau_hi)

  iter <- 0L
  while (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi > 0 && iter < 30L) {
    tau_hi <- min(tau_hi * 1.8, tau_cap)
    f_hi <- f_root(tau_hi)
    iter <- iter + 1L
    if (tau_hi >= tau_cap) break
  }

  if (is.finite(f_lo) && is.finite(f_hi) && f_lo * f_hi <= 0) {
    stats::uniroot(f_root, interval = c(tau_lo, tau_hi), tol = tol)$root
  } else {
    warning(sprintf(
      "Could not bracket tau for target_censoring=%.3f; using tau=%.3f.",
      target_censoring, tau_hi
    ))
    tau_hi
  }
}

# Pretty axis end + ticks (digitization/risk table)
compute_pretty_axis <- function(ipd, multiple = 8L){
  stopifnot(is.data.frame(ipd), all(c("time","status","arm") %in% names(ipd)))

  x_end_raw <- max(ipd$time, na.rm = TRUE)

  pretty_x_end <- function(x, multiple = 8L){
    x1 <- ceiling(as.numeric(x)[1])
    as.numeric(ceiling(x1 / multiple) * multiple)
  }

  x_end_plot <- pretty_x_end(x_end_raw, multiple = multiple)
  x_increment <- x_end_plot / multiple
  ticks <- seq(0, x_end_plot, by = x_increment)

  list(
    x_end_raw = x_end_raw,
    x_end_plot = x_end_plot,
    x_increment = x_increment,
    ticks = ticks
  )
}

# Risk table at specified times
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

# Generate one trial IPD + derived outputs
simulate_trial_one <- function(
    hr,
    n_control,
    n_treatment,
    shape = 1.0,
    mean_time = 10,
    target_censoring = 0.30,
    seed = NULL
){
  if (!is.null(seed)) set.seed(seed)

  scn <- make_scn_weibull_ph(hr = hr, shape = shape, mean_time = mean_time)

  # latent event times
  T_control <- stats::rweibull(as.integer(n_control), shape = scn$control$shape, scale = scn$control$scale)
  T_treat   <- stats::rweibull(as.integer(n_treatment), shape = scn$treat$shape,   scale = scn$treat$scale)

  # fixed uniforms for deterministic tau calibration
  U <- make_censor_uniforms_simple(n_control, n_treatment)

  tau_star <- calibrate_tau_random(
    target_censoring = target_censoring,
    T_control = T_control,
    T_treat = T_treat,
    U = U
  )

  R <- realize_ipd_random(tau_star, T_control, T_treat, U)
  ipd <- R$ipd

  axis <- compute_pretty_axis(ipd, multiple = 8L)
  risk_table <- make_risk_table(ipd, axis$ticks)

  list(
    hr_true = hr,
    ipd = ipd,
    tau_star = tau_star,
    axis = axis,
    risk_table = risk_table,
    total_events = c(Control = R$total_events_ctrl, Treatment = R$total_events_trt),
    censor_prop = R$censor_prop
  )
}

## ---- 2) Estimation + digitization + reconstruction ----

estimate_logHR_cox <- function(ipd){
  stopifnot(is.data.frame(ipd), all(c("time","status","arm") %in% names(ipd)))
  ipd$arm <- factor(ipd$arm, levels = c("Control", "Treatment"))

  out <- tryCatch({
    fit <- survival::coxph(survival::Surv(time, status) ~ arm, data = ipd, ties = "efron")
    cf  <- stats::coef(fit)
    vc  <- stats::vcov(fit)
    logHR <- unname(cf["armTreatment"])
    se    <- sqrt(unname(vc["armTreatment", "armTreatment"]))
    list(ok = TRUE, logHR = logHR, se = se)
  }, error = function(e){
    list(ok = FALSE, logHR = NA_real_, se = NA_real_, err = conditionMessage(e))
  })

  out
}

# Numeric digitization: evaluate KM step function on a grid, add noise, enforce monotonicity.
# digitize_km_numeric <- function(
#     ipd_true,
#     n_grid = 200,
#     sd_S = 0.01,
#     x_end_plot = NULL
# ){
#   stopifnot(is.data.frame(ipd_true), all(c("time","status","arm") %in% names(ipd_true)))
#   ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))
#
#   axis <- compute_pretty_axis(ipd_true, multiple = 8L)
#   x_end <- x_end_plot %||% axis$x_end_plot
#
#   time_grid <- seq(0, x_end, length.out = as.integer(n_grid))
#
#   km_on_grid <- function(arm_label){
#     dat <- ipd_true[ipd_true$arm == arm_label, c("time","status"), drop = FALSE]
#     fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
#
#     ss <- summary(fit, times = time_grid, extend = TRUE)
#     St <- ss$surv
#
#     # Safety: replace any NA with last observation carried forward
#     if (anyNA(St)) {
#       idx_na <- which(is.na(St))
#       for (i in idx_na) {
#         St[i] <- if (i == 1L) 1 else St[i - 1L]
#       }
#     }
#
#     # Add noise (except enforce start at 1)
#     St_noisy <- St + stats::rnorm(length(St), mean = 0, sd = sd_S)
#     St_noisy[1] <- 1
#
#     # bound + enforce non-increasing
#     St_noisy <- pmin(pmax(St_noisy, 0), 1)
#     St_noisy <- cummin(St_noisy)
#     St_noisy[1] <- 1
#
#     data.table(time = time_grid, St = St_noisy)
#   }
#
#   dC <- km_on_grid("Control")
#   dT <- km_on_grid("Treatment")
#
#   dC[, curve := 1L]
#   dT[, curve := 2L]
#
#   digi_dt <- rbindlist(list(dC, dT), use.names = TRUE)
#   setorder(digi_dt, curve, time)
#
#   list(
#     data = digi_dt,
#     axis = axis,
#     time_grid = time_grid
#   )
# }

# ---- 2) Estimation + digitization + reconstruction ----

# Number of points: 100
digitize_km_numeric <- function(
    ipd_true,
    risk_table,           # 风险表（包含trisk时间点）
    n_digitize = 100,     # 数字化点个数（模拟手动打点）
    sd_S = 0.01
){
  stopifnot(is.data.frame(ipd_true), all(c("time","status","arm") %in% names(ipd_true)))
  stopifnot(is.data.frame(risk_table), "time" %in% names(risk_table))

  ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))

  # 风险表时间点（用于重建）
  trisk <- risk_table$time

  # 生成密集时间网格（用于数字化）
  x_end <- max(trisk)
  time_grid <- seq(0, x_end, length.out = n_digitize)

  km_on_grid <- function(arm_label){
    dat <- ipd_true[ipd_true$arm == arm_label, c("time","status"), drop = FALSE]
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)

    # 在密集时间点提取生存概率
    ss <- summary(fit, times = time_grid, extend = TRUE)
    St <- ss$surv

    # 处理NA
    if (anyNA(St)) {
      idx_na <- which(is.na(St))
      for (i in idx_na) {
        St[i] <- if (i == 1L) 1 else St[i - 1L]
      }
    }

    # 添加噪声（模拟读点误差）
    St_noisy <- St + stats::rnorm(length(St), mean = 0, sd = sd_S)
    St_noisy[1] <- 1

    # 边界限制 + 强制非增
    St_noisy <- pmin(pmax(St_noisy, 0), 1)
    St_noisy <- cummin(St_noisy)
    St_noisy[1] <- 1

    data.table(time = time_grid, St = St_noisy)
  }

  dC <- km_on_grid("Control")
  dT <- km_on_grid("Treatment")

  dC[, curve := 1L]
  dT[, curve := 2L]

  digi_dt <- rbindlist(list(dC, dT), use.names = TRUE)
  setorder(digi_dt, curve, time)

  list(
    data = digi_dt,          # 密集数字化点（用于重建）
    risk_table = risk_table, # 风险表（包含trisk，用于preprocess）
    time_grid = time_grid,
    trisk = trisk
  )
}


# ---- Numeric at risk table ----
# digitize_km_numeric <- function(
#     ipd_true,
#     risk_table,
#     sd_S = 0.1
# ){
#   stopifnot(is.data.frame(ipd_true), all(c("time","status","arm") %in% names(ipd_true)))
#   stopifnot(is.data.frame(risk_table), "time" %in% names(risk_table))
#
#   ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))
#
#   times <- risk_table$time   # 使用风险表的时间点（即trisk）
#
#   km_on_grid <- function(arm_label){
#     dat <- ipd_true[ipd_true$arm == arm_label, c("time","status"), drop = FALSE]
#     fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
#
#     # 在指定时间点提取生存概率
#     ss <- summary(fit, times = times, extend = TRUE)
#     St <- ss$surv
#
#     # 处理NA（可能由于extend=TRUE引入，但一般不会）
#     if (anyNA(St)) {
#       idx_na <- which(is.na(St))
#       for (i in idx_na) {
#         St[i] <- if (i == 1L) 1 else St[i - 1L]
#       }
#     }
#
#     # 添加噪声（除第一个点强制为1）
#     St_noisy <- St + stats::rnorm(length(St), mean = 0, sd = sd_S)
#     St_noisy[1] <- 1
#
#     # 边界限制 + 强制非增
#     St_noisy <- pmin(pmax(St_noisy, 0), 1)
#     St_noisy <- cummin(St_noisy)
#     St_noisy[1] <- 1
#
#     data.table(time = times, St = St_noisy)
#   }
#
#   dC <- km_on_grid("Control")
#   dT <- km_on_grid("Treatment")
#
#   dC[, curve := 1L]
#   dT[, curve := 2L]
#
#   digi_dt <- rbindlist(list(dC, dT), use.names = TRUE)
#   setorder(digi_dt, curve, time)
#
#   list(
#     data = digi_dt,
#     risk_table = risk_table,   # 保留风险表以备重建使用
#     times = times
#   )
# }
# ----

# Two-arm reconstruction using IPDfromKM only
reconstruct_ipd_IPDfromKM <- function(
    digi_dt,
    risk_table,
    total_events = NULL,
    curve_map = c(Control = 1, Treatment = 2)
){
  stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))
  stopifnot(all(c("time", "nrisk_control", "nrisk_treat") %in% names(risk_table)))

  get_curve_xy <- function(curve_id){
    d <- digi_dt[digi_dt$curve == curve_id, c("time", "St")]
    d <- d[order(d$time), ]
    d$St <- pmin(pmax(d$St, 0), 1)
    data.frame(x = d$time, y = d$St)
  }

  get_arm_risk <- function(arm){
    rt <- risk_table[order(risk_table$time), ]
    nrisk <- if (arm == "Control") rt$nrisk_control else rt$nrisk_treat
    data.frame(trisk = rt$time, nrisk = nrisk)
  }

  get_tot_events <- function(arm){
    if (is.null(total_events)) return(NULL)
    if (!is.null(names(total_events)) && arm %in% names(total_events)) return(as.integer(total_events[[arm]]))
    NULL
  }

  reconstruct_one_arm <- function(arm){
    curve_id <- unname(curve_map[[arm]])
    xy <- get_curve_xy(curve_id)
    rt <- get_arm_risk(arm)

    # totalpts = baseline n at risk (ideally at time 0)
    idx0 <- which.min(rt$trisk)
    n0 <- as.integer(rt$nrisk[idx0])

    prep <- IPDfromKM::preprocess(
      dat = xy,
      trisk = rt$trisk,
      nrisk = rt$nrisk,
      totalpts = n0,
      maxy = 1
    )

    armID <- if (arm == "Control") 0 else 1
    tot_ev <- get_tot_events(arm)

    ipd_res <- IPDfromKM::getIPD(prep = prep, armID = armID, tot.events = tot_ev)
    ipd_tbl <- ipd_res$IPD

    # Robust column handling
    time_col   <- if ("time" %in% names(ipd_tbl)) "time" else names(ipd_tbl)[1]
    status_col <- if ("status" %in% names(ipd_tbl)) "status" else names(ipd_tbl)[2]

    out <- data.frame(
      time   = ipd_tbl[[time_col]],
      status = ipd_tbl[[status_col]],
      arm    = arm,
      stringsAsFactors = FALSE
    )

    out[order(out$time), , drop = FALSE]
  }

  ipd_c <- reconstruct_one_arm("Control")
  ipd_t <- reconstruct_one_arm("Treatment")

  ipd_all <- rbind(ipd_c, ipd_t)
  ipd_all$arm <- factor(ipd_all$arm, levels = c("Control", "Treatment"))

  list(ipd = ipd_all, ipd_control = ipd_c, ipd_treat = ipd_t)
}

## ---- 3) Meta-analysis helpers ----

meta_fit_rma <- function(logHR_vec, se_vec, re_method = "REML"){
  ok <- is.finite(logHR_vec) & is.finite(se_vec) & se_vec > 0
  yi <- logHR_vec[ok]
  sei <- se_vec[ok]

  if (length(yi) < 2L) {
    return(list(ok = FALSE, logHR = NA_real_, se = NA_real_, tau2 = NA_real_, I2 = NA_real_, k = length(yi)))
  }

  out <- tryCatch({
    fit <- metafor::rma(yi = yi, sei = sei, method = re_method)
    list(
      ok = TRUE,
      logHR = as.numeric(fit$b),
      se = fit$se,
      tau2 = fit$tau2,
      I2 = fit$I2,
      k = fit$k
    )
  }, error = function(e){
    list(ok = FALSE, logHR = NA_real_, se = NA_real_, tau2 = NA_real_, I2 = NA_real_, k = length(yi), err = conditionMessage(e))
  })

  out
}

pool_ipd_estimate <- function(ipd_list_true){
  stopifnot(length(ipd_list_true) >= 1)

  ipd_all <- rbindlist(lapply(seq_along(ipd_list_true), function(i){
    dt <- as.data.table(ipd_list_true[[i]])
    dt[, trial_id := i]
    dt
  }))

  ipd_all[, arm := factor(arm, levels = c("Control", "Treatment"))]
  ipd_all[, trial_id := factor(trial_id)]

  out <- tryCatch({
    fit <- survival::coxph(survival::Surv(time, status) ~ arm + strata(trial_id), data = ipd_all, ties = "efron")
    cf <- stats::coef(fit)
    vc <- stats::vcov(fit)
    logHR <- unname(cf["armTreatment"])
    se <- sqrt(unname(vc["armTreatment", "armTreatment"]))
    list(ok = TRUE, logHR = logHR, se = se)
  }, error = function(e){
    list(ok = FALSE, logHR = NA_real_, se = NA_real_, err = conditionMessage(e))
  })

  out
}

## ---- 4) One replication for one scenario ----

# simulate_one_rep_meta <- function(
#     K = 10,
#     mu_HR = 0.75,
#     tau_logHR = 0.05,
#     n_control = 250,
#     n_treatment = 250,
#     target_censoring = 0.30,
#     p_miss = 0.6,
#     weib_shape = 1.0,
#     weib_mean_time = 10,
#     digitize_n_grid = 200,
#     digitize_sd_S = 0.01,
#     re_method = "REML",
#     seed = NULL
# ){
#   if (!is.null(seed)) set.seed(seed)
#
#   # trial-specific true effects
#   logHR_k <- stats::rnorm(K, mean = log(mu_HR), sd = tau_logHR)
#   HR_k <- exp(logHR_k)
#
#   # storage
#   per_trial <- vector("list", K)
#   ipd_true_list <- vector("list", K)
#
#   for (k in seq_len(K)) {
#     tr <- simulate_trial_one(
#       hr = HR_k[k],
#       n_control = n_control,
#       n_treatment = n_treatment,
#       shape = weib_shape,
#       mean_time = weib_mean_time,
#       target_censoring = target_censoring,
#       seed = NULL
#     )
#
#     ipd_true_list[[k]] <- tr$ipd
#
#     # "reported" HR from true IPD
#     est_true <- estimate_logHR_cox(tr$ipd)
#
#     # digitize + reconstruct
#     digi <- digitize_km_numeric(tr$ipd, n_grid = digitize_n_grid, sd_S = digitize_sd_S)
#
#     rec <- tryCatch({
#       reconstruct_ipd_IPDfromKM(
#         digi_dt = digi$data,
#         risk_table = tr$risk_table,
#         total_events = NULL,
#         curve_map = c(Control = 1, Treatment = 2)
#       )
#     }, error = function(e){
#       list(ipd = NULL, err = conditionMessage(e))
#     })
#
#     if (is.null(rec$ipd)) {
#       est_rec <- list(ok = FALSE, logHR = NA_real_, se = NA_real_, err = rec$err %||% "reconstruction_failed")
#       rec_ok <- FALSE
#     } else {
#       est_rec <- estimate_logHR_cox(rec$ipd)
#       rec_ok <- isTRUE(est_rec$ok)
#     }
#
#     per_trial[[k]] <- data.table(
#       trial = k,
#       HR_gen = HR_k[k],
#       logHR_gen = logHR_k[k],
#       censor_prop = tr$censor_prop,
#       tot_ev_ctrl = tr$total_events["Control"],
#       tot_ev_trt  = tr$total_events["Treatment"],
#       logHR_true_hat = est_true$logHR,
#       se_true_hat    = est_true$se,
#       true_ok        = isTRUE(est_true$ok),
#       logHR_recon_hat = est_rec$logHR,
#       se_recon_hat    = est_rec$se,
#       recon_ok        = rec_ok
#     )
#   }
#
#   per_study_dt <- rbindlist(per_trial)
#
#   # impose random missingness on reported HR
#   is_missing <- stats::rbinom(K, size = 1, prob = p_miss) == 1
#   per_study_dt[, is_missing := is_missing]
#
#   # Assemble per-method meta inputs
#   # AD-only: keep only studies with reported HR (true IPD estimate)
#   idx_ad <- which(!per_study_dt$is_missing & is.finite(per_study_dt$logHR_true_hat) & is.finite(per_study_dt$se_true_hat))
#
#   # Hybrid: reported uses true HR; missing uses reconstructed HR
#   logHR_hyb <- ifelse(per_study_dt$is_missing, per_study_dt$logHR_recon_hat, per_study_dt$logHR_true_hat)
#   se_hyb    <- ifelse(per_study_dt$is_missing, per_study_dt$se_recon_hat,    per_study_dt$se_true_hat)
#
#   # RM: all reconstructed
#   logHR_rm <- per_study_dt$logHR_recon_hat
#   se_rm    <- per_study_dt$se_recon_hat
#
#   # Pooled IPD (truth)
#   pooled <- pool_ipd_estimate(ipd_true_list)
#
#   # Fit meta models
#   fit_ad  <- meta_fit_rma(per_study_dt$logHR_true_hat[idx_ad], per_study_dt$se_true_hat[idx_ad], re_method = re_method)
#   fit_hyb <- meta_fit_rma(logHR_hyb, se_hyb, re_method = re_method)
#   fit_rm  <- meta_fit_rma(logHR_rm,  se_rm,  re_method = re_method)
#
#   per_method_dt <- rbindlist(list(
#     data.table(method = "AD_only",   ok = fit_ad$ok,  logHR_hat = fit_ad$logHR,  se_hat = fit_ad$se,  tau2_hat = fit_ad$tau2,  I2_hat = fit_ad$I2,  k_used = fit_ad$k),
#     data.table(method = "Hybrid",    ok = fit_hyb$ok, logHR_hat = fit_hyb$logHR, se_hat = fit_hyb$se, tau2_hat = fit_hyb$tau2, I2_hat = fit_hyb$I2, k_used = fit_hyb$k),
#     data.table(method = "RM",        ok = fit_rm$ok,  logHR_hat = fit_rm$logHR,  se_hat = fit_rm$se,  tau2_hat = fit_rm$tau2,  I2_hat = fit_rm$I2,  k_used = fit_rm$k),
#     data.table(method = "Pooled_IPD", ok = pooled$ok, logHR_hat = pooled$logHR, se_hat = pooled$se, tau2_hat = NA_real_, I2_hat = NA_real_, k_used = K)
#   ), use.names = TRUE)
#
#   per_method_dt[, HR_hat := exp(logHR_hat)]
#
#   list(
#     per_study_dt = per_study_dt,
#     per_method_dt = per_method_dt
#   )
# }



simulate_one_rep_meta <- function(
    K = 10,
    mu_HR = 0.75,
    tau_logHR = 0.05,
    n_control = 250,
    n_treatment = 250,
    target_censoring = 0.30,
    p_miss = 0.6,
    weib_shape = 1.0,
    weib_mean_time = 10,
    digitize_n_grid = 200,
    digitize_sd_S = 0.04,
    re_method = "REML",
    seed = NULL
){
  if (!is.null(seed)) set.seed(seed)

  # trial-specific true effects
  logHR_k <- stats::rnorm(K, mean = log(mu_HR), sd = tau_logHR)
  HR_k <- exp(logHR_k)

  # storage
  per_trial <- vector("list", K)
  ipd_true_list <- vector("list", K)

  for (k in seq_len(K)) {
    trial_seed <- if (!is.null(seed)) seed + k else NULL   # 为每个试验设置不同种子
    tr <- simulate_trial_one(
      hr = HR_k[k],
      n_control = n_control,
      n_treatment = n_treatment,
      shape = weib_shape,
      mean_time = weib_mean_time,
      target_censoring = target_censoring,
      seed = trial_seed    # 传入种子
    )

    ipd_true_list[[k]] <- tr$ipd

    # "reported" HR from true IPD
    est_true <- estimate_logHR_cox(tr$ipd)

    # digitize + reconstruct
    digi <- digitize_km_numeric(tr$ipd, tr$risk_table)

    rec <- tryCatch({
      reconstruct_ipd_IPDfromKM(
        digi_dt = digi$data,
        risk_table = tr$risk_table,
        total_events = NULL,
        curve_map = c(Control = 1, Treatment = 2)
      )
    }, error = function(e){
      list(ipd = NULL, err = conditionMessage(e))
    })

    if (is.null(rec$ipd)) {
      est_rec <- list(ok = FALSE, logHR = NA_real_, se = NA_real_, err = rec$err %||% "reconstruction_failed")
      rec_ok <- FALSE
    } else {
      est_rec <- estimate_logHR_cox(rec$ipd)
      rec_ok <- isTRUE(est_rec$ok)
    }

    per_trial[[k]] <- data.table(
      trial = k,
      HR_gen = HR_k[k],
      logHR_gen = logHR_k[k],
      censor_prop = tr$censor_prop,
      tot_ev_ctrl = tr$total_events["Control"],
      tot_ev_trt  = tr$total_events["Treatment"],
      logHR_true_hat = est_true$logHR,
      se_true_hat    = est_true$se,
      true_ok        = isTRUE(est_true$ok),
      logHR_recon_hat = est_rec$logHR,
      se_recon_hat    = est_rec$se,
      recon_ok        = rec_ok
    )
  }

  per_study_dt <- rbindlist(per_trial)

  # impose random missingness on reported HR
  is_missing <- stats::rbinom(K, size = 1, prob = p_miss) == 1
  per_study_dt[, is_missing := is_missing]

  # Assemble per-method meta inputs
  # AD-only: keep only studies with reported HR (true IPD estimate)
  idx_ad <- which(!per_study_dt$is_missing & is.finite(per_study_dt$logHR_true_hat) & is.finite(per_study_dt$se_true_hat))

  # Hybrid: reported uses true HR; missing uses reconstructed HR
  # 在 per_study_dt 中已有 recon_ok 列
  # 构造Hybrid的logHR和se：对于缺失的试验，仅当重建成功时才使用重建值，否则仍为缺失
  logHR_hyb <- per_study_dt$logHR_true_hat   # 默认用真值
  se_hyb    <- per_study_dt$se_true_hat

  # 对于缺失的试验
  miss_idx <- which(per_study_dt$is_missing)
  for (i in miss_idx) {
    if (per_study_dt$recon_ok[i]) {
      logHR_hyb[i] <- per_study_dt$logHR_recon_hat[i]
      se_hyb[i]    <- per_study_dt$se_recon_hat[i]
    } else {
      logHR_hyb[i] <- NA_real_   # 重建失败，视为缺失
      se_hyb[i]    <- NA_real_
    }
  }

  # RM方法：仅使用重建成功的试验（重建失败的设为NA）
  logHR_rm <- ifelse(per_study_dt$recon_ok, per_study_dt$logHR_recon_hat, NA_real_)
  se_rm    <- ifelse(per_study_dt$recon_ok, per_study_dt$se_recon_hat, NA_real_)

  # Pooled IPD (truth)
  pooled <- pool_ipd_estimate(ipd_true_list)

  # Fit meta models
  fit_ad  <- meta_fit_rma(per_study_dt$logHR_true_hat[idx_ad], per_study_dt$se_true_hat[idx_ad], re_method = re_method)
  fit_hyb <- meta_fit_rma(logHR_hyb, se_hyb, re_method = re_method)
  fit_rm  <- meta_fit_rma(logHR_rm,  se_rm,  re_method = re_method)

  per_method_dt <- rbindlist(list(
    data.table(method = "AD_only",   ok = fit_ad$ok,  logHR_hat = fit_ad$logHR,  se_hat = fit_ad$se,  tau2_hat = fit_ad$tau2,  I2_hat = fit_ad$I2,  k_used = fit_ad$k),
    data.table(method = "Hybrid",    ok = fit_hyb$ok, logHR_hat = fit_hyb$logHR, se_hat = fit_hyb$se, tau2_hat = fit_hyb$tau2, I2_hat = fit_hyb$I2, k_used = fit_hyb$k),
    data.table(method = "RM",        ok = fit_rm$ok,  logHR_hat = fit_rm$logHR,  se_hat = fit_rm$se,  tau2_hat = fit_rm$tau2,  I2_hat = fit_rm$I2,  k_used = fit_rm$k),
    data.table(method = "Pooled_IPD", ok = pooled$ok, logHR_hat = pooled$logHR, se_hat = pooled$se, tau2_hat = NA_real_, I2_hat = NA_real_, k_used = K)
  ), use.names = TRUE)

  per_method_dt[, HR_hat := exp(logHR_hat)]

  list(
    per_study_dt = per_study_dt,
    per_method_dt = per_method_dt
  )
}

## ---- 5) Run the full grid ----

build_scenario_grid <- function(
    p_miss_vec = c(0.3, 0.6, 0.9),
    tau_levels = list(low = 0.05, high = 0.25),
    n_levels = list(standard = 250),
    K = 10
){
  grid <- CJ(
    p_miss = as.numeric(p_miss_vec),
    tau_label = names(tau_levels),
    n_label = names(n_levels)
  )

  grid[, tau_logHR := unlist(tau_levels)[match(tau_label, names(tau_levels))]]
  grid[, n_per_arm := unlist(n_levels)[match(n_label, names(n_levels))]]
  grid[, K := as.integer(K)]

  grid[, scenario_id := sprintf("n%d_K%d_miss%02d_tau%s",
                               n_per_arm,
                               K,
                               as.integer(round(100 * p_miss)),
                               ifelse(tau_label == "low", "L", "H"))]

  grid[]
}

run_meta_simulation_paired <- function(grid, n_rep=200, seed_base=20360220L,
                                       mu_HR=0.75, digitize_n_grid=200, digitize_sd_S=0.01) {
  library(data.table)

  grid <- as.data.table(grid)
  res_methods_all <- list()
  res_studies_all <- list()
  idx <- 1L

  # 把 grid 拆成“基础世界”：不包含 p_miss，只包含 tau/n/K
  base_worlds <- unique(grid[, .(tau_label, tau_logHR, n_label, n_per_arm, K)])

  for (b in seq_len(nrow(base_worlds))) {
    bw <- base_worlds[b]
    tau_label <- bw$tau_label
    tau_logHR <- bw$tau_logHR
    n_per_arm <- bw$n_per_arm
    K <- bw$K
    n_label <- bw$n_label

    for (r in seq_len(n_rep)) {
      # 固定一个 base world + rep 的种子 -> 固定一套 trials + digitize + recon
      seed_world <- seed_base + 100000L*b + r
      one <- simulate_one_rep_meta(
        K=K, mu_HR=mu_HR, tau_logHR=tau_logHR,
        n_control=n_per_arm, n_treatment=n_per_arm,
        p_miss=0.0,  # 先生成完整世界（missing 暂时不用）
        digitize_n_grid=digitize_n_grid, digitize_sd_S=digitize_sd_S,
        seed=seed_world
      )

      # trial-level 固定
      stud <- as.data.table(one$per_study_dt)
      stud[, `:=`(tau_label=tau_label, tau_logHR=tau_logHR, n_label=n_label, n_per_arm=n_per_arm, K=K, rep=r, seed=seed_world)]

      # 固定一组 u_k，让 missing 集合随 p_miss 嵌套
      set.seed(seed_world + 9999L)
      u <- runif(K)
      stud[, u_missing := u]

      # 对 grid 里所有 p_miss 做派生
      pm_values <- sort(unique(grid$p_miss))
      for (pm in pm_values) {
        stud_pm <- copy(stud)
        stud_pm[, `:=`(p_miss=pm, is_missing = (u_missing < pm))]

        mu_logHR <- log(mu_HR)

        # ---- AD-only: only non-missing TRUE HR ----
        ad <- stud_pm[is_missing==FALSE & true_ok==TRUE & is.finite(logHR_true_hat) & is.finite(se_true_hat) & se_true_hat>0]
        fit_ad <- meta_fit_rma(ad$logHR_true_hat, ad$se_true_hat)
        ad_row <- data.table(
          method = "AD_only",
          ok = fit_ad$ok,
          logHR_hat = fit_ad$logHR,
          se_hat = fit_ad$se,
          tau2_hat = fit_ad$tau2,
          I2_hat = fit_ad$I2,
          k_used = fit_ad$k
        )

        # ---- Hybrid: non-missing TRUE HR + missing RECON HR ----
        hy <- rbind(
          stud_pm[is_missing==FALSE & true_ok==TRUE, .(logHR=logHR_true_hat, se=se_true_hat)],
          stud_pm[is_missing==TRUE  & recon_ok==TRUE, .(logHR=logHR_recon_hat, se=se_recon_hat)]
        )
        fit_hy <- meta_fit_rma(hy$logHR, hy$se)
        hy_row <- data.table(
          method = "Hybrid",
          ok = fit_hy$ok,
          logHR_hat = fit_hy$logHR,
          se_hat = fit_hy$se,
          tau2_hat = fit_hy$tau2,
          I2_hat = fit_hy$I2,
          k_used = fit_hy$k
        )

        # ---- RM: all recon (should be constant across pm) ----
        rm <- stud_pm[recon_ok==TRUE & is.finite(logHR_recon_hat) & is.finite(se_recon_hat) & se_recon_hat>0]
        fit_rm <- meta_fit_rma(rm$logHR_recon_hat, rm$se_recon_hat)
        rm_row <- data.table(
          method = "RM",
          ok = fit_rm$ok,
          logHR_hat = fit_rm$logHR,
          se_hat = fit_rm$se,
          tau2_hat = fit_rm$tau2,
          I2_hat = fit_rm$I2,
          k_used = fit_rm$k
        )

        # ---- Pooled IPD: reuse from one (already constant) ----
        pool_row <- one$per_method_dt[method=="Pooled_IPD", .(
          method, ok, logHR_hat, se_hat,
          tau2_hat = NA_real_, I2_hat = NA_real_, k_used = K
        )]

        methods_pm <- rbindlist(list(ad_row, hy_row, rm_row, pool_row), fill=TRUE)

        methods_pm[, `:=`(
          scenario_id = sprintf("n%d_K%d_miss%02d_tau%s", n_per_arm, K, as.integer(pm*100), ifelse(tau_label=="low","L","H")),
          p_miss = pm, tau_label=tau_label, tau_logHR=tau_logHR,
          n_label=n_label, n_per_arm=n_per_arm, K=K, rep=r, seed=seed_world
        )]
        methods_pm[, err := logHR_hat - mu_logHR]
        methods_pm[, HR_hat := exp(logHR_hat)]

        res_methods_all[[idx]] <- methods_pm
        res_studies_all[[idx]] <- stud_pm
        idx <- idx + 1L
      }
    }
  }

  list(
    res_methods_long = rbindlist(res_methods_all, fill=TRUE),
    res_studies_long = rbindlist(res_studies_all, fill=TRUE)
  )
}



## ---- 6) Summaries + plotting ----
summarise_performance <- function(res_methods_long, mu_HR = 0.75, alpha = 0.05){
  stopifnot(is.data.table(res_methods_long))
  target <- log(mu_HR)
  z <- qnorm(1 - alpha/2)

  res_methods_long[, err := logHR_hat - target]
  res_methods_long[, ci_lower := logHR_hat - z * se_hat]
  res_methods_long[, ci_upper := logHR_hat + z * se_hat]
  res_methods_long[, cover := (ci_lower <= target) & (ci_upper >= target)]

  perf <- res_methods_long[, .(
    n_rep = .N,
    fail_rate = mean(!is.finite(logHR_hat)),
    mean_k_used = mean(k_used, na.rm = TRUE),
    bias = mean(err, na.rm = TRUE),
    abs_bias = mean(abs(err), na.rm = TRUE),
    rmse = sqrt(mean(err^2, na.rm = TRUE)),
    coverage = mean(cover, na.rm = TRUE)    # 新增
  ), by = .(scenario_id, p_miss, tau_label, tau_logHR, n_label, n_per_arm, K, method)]

  perf[]
}

# summarise_performance <- function(res_methods_long, mu_HR = 0.75){
#   stopifnot(is.data.table(res_methods_long))
#   target <- log(mu_HR)
#
#   res_methods_long[, err := logHR_hat - target]
#
#   perf <- res_methods_long[, .(
#     n_rep = .N,
#     fail_rate = mean(!is.finite(logHR_hat)),
#     mean_k_used = mean(k_used, na.rm = TRUE),
#     bias = mean(err, na.rm = TRUE),
#     abs_bias = mean(abs(err), na.rm = TRUE),
#     rmse = sqrt(mean(err^2, na.rm = TRUE))
#   ), by = .(scenario_id, p_miss, tau_label, tau_logHR, n_label, n_per_arm, K, method)]
#
#   perf[]
# }

plot_bias_rmse <- function(perf_dt, out_dir = NULL){
  stopifnot(is.data.table(perf_dt))

  # bias plot
  p_bias <- ggplot(perf_dt, aes(x = p_miss, y = bias, color = method, group = method)) +
    geom_line() +
    geom_point(size = 2) +
    facet_grid(tau_label ~ n_label, scales = "free_y") +
    labs(x = "Missing HR proportion", y = "Bias (logHR scale)", title = "Meta-estimator bias vs missingness")

  # rmse plot
  p_rmse <- ggplot(perf_dt, aes(x = p_miss, y = rmse, color = method, group = method)) +
    geom_line() +
    geom_point(size = 2) +
    facet_grid(tau_label ~ n_label, scales = "free_y") +
    labs(x = "Missing HR proportion", y = "RMSE (logHR scale)", title = "Meta-estimator RMSE vs missingness")

  if (!is.null(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = file.path(out_dir, "bias_by_missing.png"), plot = p_bias, width = 9, height = 6, dpi = 200)
    ggsave(filename = file.path(out_dir, "rmse_by_missing.png"), plot = p_rmse, width = 9, height = 6, dpi = 200)
  }

  list(bias_plot = p_bias, rmse_plot = p_rmse)
}

## ---- 7) Main (runs only when executed via Rscript) ----

if (sys.nframe() == 0) {
  # ---- user-tunable knobs ----
  settings <- list(
    mu_HR = 0.75,
    K = 10,
    n_rep = 200,
    target_censoring = 0.30,
    p_miss_vec = c(0, 0.3, 0.6, 0.9),
    tau_levels = list(low = 0.05, high = 0.25),
    n_levels = list(standard = 250),
    weib_shape = 1.0,
    weib_mean_time = 10,
    digitize_n_grid = 200,
    digitize_sd_S = 0.01,
    re_method = "REML",
    seed_base = 20260219
  )

  # Output folder
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(getwd(), paste0("meta_sim_v3_output_", stamp))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Build grid
  grid <- build_scenario_grid(
    p_miss_vec = settings$p_miss_vec,
    tau_levels = settings$tau_levels,
    n_levels = settings$n_levels,
    K = settings$K
  )

  # Run sim
  sim <- run_meta_simulation(
    scenario_grid = grid,
    n_rep = settings$n_rep,
    seed_base = settings$seed_base,
    mu_HR = settings$mu_HR,
    target_censoring = settings$target_censoring,
    weib_shape = settings$weib_shape,
    weib_mean_time = settings$weib_mean_time,
    digitize_n_grid = settings$digitize_n_grid,
    digitize_sd_S = settings$digitize_sd_S,
    re_method = settings$re_method,
    verbose = TRUE
  )

  # Save raw outputs
  fwrite(sim$res_methods_long, file.path(out_dir, "res_methods_long.csv"))
  fwrite(sim$res_studies_long, file.path(out_dir, "res_studies_long.csv"))

  # Summaries
  perf <- summarise_performance(sim$res_methods_long, mu_HR = settings$mu_HR)
  fwrite(perf, file.path(out_dir, "perf_summary.csv"))

  # Plots
  plt <- plot_bias_rmse(perf, out_dir = out_dir)

  # Print a quick console summary
  message("\nDone. Outputs saved to: ", out_dir)
  print(perf[order(tau_label, n_label, method, p_miss)])
}
