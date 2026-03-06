#==============================================================#
# Digitize KM plot using SurvdigitizeR
#==============================================================#

# Digitize a KM plot image into (time, survival) points.
#
# Returns a list with:
#   $data   : data.frame/data.table with columns time, St, curve (+ optional arm)
#   $plot   : path to saved QA PNG (or NA)
#   $meta   : key digitization arguments

digitize_km_event_times <- function(
    ipd_true,
    sd_S = 0.04,
    seed = NULL,
    risk_table = NULL,      # NEW: optional risk table; if provided, we include its time grid in digitized points
    sd_t = 0.00,            # NEW: optional time-axis digitization jitter (set >0 to stress-test time calibration)
    p_mark = 1.0,           # NEW: proportion of true censorings that are "visibly marked" on the plot
    sd_mark_t = 0.00,       # NEW: jitter for censor-mark times (set >0 to mimic reading marks from a figure)
    mark_merge_eps = 0.00   # NEW: merge marks closer than this (0 = no merging)
) {
  stopifnot(is.data.frame(ipd_true))
  required_cols <- c("time", "status", "arm")
  stopifnot(all(required_cols %in% names(ipd_true)))

  ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))

  if (!is.null(seed)) set.seed(seed)

  # ---- helper: build "digitized" (time, St) for one arm
  # We mimic a digitizer by sampling the *true* KM step function on a chosen time grid,
  # then adding noise to St (and optionally to time), and finally enforcing monotonicity.
  km_digitize_onearm <- function(arm_label) {
    dat <- ipd_true[ipd_true$arm == arm_label, c("time", "status")]
    fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)

    # 1) base grid: drop-bottom times (event times where n.event>0)
    event_times <- fit$time[fit$n.event > 0]

    # 2) ALSO include risk-table time grid (requested): ensures digitized points include all risk-table times
    rt_times <- numeric(0)
    if (!is.null(risk_table) && is.data.frame(risk_table) && "time" %in% names(risk_table)) {
      rt_times <- as.numeric(risk_table$time)
      rt_times <- rt_times[is.finite(rt_times)]
    }

    # 3) union + ensure includes time 0
    times_grid <- sort(unique(c(0, event_times, rt_times)))
    times_grid <- times_grid[is.finite(times_grid) & times_grid >= 0]

    # Evaluate the true KM step function at the grid times
    S_true <- summary(fit, times = times_grid, extend = TRUE)$surv
    S_true[1] <- 1  # enforce S(0)=1

    # Add noise to time (optional) and survival
    if (sd_t > 0) {
      times_grid <- pmax(0, times_grid + rnorm(length(times_grid), 0, sd_t))
    }

    S_noisy <- S_true + rnorm(length(S_true), 0, sd_S)
    S_noisy[1] <- 1

    # Clamp and enforce non-increasing survival
    S_noisy <- pmax(pmin(S_noisy, 1), 0)

    # Order by (possibly jittered) time; keep the *minimum* St per time to represent bottoms of drops
    ord <- order(times_grid)
    times_grid <- times_grid[ord]
    S_noisy <- S_noisy[ord]

    # Collapse duplicated times by taking MIN St (bottom of any vertical drop)
    tmp <- data.frame(time = times_grid, St = S_noisy)
    tmp <- aggregate(St ~ time, data = tmp, FUN = min)
    tmp <- tmp[order(tmp$time), ]
    tmp$St <- cummin(tmp$St)

    tmp
  }

  # ---- helper: generate "marked censoring times" (KMtoIPD expects a numeric vector of times)
  censor_marks_onearm <- function(arm_label) {
    dat <- ipd_true[ipd_true$arm == arm_label, c("time", "status")]
    cens_times <- dat$time[dat$status == 0]
    cens_times <- as.numeric(cens_times[is.finite(cens_times) & cens_times >= 0])

    if (length(cens_times) == 0) return(numeric(0))

    # Subsample to mimic incomplete visibility of censor marks on published figures
    p_mark_use <- max(min(as.numeric(p_mark)[1], 1), 0)
    if (p_mark_use < 1) {
      keep <- stats::runif(length(cens_times)) < p_mark_use
      cens_times <- cens_times[keep]
    }

    # Add small jitter (optional) to mimic reading marks from an image
    if (sd_mark_t > 0 && length(cens_times) > 0) {
      cens_times <- pmax(0, cens_times + rnorm(length(cens_times), 0, sd_mark_t))
    }

    cens_times <- sort(cens_times)

    # Optionally merge very-close marks (to mimic overplotting where multiple marks look like one)
    if (mark_merge_eps > 0 && length(cens_times) > 1) {
      out <- cens_times[1]
      for (i in 2:length(cens_times)) {
        if (cens_times[i] - out[length(out)] > mark_merge_eps) out <- c(out, cens_times[i])
      }
      cens_times <- out
    }

    cens_times
  }

  dt_control <- km_digitize_onearm("Control")
  dt_treat   <- km_digitize_onearm("Treatment")

  dt_control$curve <- 1L   # Control 对应 curve 1
  dt_treat$curve   <- 2L   # Treatment 对应 curve 2

  digi_dt <- rbind(dt_control, dt_treat)
  digi_dt <- digi_dt[order(digi_dt$curve, digi_dt$time), ]

  # --- censor marks (numeric vector per arm; KMtoIPD format)
  cens_mark <- list(
    Control = censor_marks_onearm("Control"),
    Treatment = censor_marks_onearm("Treatment")
  )

  # quality check
  stopifnot(!anyNA(digi_dt$St))
  stopifnot(all(is.finite(digi_dt$time)))
  stopifnot(all(digi_dt$St >= 0 & digi_dt$St <= 1))

  list(
    data = digi_dt,
    cens_mark = cens_mark,          # NEW: expose censor marks for KMtoIPD reconstruction
    plot = NA_character_,
    meta = list(
      method = "event_times_plus_risktable_numeric",  # UPDATED label
      sd_S = sd_S,
      sd_t = sd_t,
      p_mark = p_mark,
      sd_mark_t = sd_mark_t,
      mark_merge_eps = mark_merge_eps,
      seed = seed
    )
  )
}

# Grid base ----
# digitize_km_event_times2 <- function(
#     ipd_true,
#     sd_S = 0.04,
#     sd_t = 0.01,
#     n_grid = 200,
#     merge_tol = 0,
#     seed = NULL
# ) {
#   stopifnot(is.data.frame(ipd_true))
#   required_cols <- c("time", "status", "arm")
#   stopifnot(all(required_cols %in% names(ipd_true)))
#
#   ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))
#
#   if (!is.null(seed)) set.seed(seed)
#
#   # ------------------------------------------------------------------
#   # v2 digitization (numeric surrogate):
#   #   - generate *dense* (time, S(t)) points on a grid (more like SurvDigitizeR)
#   #   - add noise to BOTH time and survival
#   #   - always output censor marks (per arm), with optional overlap merging
#   # Rationale:
#   #   - KMtoIPD requires drop-bottom coordinates, which we will later *extract*
#   #     from these dense noisy points (closer to real digitization).
#   #   - other reconstruction methods also benefit from a dense digitized curve.
#   # ------------------------------------------------------------------
#
#   # Determine a shared administrative end time (used for clamping noisy times)
#   t_end <- max(ipd_true$time, na.rm = TRUE)
#   stopifnot(is.finite(t_end) && t_end > 0)
#
#   # If sd_t is not supplied, pick an "auto" value tied to the grid spacing.
#   # NOTE: this ensures every simulation has *some* time-axis digitization noise.
#   if (is.null(sd_t)) {
#     # comment: auto time noise ≈ 1/3 of average grid spacing
#     sd_t <- (t_end / max(as.integer(n_grid), 2L)) / 3
#   }
#   sd_t <- as.numeric(sd_t)[1]
#   sd_S <- as.numeric(sd_S)[1]
#   n_grid <- as.integer(n_grid)[1]
#   merge_tol <- as.numeric(merge_tol)[1]
#
#   stopifnot(is.finite(sd_t) && sd_t >= 0)
#   stopifnot(is.finite(sd_S) && sd_S >= 0)
#   stopifnot(is.finite(merge_tol) && merge_tol >= 0)
#   stopifnot(n_grid >= 10)
#
#   if (sd_t == 0 && sd_S == 0) {
#     stop("digitize_km_event_times(): sd_t and sd_S cannot both be 0 (each simulation must include noise).")
#   }
#
#   # helper: optionally merge censor marks that are too close (simulate overlap)
#   merge_close_times <- function(x, tol) {
#     x <- sort(x)
#     if (length(x) <= 1 || tol <= 0) return(x)
#     keep <- c(TRUE, diff(x) > tol)
#     x[keep]
#   }
#
#   km_dense_noisy <- function(arm_label) {
#     dat <- ipd_true[ipd_true$arm == arm_label, c("time", "status")]
#     fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
#
#     # --- construct a dense time grid (more like SurvDigitizeR output)
#     grid_base <- seq(0, t_end, length.out = n_grid)
#     # Include true event times so that drops remain observable even after thinning
#     ev_times <- fit$time[fit$n.event > 0]
#     times <- sort(unique(c(0, grid_base, ev_times)))
#
#     # --- KM survival on the grid (right-continuous step function)
#     S_grid <- summary(fit, times = times, extend = TRUE)$surv
#     # Ensure start at (0, 1)
#     if (length(times) == 0 || times[1] != 0) {
#       times <- c(0, times)
#       S_grid <- c(1, S_grid)
#     } else {
#       S_grid[1] <- 1
#     }
#
#     # --- add digitization noise to time and survival
#     times_noisy <- times + rnorm(length(times), 0, sd_t)
#     # clamp to plotting window
#     times_noisy <- pmin(pmax(times_noisy, 0), t_end)
#     # enforce time order after jitter
#     ord <- order(times_noisy)
#     times_noisy <- times_noisy[ord]
#     S_grid <- S_grid[ord]
#
#     S_noisy <- S_grid + rnorm(length(S_grid), 0, sd_S)
#     # clamp & force monotone non-increasing (stabilize noisy digitization)
#     S_noisy <- pmax(pmin(S_noisy, 1), 0)
#     S_noisy[1] <- 1
#
#     # S_noisy <- cummin(S_noisy)
#     # New
#     # IMPORTANT FIX: avoid downward bias from cummin()
#     # Use isotonic regression to enforce monotone non-increasing with minimal L2 change
#     iso <- stats::isoreg(seq_along(S_noisy), -S_noisy)   # increasing fit on (-S)
#     S_noisy <- -iso$yf
#     S_noisy[1] <- 1
#     S_noisy <- pmax(pmin(S_noisy, 1), 0)
#
#
#     data.frame(time = times_noisy, St = S_noisy)
#   }
#
#   # --- digitized dense curve points
#   dt_control <- km_dense_noisy("Control")
#   dt_treat   <- km_dense_noisy("Treatment")
#
#   dt_control$curve <- 1L   # Control 对应 curve 1
#   dt_treat$curve   <- 2L   # Treatment 对应 curve 2
#
#   digi_dt <- rbind(dt_control, dt_treat)
#   digi_dt <- digi_dt[order(digi_dt$curve, digi_dt$time), ]
#
#   # --- censor marks (always output)
#   # NOTE: these are *marked* censoring times, which are key optional inputs for KMtoIPD.
#   get_censor_marks <- function(arm_label) {
#     dat <- ipd_true[ipd_true$arm == arm_label, c("time", "status")]
#     cens <- dat$time[dat$status == 0]
#     if (length(cens) == 0) return(numeric(0))
#     cens_noisy <- cens + rnorm(length(cens), 0, sd_t)  # use same time-axis noise
#     cens_noisy <- pmin(pmax(cens_noisy, 0), t_end)
#     cens_noisy <- merge_close_times(sort(cens_noisy), tol = merge_tol)
#     cens_noisy
#   }
#   censor_marks <- list(
#     Control = get_censor_marks("Control"),
#     Treatment = get_censor_marks("Treatment")
#   )
#   # also store by curve id for convenience
#   censor_marks_by_curve <- list(
#     `1` = censor_marks$Control,
#     `2` = censor_marks$Treatment
#   )
#
#   # 最终质量检查
#   stopifnot(!anyNA(digi_dt$St))
#   stopifnot(all(digi_dt$St >= 0 & digi_dt$St <= 1))
#
#   list(
#     data = digi_dt,
#     plot = NA_character_,
#     censor_marks = censor_marks,              # <-- NEW: per-arm censor marks
#     censor_marks_by_curve = censor_marks_by_curve,  # <-- NEW: per-curve censor marks
#     meta = list(
#       method = "dense_grid_numeric_v2",
#       sd_S = sd_S,
#       sd_t = sd_t,
#       n_grid = n_grid,
#       merge_tol = merge_tol,
#       t_end = t_end,
#       n_pts_control = nrow(dt_control),
#       n_pts_treatment = nrow(dt_treat),
#       n_cens_control = length(censor_marks$Control),
#       n_cens_treatment = length(censor_marks$Treatment),
#       seed = seed
#     )
#   )
# }

# Event base----
# digitize_km_event_times_v1 <- function(
#     ipd_true,
#     sd_S = 0.04,
#     seed = NULL
# ) {
#   stopifnot(is.data.frame(ipd_true))
#   required_cols <- c("time", "status", "arm")
#   stopifnot(all(required_cols %in% names(ipd_true)))
#
#   ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))
#
#   if (!is.null(seed)) set.seed(seed)
#
#   km_event_noisy <- function(arm_label) {
#     dat <- ipd_true[ipd_true$arm == arm_label, c("time", "status")]
#     fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
#
#     # 提取事件时间点和对应的生存概率
#     event_times <- fit$time[fit$n.event > 0]
#     S_at_events <- summary(fit, times = event_times)$surv
#
#     # 确保包含时间 0 且 S=1
#     if (event_times[1] > 0) {
#       event_times <- c(0, event_times)
#       S_at_events <- c(1, S_at_events)
#     }
#
#     # 添加噪声
#     S_noisy <- S_at_events + rnorm(length(S_at_events), 0, sd_S)
#     S_noisy[1] <- 1  # 时间 0 强制为 1
#
#     # 限制在 [0,1] 并强制非增
#     S_noisy <- pmax(pmin(S_noisy, 1), 0)
#     for (i in length(S_noisy):2) {
#       if (S_noisy[i] > S_noisy[i-1]) S_noisy[i] <- S_noisy[i-1]
#     }
#
#     data.frame(time = event_times, St = S_noisy)
#   }
#
#   dt_control <- km_event_noisy("Control")
#   dt_treat   <- km_event_noisy("Treatment")
#
#   dt_control$curve <- 1L   # Control 对应 curve 1
#   dt_treat$curve   <- 2L   # Treatment 对应 curve 2
#
#   digi_dt <- rbind(dt_control, dt_treat)
#   digi_dt <- digi_dt[order(digi_dt$curve, digi_dt$time), ]
#
#   # 最终质量检查
#   stopifnot(!anyNA(digi_dt$St))
#   stopifnot(all(digi_dt$St >= 0 & digi_dt$St <= 1))
#
#   list(
#     data = digi_dt,
#     plot = NA_character_,
#     meta = list(
#       method = "event_times_numeric",
#       sd_S = sd_S,
#       n_control = nrow(dt_control),
#       n_treatment = nrow(dt_treat),
#       seed = seed
#     )
#   )
# }

# digitize_km_image ----
digitize_km_image <- function(
    img_path,
    x_end,
    mode = c("remove_censor_marks", "keep_censor_marks"),
    num_curves = 2,
    x_increment = NULL,
    y_increment = 0.2,
    bg_lightness = 0.4,
    attempt_OCR = FALSE,
    word_sensitivity = 30,
    y_text_vertical = FALSE,
    nr_neighbors = 20,
    enhance = FALSE,
    impute_size = 0,
    line_censoring = FALSE,
    seed = 2025,
    make_plot = TRUE,
    plot_out = sub("\\.png$", "_digitized_points.png", img_path),
    # Optional: provide truth IPD for deterministic arm alignment in simulation
    truth_ipd = NULL,
    arm_var = "arm",
    arm_levels = c("Control", "Treatment"),
    # ... existing args ...,
    enforce_hr_lt1 = TRUE,   # set TRUE in simulation where true HR < 1 always
    swap_tol = 0.02,          # tolerance to avoid swapping when effect ~ 0
    swap_grid_n = 200,        # grid points for estimating logHR direction
    swap_eps = 1e-6           # numerical safeguard for log()
){
  mode <- match.arg(mode)

  if (!requireNamespace("SurvdigitizeR", quietly = TRUE)) {
    stop("Package 'SurvdigitizeR' is required for digitize_km_image().")
  }

  stopifnot(is.character(img_path), length(img_path) == 1L, file.exists(img_path))
  stopifnot(is.finite(x_end), x_end > 0)

  set.seed(seed)

  if (is.null(x_increment)) x_increment <- x_end / 8
  x_increment <- as.numeric(x_increment)[1]
  y_increment <- as.numeric(y_increment)[1]

  censoring_flag <- (mode == "remove_censor_marks")

  # 初始化 meta 列表
  meta <- list(
    img_path = img_path,
    x_end = x_end,
    x_increment = x_increment,
    y_increment = y_increment,
    mode = mode,
    num_curves = num_curves,
    bg_lightness = bg_lightness,
    attempt_OCR = attempt_OCR,
    word_sensitivity = word_sensitivity,
    y_text_vertical = y_text_vertical,
    nr_neighbors = nr_neighbors,
    enhance = enhance,
    impute_size = impute_size,
    line_censoring = line_censoring,
    seed = seed
  )

  out <- tryCatch(
    SurvdigitizeR::survival_digitize(
      img_path = img_path,
      bg_lightness = as.numeric(bg_lightness)[1],
      attempt_OCR = as.logical(attempt_OCR)[1],
      word_sensitivity = as.numeric(word_sensitivity)[1],
      num_curves = as.integer(num_curves)[1],
      censoring = as.logical(censoring_flag)[1],
      x_start = 0,
      x_end = as.numeric(x_end)[1],
      x_increment = x_increment,
      y_start = 0,
      y_end = 1,
      y_increment = y_increment,
      y_text_vertical = as.logical(y_text_vertical)[1],
      nr_neighbors = as.integer(nr_neighbors)[1],
      enhance = as.logical(enhance)[1],
      impute_size = as.integer(impute_size)[1],
      line_censoring = as.logical(line_censoring)[1]
    ),
    error = function(e){
      warning(sprintf("[SurvdigitizeR failed] %s", e$message))
      return(NULL)
    }
  )
  if (is.null(out)) return(NULL)

  # 统一转换为 data.table（你已加载 data.table，直接使用）
  dt <- data.table::as.data.table(out)
  if ("times" %in% names(dt)) data.table::setnames(dt, "times", "time")
  if (!all(c("time", "St", "curve") %in% names(dt))) return(NULL)

  # 基本清理
  dt <- dt[is.finite(time) & is.finite(St)]
  dt <- dt[time >= 0 & time <= x_end]
  dt[, St := pmax(pmin(St, 1), 0)]
  data.table::setorder(dt, curve, time)
  dt[, St := cummin(St), by = curve]

  # 可选：强制曲线 1 = 对照，曲线 2 = 治疗（当真实 HR < 1 时）
  if (isTRUE(enforce_hr_lt1) && "curve" %in% names(dt) &&
      length(unique(dt$curve)) >= 2) {

    # 局部函数：从数字化 KM 估计 logHR 方向
    estimate_logHR_from_KM <- function(dt, ctrl_curve = 1, trt_curve = 2,
                                       grid_n = 200, eps = 1e-6) {
      prep_one <- function(df) {
        df <- df[is.finite(df$time) & is.finite(df$St), .(time, St)]
        df <- df[order(time)]
        # 去重（保留每个时间的最后一个 St）
        df <- df[, .(St = St[.N]), by = time]
        # 添加 (0,1) 如果缺失
        if (!any(df$time == 0)) df <- rbind(data.table(time = 0, St = 1), df)
        df <- df[order(time)]
        df[, St := pmin(1, pmax(0, St))]
        df[, St := cummin(St)]
        df
      }

      cdf <- prep_one(dt[curve == ctrl_curve])
      tdf <- prep_one(dt[curve == trt_curve])

      if (nrow(cdf) < 3 || nrow(tdf) < 3) return(NA_real_)

      tmax <- min(max(cdf$time), max(tdf$time))
      if (!is.finite(tmax) || tmax <= 0) return(NA_real_)

      tg <- seq(0, tmax, length.out = grid_n)

      # 用 step 函数近似
      Sc <- stats::approx(cdf$time, cdf$St, xout = tg, method = "constant", f = 0, rule = 2)$y
      St <- stats::approx(tdf$time, tdf$St, xout = tg, method = "constant", f = 0, rule = 2)$y

      Sc <- pmin(1 - eps, pmax(eps, Sc))
      St <- pmin(1 - eps, pmax(eps, St))

      Hc <- -log(Sc)
      Ht <- -log(St)

      ok <- is.finite(Hc) & is.finite(Ht) & (Hc > 0) & (Ht > 0)
      if (!any(ok)) return(NA_real_)

      beta_t <- log(Ht[ok] / Hc[ok])
      beta_t <- beta_t[is.finite(beta_t)]
      if (!length(beta_t)) return(NA_real_)

      stats::median(beta_t)
    }

    curves <- sort(unique(dt$curve))
    c1 <- curves[1]
    c2 <- curves[2]

    beta_hat <- estimate_logHR_from_KM(dt, ctrl_curve = c1, trt_curve = c2,
                                       grid_n = swap_grid_n, eps = swap_eps)

    swapped <- FALSE
    if (is.finite(beta_hat) && beta_hat > swap_tol) {
      # 交换曲线标签，使治疗组风险比 < 1（即治疗组曲线应在对照组下方）
      dt$curve <- ifelse(dt$curve == c1, c2,
                         ifelse(dt$curve == c2, c1, dt$curve))
      swapped <- TRUE
    }

    # 记录交换信息到 meta
    meta$curve_swapped <- swapped
    meta$digitized_logHR_hat <- beta_hat
  }

  #  可选生成 QA 图
  plot_path <- NA_character_
  if (isTRUE(make_plot)) {
    plot_path <- plot_out
    dir.create(dirname(plot_path), recursive = TRUE, showWarnings = FALSE)

    ok <- TRUE
    tryCatch({
      grDevices::png(plot_path, width = 1400, height = 900, res = 150)
    }, error = function(e){
      ok <<- FALSE
      warning(sprintf("Could not start png device: %s", e$message))
    })

    if (ok) {
      on.exit(grDevices::dev.off(), add = TRUE)

      graphics::plot(
        0, 0, type = "n",
        xlim = c(0, x_end), ylim = c(0, 1),
        xlab = "Time", ylab = "Survival",
        main = "Digitized KM points"
      )
      graphics::grid()

      curves <- sort(unique(dt$curve))
      cols <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e")

      for (i in seq_along(curves)) {
        dd <- dt[curve == curves[i]]
        graphics::points(dd$time, dd$St, pch = 16, cex = 0.7, col = cols[(i - 1) %% length(cols) + 1])
        graphics::lines(dd$time, dd$St, col = cols[(i - 1) %% length(cols) + 1], lwd = 1)
      }

      # 如果后续有 arm 信息可添加图例（此处未用，可保留）
      # if ("arm" %in% names(dt)) { ... }
    } else {
      plot_path <- NA_character_
    }
  }

  # 返回最终列表
  list(
    data = dt,
    plot = plot_path,
    meta = meta
  )
}
#
# ---- digitize_km_dense ----
#
# digitize_km_dense <- function(
#     ipd_true,
#     x_end,
#     n_points = 150,
#     sd_S = 0.04,
#     seed = NULL
# ) {
#   stopifnot(is.data.frame(ipd_true))
#   required_cols <- c("time", "status", "arm")
#   stopifnot(all(required_cols %in% names(ipd_true)))
#   stopifnot(is.numeric(x_end), length(x_end) == 1, x_end > 0)
#   stopifnot(is.numeric(n_points), n_points >= 2)
#
#   ipd_true$arm <- factor(ipd_true$arm, levels = c("Control", "Treatment"))
#   times <- seq(0, x_end, length.out = n_points)
#
#   if (!is.null(seed)) set.seed(seed)
#
#   km_noisy <- function(arm_label) {
#     dat <- ipd_true[ipd_true$arm == arm_label, c("time", "status")]
#     fit <- survival::survfit(survival::Surv(time, status) ~ 1, data = dat)
#
#     ss <- summary(fit, times = times, extend = TRUE)
#     St <- ss$surv
#
#     # 处理可能的 NA（理论上 extend=TRUE 不会产生 NA，但以防万一）
#     if (anyNA(St)) {
#       St <- zoo::na.locf(St, na.rm = FALSE)
#       St[is.na(St)] <- 1
#     }
#
#     St_noisy <- St + rnorm(length(St), 0, sd_S)
#     St_noisy[1] <- 1                     # 强制时间 0 处为 1
#     St_noisy <- pmax(pmin(St_noisy, 1), 0)   # 限制在 [0,1]
#     St_noisy <- cummin(St_noisy)              # 强制非增
#     St_noisy[1] <- 1                          # cummin 可能改变第一个点，重新设为 1
#
#     data.frame(time = times, St = St_noisy)
#   }
#
#   dt_control <- km_noisy("Control")
#   dt_treat   <- km_noisy("Treatment")
#
#   dt_control$curve <- 1L   # Control 对应 curve 1
#   dt_treat$curve   <- 2L   # Treatment 对应 curve 2
#
#   digi_dt <- rbind(dt_control, dt_treat)
#   digi_dt <- digi_dt[order(digi_dt$curve, digi_dt$time), ]
#
#   # 最终质量检查
#   stopifnot(!anyNA(digi_dt$St))
#   stopifnot(all(digi_dt$St >= 0 & digi_dt$St <= 1))
#
#   list(
#     data = digi_dt,
#     plot = NA_character_,
#     meta = list(
#       method = "dense_numeric",
#       x_end = x_end,
#       n_points = n_points,
#       sd_S = sd_S,
#       seed = seed
#     )
#   )
# }




#--------------------------------------------------------------#
# ---- Minimal usage example ----
#--------------------------------------------------------------#

library(here)

# 3) Digitize the KM image produced by simulate_trial
digi <- digitize_km_image(
  img_path = res3$km_png,
  x_end = res3$axis$x_end_plot,
  x_increment = res3$axis$x_increment,
  mode = "keep_censor_marks",
  num_curves = 2,
  make_plot = TRUE
)

# 新调用（使用数值模拟）
digi <- digitize_km_numeric(
  ipd_true = res$ipd,
  risk_table = res$risk_table,
  sd_S = 0.04,          # 可根据需要调整噪声大小
  seed = seed_run + 1L
)

digi$data
digi$plot
names(digi)
