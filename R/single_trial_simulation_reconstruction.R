# ============================================================
# Two-arm reconstruction from SurvDigitizeR digitization output
# Methods: IPDfromKM, kmdata::ipd, KMtoIPD
# ============================================================


reconstruct_ipd_twoarm <- function(
    digi,
    risk_table,
    method = c("IPDfromKM", "kmdata", "KMtoIPD"),
    curve_map = c(Control = 1, Treatment = 2),
    total_events = NULL,        # optional: list(Control=..., Treatment=...) or named numeric

    # NEW (2026-02): KMtoIPD scenario switch implemented here so benchmark can test both
    kmtoipd_scenario = c("no_marks", "with_marks"),

    # NEW (2026-02): kmdata stability knob (TRUE/number uses interpolate; fallback retry is automatic)
    kmdata_interpolate = FALSE
){
  method <- match.arg(method)
  if (method == "KMtoIPD") kmtoipd_scenario <- match.arg(kmtoipd_scenario)

  # --- accept digi as list (digitize_km_image output) or directly as data.frame/data.table
  digi_dt <- if (is.list(digi) && !is.null(digi$data)) digi$data else digi
  stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))

  # --- optional: censor marks for KMtoIPD (added by digitize_km_event_times)
  cens_mark <- NULL
  if (is.list(digi) && !is.null(digi$cens_mark)) cens_mark <- digi$cens_mark

  # --- risk table columns (fixed by your simulate_trial output)
  stopifnot(all(c("time", "nrisk_control", "nrisk_treat") %in% names(risk_table)))

  # --- helper: pull per-curve points and rename to x,y for IPDfromKM/kmdata
  get_curve_xy <- function(curve_id){
    d <- digi_dt[digi_dt$curve == curve_id, c("time", "St")]
    d <- d[order(d$time), ]
    d$St <- pmin(pmax(d$St, 0), 1)
    data.frame(x = d$time, y = d$St)
  }

  # --- helper: risk table for an arm
  get_arm_risk <- function(arm){
    rt <- risk_table[order(risk_table$time), ]
    nrisk <- if (arm == "Control") rt$nrisk_control else rt$nrisk_treat
    data.frame(trisk = rt$time, nrisk = nrisk)
  }

  # --- helper: total events (optional)
  get_tot_events <- function(arm){
    if (is.null(total_events)) return(NULL)
    if (is.list(total_events) && !is.null(total_events[[arm]])) return(as.integer(total_events[[arm]]))
    if (!is.null(names(total_events)) && arm %in% names(total_events)) return(as.integer(total_events[[arm]]))
    NULL
  }

  # --- helper: KMtoIPD needs "bottoms of drops": t (times), S (survival at bottoms)
  extract_drop_bottoms <- function(curve_id){
    d <- digi_dt[digi_dt$curve == curve_id, c("time", "St")]
    d <- d[order(d$time), ]
    d$St <- pmin(pmax(d$St, 0), 1)

    # For each unique time, take the MIN survival (bottom of any vertical drop)
    agg <- aggregate(St ~ time, data = d, FUN = min)
    agg <- agg[order(agg$time), ]

    # Identify times where survival decreases (drop bottoms are at the new, lower level)
    eps <- 1e-10
    idx <- which(diff(agg$St) < -eps) + 1L
    t_drop <- agg$time[idx]
    S_drop <- agg$St[idx]

    # Make sure survival is non-increasing (light stabilization)
    S_drop <- cummin(S_drop)

    list(t = t_drop, S = S_drop)
  }

  # --- reconstruct per arm
  reconstruct_one_arm <- function(arm){
    curve_id <- unname(curve_map[[arm]])
    xy <- get_curve_xy(curve_id)
    rt <- get_arm_risk(arm)

    # baseline n at (earliest) risk-table time
    n0 <- rt$nrisk[which.min(rt$trisk)]
    n0 <- as.integer(n0)
    tot_ev <- get_tot_events(arm)

    if (method == "IPDfromKM") {
      if (!requireNamespace("IPDfromKM", quietly = TRUE)) {
        stop("Package 'IPDfromKM' is required for method='IPDfromKM'.")
      }
      prep <- IPDfromKM::preprocess(
        dat = xy,
        trisk = rt$trisk,
        nrisk = rt$nrisk,
        totalpts = n0,
        maxy = 1
      )
      armID <- if (arm == "Control") 0 else 1
      ipd_res <- IPDfromKM::getIPD(prep = prep, armID = armID, tot.events = tot_ev)
      ipd_tbl <- ipd_res$IPD

      out <- data.frame(
        time   = ipd_tbl$time,
        status = ipd_tbl$status,
        arm    = arm,
        stringsAsFactors = FALSE
      )
      return(out[order(out$time), ])

    } else if (method == "kmdata") {
      if (!requireNamespace("kmdata", quietly = TRUE)) {
        stop("Package 'kmdata' is required for method='kmdata'.")
      }

      # NEW (2026-02): pre-clean inputs to avoid kmdata::ipd NA/0 edge cases
      xy2 <- xy
      xy2 <- xy2[is.finite(xy2$x) & is.finite(xy2$y) & xy2$x >= 0, , drop = FALSE]
      if (!nrow(xy2)) stop("kmdata: empty digitized curve after filtering")
      # per-time keep MIN prob (bottom) to avoid vertical-segment duplicates
      xy2 <- aggregate(y ~ x, data = xy2, FUN = min)
      xy2 <- xy2[order(xy2$x), , drop = FALSE]
      if (!any(xy2$x == 0)) xy2 <- rbind(data.frame(x = 0, y = 1), xy2)
      xy2 <- xy2[order(xy2$x), , drop = FALSE]
      xy2$y <- pmin(pmax(xy2$y, 0), 1)
      xy2$y <- cummin(xy2$y)
      # avoid exact 0/1 leading to division issues inside kmdata
      eps <- 1e-6
      xy2$y <- pmin(1 - eps, pmax(eps, xy2$y))

      rt2 <- rt
      rt2 <- rt2[is.finite(rt2$trisk) & is.finite(rt2$nrisk), , drop = FALSE]
      rt2 <- rt2[rt2$nrisk > 0, , drop = FALSE]  # critical for kmdata stability
      rt2 <- rt2[order(rt2$trisk), , drop = FALSE]
      # de-duplicate trisk
      rt2 <- rt2[!duplicated(rt2$trisk), , drop = FALSE]
      if (!nrow(rt2)) stop("kmdata: empty risk table after filtering nrisk>0")

      # detect whether kmdata::ipd supports 'interpolate'
      has_interp <- FALSE
      ff <- tryCatch(formals(kmdata::ipd), error = function(e) NULL)
      if (!is.null(ff)) has_interp <- ("interpolate" %in% names(ff))

      interp_val <- NULL
      if (isTRUE(kmdata_interpolate)) interp_val <- 500
      if (is.numeric(kmdata_interpolate) && length(kmdata_interpolate) == 1 && is.finite(kmdata_interpolate)) {
        interp_val <- as.integer(kmdata_interpolate)
      }

      call_ipd <- function(interp = NULL) {
        args <- list(
          time = xy2$x,
          prob = xy2$y,
          t.atrisk = rt2$trisk,
          n.atrisk = rt2$nrisk,
          arm = arm
        )
        if (!is.null(interp) && has_interp) args$interpolate <- interp
        do.call(kmdata::ipd, args)
      }

      ipd <- tryCatch(call_ipd(interp = interp_val), error = function(e) {
        # NEW: fallback retry with interpolate=500 if available
        if (has_interp) {
          return(call_ipd(interp = 500L))
        }
        stop(e)
      })

      out <- data.frame(time = ipd$time, status = ipd$event, arm = arm, stringsAsFactors = FALSE)
      return(out[order(out$time), ])

    } else { # KMtoIPD
      if (!requireNamespace("KMtoIPD", quietly = TRUE)) {
        stop("Package 'KMtoIPD' is required for method='KMtoIPD'.")
      }

      # tail time (study end) -> added to cens.t (per our benchmark convention)
      t_end <- max(risk_table$time[is.finite(risk_table$time)], na.rm = TRUE)
      if (!is.finite(t_end) || t_end <= 0) {
        # fallback to last observed time
        t_end <- max(xy$x[is.finite(xy$x)], na.rm = TRUE)
      }

      drops <- extract_drop_bottoms(curve_id)

      # Scenario A/B: censor marks
      cens_vec <- numeric(0)
      if (kmtoipd_scenario == "with_marks" && !is.null(cens_mark) && arm %in% names(cens_mark)) {
        cens_vec <- cens_mark[[arm]]
        cens_vec <- cens_vec[is.finite(cens_vec) & cens_vec >= 0]
      }

      cens_use <- sort(unique(c(cens_vec, t_end)))

      # Avoid exact ties between censor marks and drop times: ensure events happen before censoring
      if (length(drops$t) > 0 && length(cens_use) > 0) {
        tie <- cens_use %in% drops$t
        cens_use[tie] <- cens_use[tie] + 1e-8
      }

      # Edge case: no drops (no events) -> all censored at t_end
      if (length(drops$t) == 0L) {
        out <- data.frame(time = rep(t_end, n0), status = 0L, arm = arm, stringsAsFactors = FALSE)
        return(out)
      }

      ipd <- KMtoIPD::getIPD(
        n = n0,
        t = drops$t,
        S = drops$S,
        cens.t = cens_use
      )

      out <- data.frame(time = ipd$t, status = ipd$event, arm = arm, stringsAsFactors = FALSE)
      return(out[order(out$time), ])
    }
  }

  ipd_c <- reconstruct_one_arm("Control")
  ipd_t <- reconstruct_one_arm("Treatment")

  ipd_all <- rbind(ipd_c, ipd_t)
  ipd_all$arm <- factor(ipd_all$arm, levels = c("Control", "Treatment"))

  list(
    method = method,
    kmtoipd_scenario = if (method == "KMtoIPD") kmtoipd_scenario else NA_character_,
    curve_map = curve_map,
    ipd = ipd_all,
    ipd_control = ipd_c,
    ipd_treat = ipd_t
  )
}


sanity_plot_overlay <- function(truth_ipd,
                                digi_dt,
                                rec_list,
                                rec_names = names(rec_list),
                                curve_map = c(Control = 1, Treatment = 2),
                                out_file = NULL,
                                width = 2000, height = 700, res = 150) {

  stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))
  stopifnot(all(c("time", "status", "arm") %in% names(truth_ipd)))

  if (is.null(rec_names)) rec_names <- paste0("rec", seq_along(rec_list))
  if (length(rec_names) != length(rec_list)) rec_names <- paste0("rec", seq_along(rec_list))

  # --- open device if requested ---
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(out_file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, length(rec_list)), mar = c(4, 4, 3, 1))

  # ---------- style: colors by arm ----------
  col_arm <- c(Control = "dodgerblue3", Treatment = "firebrick3")

  # ---------- helper: survfit -> step data including (0,1) ----------
  survfit_to_step <- function(fit, arm_level, x_end = NULL) {
    s <- summary(fit)
    strata <- s$strata
    tt <- s$time
    ss <- s$surv

    # strata label looks like "arm=Control"
    want <- paste0("arm=", arm_level)
    idx <- which(strata == want)
    if (length(idx) == 0) stop("Cannot find strata ", want, " in survfit.")

    t <- tt[idx]
    y <- ss[idx]

    # Ensure starts at (0,1)
    if (length(t) == 0 || t[1] > 0) {
      t <- c(0, t)
      y <- c(1, y)
    } else {
      # if first time is 0 but surv not 1, force it
      t[1] <- 0
      y[1] <- 1
    }

    # Optionally extend to x_end
    if (!is.null(x_end) && is.finite(x_end) && tail(t, 1) < x_end) {
      t <- c(t, x_end)
      y <- c(y, tail(y, 1))
    }

    list(time = t, surv = y)
  }

  # ---------- truth KM ----------
  truth_ipd$arm <- factor(truth_ipd$arm, levels = c("Control", "Treatment"))
  fit_truth <- survival::survfit(survival::Surv(time, status) ~ arm, data = truth_ipd)

  # max x across everything (for consistent panels)
  get_ipd_time_max <- function(x) {
    ipd <- if (is.list(x) && !is.null(x$ipd)) x$ipd else x
    max(ipd$time, na.rm = TRUE)
  }
  x_max <- max(
    truth_ipd$time,
    digi_dt$time,
    vapply(rec_list, get_ipd_time_max, numeric(1)),
    na.rm = TRUE
  )

  # ---------- digitized points by curve ----------
  add_digitized <- function(curve_id, arm_label) {
    dd <- digi_dt[digi_dt$curve == curve_id, , drop = FALSE]
    dd <- dd[order(dd$time), , drop = FALSE]
    # points style: Control filled, Treatment open
    if (arm_label == "Control") {
      points(dd$time, dd$St, pch = 16, cex = 0.30, col = col_arm[arm_label])
    } else {
      points(dd$time, dd$St, pch = 1,  cex = 0.35, col = col_arm[arm_label])
    }
  }

  # ---------- per method panels ----------
  for (j in seq_along(rec_list)) {
    rec <- rec_list[[j]]
    rec_ipd <- if (is.list(rec) && !is.null(rec$ipd)) rec$ipd else rec
    stopifnot(all(c("time", "status", "arm") %in% names(rec_ipd)))

    rec_ipd$arm <- factor(rec_ipd$arm, levels = c("Control", "Treatment"))
    fit_rec <- survival::survfit(survival::Surv(time, status) ~ arm, data = rec_ipd)

    plot(NA, xlim = c(0, x_max), ylim = c(0, 1),
         xlab = "Time", ylab = "Survival",
         main = rec_names[j])

    # Truth (solid)
    for (arm in c("Control", "Treatment")) {
      st <- survfit_to_step(fit_truth, arm_level = arm, x_end = x_max)
      lines(st$time, st$surv, type = "s", lty = 1, lwd = 2, col = col_arm[arm])
    }

    # Digitized points
    add_digitized(curve_map[["Control"]], "Control")
    add_digitized(curve_map[["Treatment"]], "Treatment")

    # Reconstructed (dashed)
    for (arm in c("Control", "Treatment")) {
      st <- survfit_to_step(fit_rec, arm_level = arm, x_end = x_max)
      lines(st$time, st$surv, type = "s", lty = 2, lwd = 2, col = col_arm[arm])
    }

    legend("topright",
           legend = c("Truth KM", "Digitized points", "Reconstructed KM",
                      "Control (blue)", "Treatment (red)"),
           lty = c(1, NA, 2, NA, NA),
           pch = c(NA, 16, NA, NA, NA),
           col = c("black", "black", "black", col_arm["Control"], col_arm["Treatment"]),
           bty = "n", cex = 0.8)
  }

  invisible(out_file)
}


# --------------------------
# Minimal test example
if (FALSE) {

  # --------------------------
  out_png <- file.path("C:/Amber/JHU/km2ipd/Benchmark/test",
                       "sanity_overlay_W_shape0p6_HR085_2026.png")


  rec1 <- reconstruct_ipd_twoarm(digi, res3$risk_table, method = "IPDfromKM",
                                 curve_map = c(Control = 1, Treatment = 2))
  rec2 <- reconstruct_ipd_twoarm(digi, res3$risk_table, method = "kmdata",
                                 total_events = c(Control = res3$totals$total_events_ctrl,
                                                  Treatment = res3$totals$total_events_trt),
                                 curve_map = c(Control = 1, Treatment = 2))
  rec3 <- reconstruct_ipd_twoarm(digi, res3$risk_table, method = "KMtoIPD",
                                 curve_map = c(Control = 1, Treatment = 2))



  sanity_plot_overlay(
    truth_ipd = res3$ipd,
    digi_dt = digi$data,
    rec_list = list(rec1 = rec1, rec2 = rec2, rec3 = rec3),
    rec_names = c("IPDfromKM", "kmdata::ipd", "KMtoIPD"),
    curve_map = c(Control = 1, Treatment = 2),  # 如果你发现反了就 swap
    out_file = out_png
  )

}
