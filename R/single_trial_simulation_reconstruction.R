# ============================================================
# Two-arm reconstruction from SurvDigitizeR digitization output
# Methods: IPDfromKM, kmdata::ipd, KMtoIPD
# ============================================================


reconstruct_ipd_twoarm <- function(
    digi,
    risk_table,
    method = c("IPDfromKM", "kmdata", "KMtoIPD_useMarks", "KMtoIPD_noMarks"),
    curve_map = c(Control = 1, Treatment = 2),
    total_events = NULL         # optional: list(Control=..., Treatment=...) or named numeric
){
  method <- match.arg(method)

  # --- accept digi as list (digitize_km_* output) or directly as data.frame/data.table
  digi_dt <- if (is.list(digi) && !is.null(digi$data)) digi$data else digi
  stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))

  # --- NEW: optional censor marks from digitization (KMtoIPD can incorporate these)
  # digitize_km_event_times() now outputs $censor_marks and $censor_marks_by_curve.
  digi_cens_marks <- NULL
  if (is.list(digi) && !is.null(digi$censor_marks)) {
    digi_cens_marks <- digi$censor_marks
  }
  digi_cens_marks_by_curve <- NULL
  if (is.list(digi) && !is.null(digi$censor_marks_by_curve)) {
    digi_cens_marks_by_curve <- digi$censor_marks_by_curve
  }

  # --- risk table columns (fixed by your simulate_trial output)
  stopifnot(all(c("time", "nrisk_control", "nrisk_treat") %in% names(risk_table)))

  # --- helper: pull per-curve points and rename to x,y for IPDfromKM
  get_curve_xy <- function(curve_id){
    d <- digi_dt[digi_dt$curve == curve_id, c("time", "St")]
    d <- d[order(d$time), ]
    # minimal bounding for safety (not "cleaning" shape)
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

  # --- NEW helper: fetch censor marks for a given arm/curve
  get_arm_censor_marks <- function(arm, curve_id){
    # Prefer per-arm marks if present
    if (is.list(digi_cens_marks) && !is.null(digi_cens_marks[[arm]])) {
      return(as.numeric(digi_cens_marks[[arm]]))
    }
    # Fallback: per-curve marks if present
    if (is.list(digi_cens_marks_by_curve)) {
      key <- as.character(curve_id)
      if (!is.null(digi_cens_marks_by_curve[[key]])) {
        return(as.numeric(digi_cens_marks_by_curve[[key]]))
      }
    }
    numeric(0)
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

    # Make sure survival is non-increasing (very light stabilization)
    S_drop <- cummin(S_drop)

    list(t = t_drop, S = S_drop)
  }

  # --- reconstruct per arm
  reconstruct_one_arm <- function(arm){
    curve_id <- unname(curve_map[[arm]])
    xy <- get_curve_xy(curve_id)
    rt <- get_arm_risk(arm)

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
      ipd_tbl <- ipd_res$IPD  # <-- the 3-column IPD table lives here

      # ipd_tbl columns are typically: time, status, armID
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
      ipd <- kmdata::ipd(
        time = xy$x,
        prob = xy$y,
        t.atrisk = rt$trisk,
        n.atrisk = rt$nrisk,
        arm = arm
      )
      out <- data.frame(time = ipd$time, status = ipd$event, arm = arm, stringsAsFactors = FALSE)
      return(out[order(out$time), ])

    } else { # KMtoIPD (two scenarios)
      if (!requireNamespace("KMtoIPD", quietly = TRUE)) {
        stop("Package 'KMtoIPD' is required for method starting with 'KMtoIPD_...'.")
      }

      # --- Extract the drop-bottom coordinates required by Rogula 2022.
      # NOTE: digi_dt is dense noisy digitized points; we extract the bottoms of drops here.
      drops <- extract_drop_bottoms(curve_id)

      # --- IMPORTANT tail handling (per Rogula 2022):
      # If censor marks are not provided, the end of the KM curve can be included in cens.t.
      # Here we always add an administrative end time (from risk table) to cens.t to
      # stabilize the curve tail (horizontal end segments).
      t_end <- max(risk_table$time, na.rm = TRUE)
      if (!is.finite(t_end) || t_end <= 0) {
        stop("risk_table$time must contain a positive finite follow-up end time for KMtoIPD tail handling.")
      }

      # --- Scenario A/B as separate method names:
      #   KMtoIPD_useMarks: incorporate digitized censor marks + add t_end
      #   KMtoIPD_noMarks : ignore digitized marks; use only t_end
      if (method == "KMtoIPD_useMarks") {
        # NEW: incorporate marked censoring times from digitization
        cens_in <- get_arm_censor_marks(arm = arm, curve_id = curve_id)
        cens_in <- cens_in[is.finite(cens_in)]
        cens_in <- pmin(pmax(cens_in, 0), t_end)
        cens_use <- sort(unique(c(cens_in, t_end)))
      } else {
        cens_use <- t_end
      }

      # If digitization is extremely noisy or the true data has ~0 events,
      # we may fail to detect any drops. In that case, return a sensible
      # all-censored IPD (consistent with a flat KM curve).
      if (length(drops$t) == 0L) {
        cens_vec <- as.numeric(cens_use)
        cens_vec <- cens_vec[is.finite(cens_vec)]
        if (length(cens_vec) > n0) cens_vec <- cens_vec[seq_len(n0)]
        if (length(cens_vec) < n0) cens_vec <- c(cens_vec, rep(t_end, n0 - length(cens_vec)))
        ipd <- data.frame(t = cens_vec, event = rep(0L, length(cens_vec)))
      } else {
        ipd <- KMtoIPD::getIPD(
          n = n0,
          t = drops$t,
          S = drops$S,
          cens.t = cens_use
        )
      }
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
                                width = 2500, height = 700, res = 150,
                                digi_max_points = 600,        # FIX: 抽样点数，避免太密看不清
                                digi_alpha = 0.35             # FIX: 点透明度
) {

  stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))
  stopifnot(all(c("time", "status", "arm") %in% names(truth_ipd)))

  if (is.null(rec_names)) rec_names <- paste0("rec", seq_along(rec_list))
  if (length(rec_names) != length(rec_list)) rec_names <- paste0("rec", seq_along(rec_list))

  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(out_file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, length(rec_list)), mar = c(4, 4, 3, 1))

  col_arm <- c(Control = "dodgerblue3", Treatment = "firebrick3")

  survfit_to_step <- function(fit, arm_level, x_end = NULL) {
    s <- summary(fit)
    want <- paste0("arm=", arm_level)
    idx <- which(s$strata == want)
    if (length(idx) == 0) stop("Cannot find strata ", want, " in survfit.")
    t <- s$time[idx]
    y <- s$surv[idx]

    if (length(t) == 0 || t[1] > 0) {
      t <- c(0, t); y <- c(1, y)
    } else {
      t[1] <- 0; y[1] <- 1
    }
    if (!is.null(x_end) && is.finite(x_end) && tail(t, 1) < x_end) {
      t <- c(t, x_end)
      y <- c(y, tail(y, 1))
    }
    list(time = t, surv = y)
  }

  truth_ipd$arm <- factor(truth_ipd$arm, levels = c("Control", "Treatment"))
  fit_truth <- survival::survfit(survival::Surv(time, status) ~ arm, data = truth_ipd)

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

  # FIX: digitized 点抽样 + 透明度 + 用 arm 颜色
  add_digitized <- function(curve_id, arm_label) {
    dd <- digi_dt[digi_dt$curve == curve_id, , drop = FALSE]
    dd <- dd[order(dd$time), , drop = FALSE]

    if (nrow(dd) > digi_max_points) {
      keep <- unique(round(seq(1, nrow(dd), length.out = digi_max_points)))
      dd <- dd[keep, , drop = FALSE]
    }

    col_dot <- grDevices::adjustcolor(col_arm[arm_label], alpha.f = digi_alpha)
    if (arm_label == "Control") {
      points(dd$time, dd$St, pch = 16, cex = 0.55, col = "pink")  # FIX: 更大
    } else {
      points(dd$time, dd$St, pch = 1,  cex = 0.65, col = "grey")  # FIX: 更大
    }
  }

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

    # Reconstructed (dashed)  --- FIX: 先画 reconstruction，再画 points
    for (arm in c("Control", "Treatment")) {
      st <- survfit_to_step(fit_rec, arm_level = arm, x_end = x_max)
      lines(st$time, st$surv, type = "s", lty = 2, lwd = 2, col = col_arm[arm])
    }

    # Digitized points (draw last so they are visible)
    add_digitized(curve_map[["Control"]], "Control")
    add_digitized(curve_map[["Treatment"]], "Treatment")

    legend("topright",
           legend = c("Truth (solid)", "Reconstructed (dashed)", "Digitized points"),
           lty = c(1, 2, NA),
           pch = c(NA, NA, 16),
           col = c("black", "black", "black"),
           bty = "n", cex = 0.85)
  }

  invisible(out_file)
}



# reconstruct_ipd_twoarm <- function(
#     digi,
#     risk_table,
#     method = c("IPDfromKM", "kmdata", "KMtoIPD"),
#     curve_map = c(Control = 1, Treatment = 2),
#     total_events = NULL         # optional: list(Control=..., Treatment=...) or named numeric
# ){
#   method <- match.arg(method)
#
#   # --- accept digi as list (digitize_km_image output) or directly as data.frame/data.table
#   digi_dt <- if (is.list(digi) && !is.null(digi$data)) digi$data else digi
#   stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))
#
#   # --- risk table columns (fixed by your simulate_trial output)
#   stopifnot(all(c("time", "nrisk_control", "nrisk_treat") %in% names(risk_table)))
#
#   # --- helper: pull per-curve points and rename to x,y for IPDfromKM
#   get_curve_xy <- function(curve_id){
#     d <- digi_dt[digi_dt$curve == curve_id, c("time", "St")]
#     d <- d[order(d$time), ]
#     # minimal bounding for safety (not "cleaning" shape)
#     d$St <- pmin(pmax(d$St, 0), 1)
#     data.frame(x = d$time, y = d$St)
#   }
#
#   # --- helper: risk table for an arm
#   get_arm_risk <- function(arm){
#     rt <- risk_table[order(risk_table$time), ]
#     nrisk <- if (arm == "Control") rt$nrisk_control else rt$nrisk_treat
#     data.frame(trisk = rt$time, nrisk = nrisk)
#   }
#
#   # --- helper: total events (optional)
#   get_tot_events <- function(arm){
#     if (is.null(total_events)) return(NULL)
#     if (is.list(total_events) && !is.null(total_events[[arm]])) return(as.integer(total_events[[arm]]))
#     if (!is.null(names(total_events)) && arm %in% names(total_events)) return(as.integer(total_events[[arm]]))
#     NULL
#   }
#
#   # --- helper: KMtoIPD needs "bottoms of drops": t (times), S (survival at bottoms)
#   extract_drop_bottoms <- function(curve_id){
#     d <- digi_dt[digi_dt$curve == curve_id, c("time", "St")]
#     d <- d[order(d$time), ]
#     d$St <- pmin(pmax(d$St, 0), 1)
#
#     # For each unique time, take the MIN survival (bottom of any vertical drop)
#     agg <- aggregate(St ~ time, data = d, FUN = min)
#     agg <- agg[order(agg$time), ]
#
#     # Identify times where survival decreases (drop bottoms are at the new, lower level)
#     eps <- 1e-10
#     idx <- which(diff(agg$St) < -eps) + 1L
#     t_drop <- agg$time[idx]
#     S_drop <- agg$St[idx]
#
#     # Make sure survival is non-increasing (very light stabilization)
#     S_drop <- cummin(S_drop)
#
#     list(t = t_drop, S = S_drop)
#   }
#
#   # --- reconstruct per arm
#   reconstruct_one_arm <- function(arm){
#     curve_id <- unname(curve_map[[arm]])
#     xy <- get_curve_xy(curve_id)
#     rt <- get_arm_risk(arm)
#
#     n0 <- rt$nrisk[which.min(rt$trisk)]
#     n0 <- as.integer(n0)
#     tot_ev <- get_tot_events(arm)
#
#     if (method == "IPDfromKM") {
#       if (!requireNamespace("IPDfromKM", quietly = TRUE)) {
#         stop("Package 'IPDfromKM' is required for method='IPDfromKM'.")
#       }
#       prep <- IPDfromKM::preprocess(
#         dat = xy,
#         trisk = rt$trisk,
#         nrisk = rt$nrisk,
#         totalpts = n0,
#         maxy = 1
#       )
#       armID <- if (arm == "Control") 0 else 1
#       ipd_res <- IPDfromKM::getIPD(prep = prep, armID = armID, tot.events = tot_ev)
#       ipd_tbl <- ipd_res$IPD  # <-- the 3-column IPD table lives here
#
#       # ipd_tbl columns are typically: time, status, armID
#       out <- data.frame(
#         time   = ipd_tbl$time,
#         status = ipd_tbl$status,
#         arm    = arm,
#         stringsAsFactors = FALSE
#       )
#
#       return(out[order(out$time), ])
#
#     } else if (method == "kmdata") {
#       if (!requireNamespace("kmdata", quietly = TRUE)) {
#         stop("Package 'kmdata' is required for method='kmdata'.")
#       }
#       ipd <- kmdata::ipd(
#         time = xy$x,
#         prob = xy$y,
#         t.atrisk = rt$trisk,
#         n.atrisk = rt$nrisk,
#         arm = arm
#       )
#       out <- data.frame(time = ipd$time, status = ipd$event, arm = arm, stringsAsFactors = FALSE)
#       return(out[order(out$time), ])
#
#     } else { # KMtoIPD
#       if (!requireNamespace("KMtoIPD", quietly = TRUE)) {
#         stop("Package 'KMtoIPD' is required for method='KMtoIPD'.")
#       }
#       drops <- extract_drop_bottoms(curve_id)
#       ipd <- KMtoIPD::getIPD(
#         n = n0,
#         t = drops$t,
#         S = drops$S,
#         cens.t = NA
#       )
#       out <- data.frame(time = ipd$t, status = ipd$event, arm = arm, stringsAsFactors = FALSE)
#       return(out[order(out$time), ])
#     }
#   }
#
#   ipd_c <- reconstruct_one_arm("Control")
#   ipd_t <- reconstruct_one_arm("Treatment")
#
#   ipd_all <- rbind(ipd_c, ipd_t)
#   ipd_all$arm <- factor(ipd_all$arm, levels = c("Control", "Treatment"))
#
#   list(
#     method = method,
#     curve_map = curve_map,
#     ipd = ipd_all,
#     ipd_control = ipd_c,
#     ipd_treat = ipd_t
#   )
# }


# sanity_plot_overlay <- function(truth_ipd,
#                                 digi_dt,
#                                 rec_list,
#                                 rec_names = names(rec_list),
#                                 curve_map = c(Control = 1, Treatment = 2),
#                                 out_file = NULL,
#                                 width = 2000, height = 700, res = 150) {
#
#   stopifnot(all(c("time", "St", "curve") %in% names(digi_dt)))
#   stopifnot(all(c("time", "status", "arm") %in% names(truth_ipd)))
#
#   if (is.null(rec_names)) rec_names <- paste0("rec", seq_along(rec_list))
#   if (length(rec_names) != length(rec_list)) rec_names <- paste0("rec", seq_along(rec_list))
#
#   # --- open device if requested ---
#   if (!is.null(out_file)) {
#     dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
#     grDevices::png(out_file, width = width, height = height, res = res)
#     on.exit(grDevices::dev.off(), add = TRUE)
#   }
#
#   oldpar <- par(no.readonly = TRUE)
#   on.exit(par(oldpar), add = TRUE)
#   par(mfrow = c(1, length(rec_list)), mar = c(4, 4, 3, 1))
#
#   # ---------- style: colors by arm ----------
#   col_arm <- c(Control = "dodgerblue3", Treatment = "firebrick3")
#
#   # ---------- helper: survfit -> step data including (0,1) ----------
#   survfit_to_step <- function(fit, arm_level, x_end = NULL) {
#     s <- summary(fit)
#     strata <- s$strata
#     tt <- s$time
#     ss <- s$surv
#
#     # strata label looks like "arm=Control"
#     want <- paste0("arm=", arm_level)
#     idx <- which(strata == want)
#     if (length(idx) == 0) stop("Cannot find strata ", want, " in survfit.")
#
#     t <- tt[idx]
#     y <- ss[idx]
#
#     # Ensure starts at (0,1)
#     if (length(t) == 0 || t[1] > 0) {
#       t <- c(0, t)
#       y <- c(1, y)
#     } else {
#       # if first time is 0 but surv not 1, force it
#       t[1] <- 0
#       y[1] <- 1
#     }
#
#     # Optionally extend to x_end
#     if (!is.null(x_end) && is.finite(x_end) && tail(t, 1) < x_end) {
#       t <- c(t, x_end)
#       y <- c(y, tail(y, 1))
#     }
#
#     list(time = t, surv = y)
#   }
#
#   # ---------- truth KM ----------
#   truth_ipd$arm <- factor(truth_ipd$arm, levels = c("Control", "Treatment"))
#   fit_truth <- survival::survfit(survival::Surv(time, status) ~ arm, data = truth_ipd)
#
#   # max x across everything (for consistent panels)
#   get_ipd_time_max <- function(x) {
#     ipd <- if (is.list(x) && !is.null(x$ipd)) x$ipd else x
#     max(ipd$time, na.rm = TRUE)
#   }
#   x_max <- max(
#     truth_ipd$time,
#     digi_dt$time,
#     vapply(rec_list, get_ipd_time_max, numeric(1)),
#     na.rm = TRUE
#   )
#
#   # ---------- digitized points by curve ----------
#   add_digitized <- function(curve_id, arm_label) {
#     dd <- digi_dt[digi_dt$curve == curve_id, , drop = FALSE]
#     dd <- dd[order(dd$time), , drop = FALSE]
#     # points style: Control filled, Treatment open
#     if (arm_label == "Control") {
#       points(dd$time, dd$St, pch = 16, cex = 0.30, col = col_arm[arm_label])
#     } else {
#       points(dd$time, dd$St, pch = 1,  cex = 0.35, col = col_arm[arm_label])
#     }
#   }
#
#   # ---------- per method panels ----------
#   for (j in seq_along(rec_list)) {
#     rec <- rec_list[[j]]
#     rec_ipd <- if (is.list(rec) && !is.null(rec$ipd)) rec$ipd else rec
#     stopifnot(all(c("time", "status", "arm") %in% names(rec_ipd)))
#
#     rec_ipd$arm <- factor(rec_ipd$arm, levels = c("Control", "Treatment"))
#     fit_rec <- survival::survfit(survival::Surv(time, status) ~ arm, data = rec_ipd)
#
#     plot(NA, xlim = c(0, x_max), ylim = c(0, 1),
#          xlab = "Time", ylab = "Survival",
#          main = rec_names[j])
#
#     # Truth (solid)
#     for (arm in c("Control", "Treatment")) {
#       st <- survfit_to_step(fit_truth, arm_level = arm, x_end = x_max)
#       lines(st$time, st$surv, type = "s", lty = 1, lwd = 2, col = col_arm[arm])
#     }
#
#     # Digitized points
#     add_digitized(curve_map[["Control"]], "Control")
#     add_digitized(curve_map[["Treatment"]], "Treatment")
#
#     # Reconstructed (dashed)
#     for (arm in c("Control", "Treatment")) {
#       st <- survfit_to_step(fit_rec, arm_level = arm, x_end = x_max)
#       lines(st$time, st$surv, type = "s", lty = 2, lwd = 2, col = col_arm[arm])
#     }
#
#     legend("topright",
#            legend = c("Truth KM", "Digitized points", "Reconstructed KM",
#                       "Control (blue)", "Treatment (red)"),
#            lty = c(1, NA, 2, NA, NA),
#            pch = c(NA, 16, NA, NA, NA),
#            col = c("black", "black", "black", col_arm["Control"], col_arm["Treatment"]),
#            bty = "n", cex = 0.8)
#   }
#
#   invisible(out_file)
# }


# --------------------------
# Minimal test example
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
  digi = digi,
  rec_list = list(rec1 = rec1, rec2 = rec2, rec3 = rec3),
  rec_names = c("IPDfromKM", "kmdata::ipd", "KMtoIPD"),
  curve_map = c(Control = 1, Treatment = 2),  # 如果你发现反了就 swap
  out_file = out_png
)
