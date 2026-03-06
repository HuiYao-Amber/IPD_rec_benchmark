# ============================================================
# single_trial_simulation_benchmark.R
# Benchmark runner: DGM -> digitize -> reconstruct -> metrics
# ============================================================
run_single_trial_benchmark <- function(
    design,
    n_rep = 10,
    scenario_lib = NULL,
    methods = c("IPDfromKM", "kmdata", "KMtoIPD"),
    curve_map = c(Control = 1, Treatment = 2),
    base_seed = 202500,
    out_dir = tempdir(),
    keep_files = FALSE,          # keep intermediate PNGs for debugging
    keep_failed_files = TRUE,    # if FALSE, even failed runs will delete files
    verbose = TRUE,

    # ---- digitization knobs (passed into digitize_km_image) ----
    digitize_mode = "keep_censor_marks",
    digitize_bg_lightness = 0.4,
    digitize_nr_neighbors = 20,
    digitize_sd_S = 0.04,
    digitize_make_plot = FALSE,  # QA points plot (extra PNG). Usually FALSE for large benchmark.

    # ---- reconstruction knobs ----
    kmdata_interpolate = FALSE,

    # ---- metrics knobs ----
    metrics_do_rmst = TRUE,
    metrics_grid_n = 300
){

  # Basic checks
  stopifnot(is.data.frame(design))
  req_cols <- c("scenario_id", "censoring", "target_censoring", "n_control", "n_treatment")
  miss <- setdiff(req_cols, names(design))
  if (length(miss)) stop("design is missing columns: ", paste(miss, collapse = ", "))

  # scenario library (build once)
  if (is.null(scenario_lib)) scenario_lib <- make_scenario_library(mean_time = 10)

  # output folders
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_root <- if (keep_files) out_dir else file.path(out_dir, "_tmp_km2ipd")
  if (!dir.exists(tmp_root)) dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)

  # helper: safe bind rows without extra deps
  bind_rows_safe <- function(x){
    x <- x[!vapply(x, is.null, logical(1))]
    if (!length(x)) return(data.frame())
    # make sure each is a data.frame
    x <- lapply(x, function(z) as.data.frame(z, stringsAsFactors = FALSE))
    # fill missing columns
    all_cols <- unique(unlist(lapply(x, names)))
    x <- lapply(x, function(z){
      for (cc in setdiff(all_cols, names(z))) z[[cc]] <- NA
      z <- z[, all_cols, drop = FALSE]
      z
    })
    do.call(rbind, x)
  }

  # results collectors
  rows <- vector("list", length = nrow(design) * n_rep * length(methods))
  k_row <- 0L

  if (verbose) {
    message("Benchmark starts: ", nrow(design), " settings x ", n_rep, " reps x ",
            length(methods), " methods = ", nrow(design) * n_rep * length(methods), " rows")
  }

  # Main loops
  for (i in seq_len(nrow(design))) {

    scenario_id      <- as.character(design$scenario_id[i])
    censoring        <- as.character(design$censoring[i])
    target_censoring <- as.numeric(design$target_censoring[i])
    n_control        <- as.integer(design$n_control[i])
    n_treatment      <- as.integer(design$n_treatment[i])

    for (r in seq_len(n_rep)) {

      seed_run <- base_seed + i * 10000L + r
      run_id <- paste0(scenario_id, "_", censoring, "_c", sprintf("%03d", round(100 * target_censoring)),
                       "_n", n_control, "v", n_treatment, "_rep", sprintf("%03d", r))

      run_dir <- if (keep_files) file.path(tmp_root, scenario_id, censoring, paste0("rep_", r)) else tmp_root
      if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

      prefix <- run_id

      # track files for cleanup
      km_png <- NA_character_
      digi_png <- NA_character_


      # 1) DGM (simulate + KM PNG)
      dgm_err <- NA_character_
      res <- tryCatch({
        simulate_trial(
          scenario_id = scenario_id,
          scenario_lib = scenario_lib,
          n_control = n_control,
          n_treatment = n_treatment,
          censoring = censoring,
          target_censoring = target_censoring,
          seed = seed_run,
          out_dir = run_dir,
          prefix = prefix,
          plot = "km",
          add_censor_marks = FALSE
        )
      }, error = function(e){
        dgm_err <<- e$message
        NULL
      })

      if (!is.null(res)) km_png <- res$km_png

      dgm_err <- NA_character_
      res <- tryCatch({
        simulate_trial(
          scenario_id = scenario_id,
          scenario_lib = scenario_lib,
          n_control = n_control,
          n_treatment = n_treatment,
          censoring = censoring,
          target_censoring = target_censoring,
          seed = seed_run,
          out_dir = run_dir,
          prefix = prefix,
          plot = "km",                 # 仍然生成 PNG，可用于事后检查
          add_censor_marks = FALSE
        )
      }, error = function(e) { ... })

      if (!is.null(res)) km_png <- res$km_png

      # 2) Digitization (使用数值模拟，不再检查 PNG 文件)
      digi_err <- NA_character_
      digi <- NULL
      if (!is.null(res)) {   # 只要模拟成功，就进行数字化
        digi <- tryCatch({
          digitize_km_event_times(
            ipd_true = res$ipd,
            sd_S = digitize_sd_S,
            seed = seed_run + 1L
          )
        }, error = function(e) {
          digi_err <<- e$message
          NULL
        })
      }



      # # 2) Digitization (once per run)
      #
      # digi_err <- NA_character_
      # digi <- NULL
      # if (!is.null(res) && is.character(km_png) && length(km_png) == 1L && file.exists(km_png)) {
      #
      #   digi <- tryCatch({
      #     digitize_km_image(
      #       img_path = km_png,
      #       x_end = res$axis$x_end_plot,
      #       x_increment = res$axis$x_increment,
      #       y_increment = 0.2,
      #       mode = digitize_mode,
      #       num_curves = 2,
      #       bg_lightness = digitize_bg_lightness,
      #       nr_neighbors = digitize_nr_neighbors,
      #       seed = seed_run + 1L,
      #       make_plot = isTRUE(digitize_make_plot),
      #       plot_out = file.path(run_dir, paste0(prefix, "_digitized_points.png"))
      #     )
      #   }, error = function(e){
      #     digi_err <<- e$message
      #     NULL
      #   })
      #
      #   if (!is.null(digi) && is.character(digi$plot) && length(digi$plot) == 1L && file.exists(digi$plot)) {
      #     digi_png <- digi$plot
      #   }
      # }
      # # ---- 2.5) 根据真实 IPD 强制对齐数字化曲线的臂标签 ----
      # if (!is.null(res) && !is.null(digi)) {
      #
      #   # 辅助函数：从数字化数据估计 logHR 方向（正数表示 curve2 风险更高）
      #   estimate_digi_logHR <- function(dt, c1 = 1, c2 = 2, grid_n = 200, eps = 1e-6) {
      #     # dt: data.frame with columns time, St, curve
      #     # 提取两条曲线数据
      #     cdf <- dt[dt$curve == c1, c("time", "St"), drop = FALSE]
      #     tdf <- dt[dt$curve == c2, c("time", "St"), drop = FALSE]
      #
      #     # 预处理单条曲线：去重、添加 (0,1)、保证非增
      #     prep_one <- function(df) {
      #       if (nrow(df) == 0) return(NULL)
      #       df <- df[order(df$time), , drop = FALSE]
      #       # 保留每个时间点的最后一个 St（去重）
      #       df <- df[!duplicated(df$time, fromLast = TRUE), , drop = FALSE]
      #       if (!any(df$time == 0)) {
      #         df <- rbind(data.frame(time = 0, St = 1), df)
      #       }
      #       df <- df[order(df$time), , drop = FALSE]
      #       df$St <- pmin(1, pmax(0, df$St))          # 限制在 [0,1]
      #       df$St <- cummin(df$St)                     # 强制非增
      #       df
      #     }
      #
      #     cdf <- prep_one(cdf)
      #     tdf <- prep_one(tdf)
      #
      #     if (is.null(cdf) || is.null(tdf) || nrow(cdf) < 3 || nrow(tdf) < 3) return(NA_real_)
      #
      #     tmax <- min(max(cdf$time), max(tdf$time))
      #     if (!is.finite(tmax) || tmax <= 0) return(NA_real_)
      #
      #     tg <- seq(0, tmax, length.out = grid_n)
      #
      #     # 用 step 函数近似
      #     Sc <- stats::approx(cdf$time, cdf$St, xout = tg, method = "constant", f = 0, rule = 2)$y
      #     St <- stats::approx(tdf$time, tdf$St, xout = tg, method = "constant", f = 0, rule = 2)$y
      #
      #     Sc <- pmin(1 - eps, pmax(eps, Sc))
      #     St <- pmin(1 - eps, pmax(eps, St))
      #
      #     Hc <- -log(Sc)
      #     Ht <- -log(St)
      #
      #     ok <- is.finite(Hc) & is.finite(Ht) & Hc > 0 & Ht > 0
      #     if (!any(ok)) return(NA_real_)
      #
      #     beta <- log(Ht[ok] / Hc[ok])
      #     beta <- beta[is.finite(beta)]
      #     if (length(beta) == 0) return(NA_real_)
      #
      #     stats::median(beta)
      #   }
      #
      #   # 从真实 IPD 估计 logHR 方向（使用 Cox 模型）
      #   require(survival)  # 确保 survival 包已加载
      #   true_cox <- tryCatch(
      #     survival::coxph(survival::Surv(time, status) ~ arm, data = res$ipd),
      #     error = function(e) NULL
      #   )
      #   if (!is.null(true_cox) && length(stats::coef(true_cox)) == 1) {
      #     true_beta <- stats::coef(true_cox)  # 正数表示 Treatment 风险更高
      #   } else {
      #     true_beta <- NA_real_
      #   }
      #
      #   # 从数字化数据估计 logHR 方向（基于当前 curve 标签）
      #   digi_beta <- estimate_digi_logHR(digi$data, c1 = 1, c2 = 2)
      #
      #   swapped_by_truth <- FALSE
      #   if (is.finite(digi_beta) && is.finite(true_beta)) {
      #     # 如果数字化估计的方向与真实方向相反（符号相反），则需要交换
      #     if (sign(digi_beta) != sign(true_beta)) {
      #       # 交换 curve 标签
      #       digi$data$curve <- ifelse(digi$data$curve == 1, 2,
      #                                 ifelse(digi$data$curve == 2, 1, digi$data$curve))
      #       swapped_by_truth <- TRUE
      #       # 更新 meta 信息（便于调试）
      #       if (is.null(digi$meta)) digi$meta <- list()
      #       digi$meta$curve_swapped_by_truth <- TRUE
      #       digi$meta$true_logHR <- true_beta
      #       digi$meta$digi_logHR_before_swap <- digi_beta
      #     }
      #   }
      #
      #   # 可选：记录交换情况
      #   if (swapped_by_truth && verbose) {
      #     message("Run ", run_id, ": curves swapped based on truth (true beta = ", round(true_beta, 3), ")")
      #   }
      # }


      # 3) For each method: reconstruct + metrics

      for (method in methods) {
        message("  method: ", method)   # 新增

        k_row <- k_row + 1L

        # default row skeleton
        base_row <- list(
          run_id = run_id,
          scenario_id = scenario_id,
          censoring = censoring,
          target_censoring = target_censoring,
          n_control = n_control,
          n_treatment = n_treatment,
          rep = r,
          seed = seed_run,
          method = method,
          km_png = if (keep_files) km_png else NA_character_,
          digi_png = if (keep_files && isTRUE(digitize_make_plot)) digi_png else NA_character_,
          dgm_error = dgm_err,
          digitize_error = digi_err,
          reconstruct_error = NA_character_
        )

        # if DGM or digitize failed -> metrics NA
        if (is.null(res) || is.null(digi)) {
          rows[[k_row]] <- as.data.frame(base_row, stringsAsFactors = FALSE)
          next
        }

        # reconstruct
        rec_err <- NA_character_
        rec <- tryCatch({
          # 设置超时 30 秒
          R.utils::withTimeout({
            reconstruct_ipd_twoarm(
              digi = digi,
              risk_table = res$risk_table,
              method = method,
              curve_map = curve_map,
              total_events = list(
                Control = res$totals$total_events_ctrl,
                Treatment = res$totals$total_events_trt
              )
            )
          }, timeout = 30, onTimeout = "error")
        }, error = function(e){
          # 区分超时错误和其他错误
          if (inherits(e, "TimeoutException")) {
            rec_err <<- paste("Timeout after 30 sec:", e$message)
          } else {
            rec_err <<- e$message
          }
          NULL
        })
        base_row$reconstruct_error <- rec_err

        if (is.null(rec) || is.null(rec$ipd)) {
          rows[[k_row]] <- as.data.frame(base_row, stringsAsFactors = FALSE)
          next
        }

        # metrics
        tau <- res$totals$tau_study_end
        met <- tryCatch({
          evaluate_one_method(
            true_ipd = res$ipd,
            digi_dt = digi$data,
            rec_ipd = rec$ipd,
            risk_table = res$risk_table,
            tau = tau,
            curve_map = curve_map,
            do_rmst = metrics_do_rmst,
            grid_n = metrics_grid_n
          )
        }, error = function(e){
          # metrics failure should not crash the loop
          list(metrics_error = e$message)
        })

        # pack row
        out_row <- c(
          base_row,
          list(
            tau = tau,
            dgm_censor_prop = res$totals$censor_prop,
            dgm_total_events_ctrl = res$totals$total_events_ctrl,
            dgm_total_events_trt  = res$totals$total_events_trt
          ),
          met
        )

        rows[[k_row]] <- as.data.frame(out_row, stringsAsFactors = FALSE)
      }


      # 4) Cleanup (unless keep_files)

      if (!keep_files) {
        ok_to_delete <- TRUE
        if (keep_failed_files) {
          # keep files if any failure happened in DGM/digitize
          if (!is.na(dgm_err) || !is.na(digi_err)) ok_to_delete <- FALSE
        }
        if (ok_to_delete) {
          if (is.character(km_png) && file.exists(km_png)) file.remove(km_png)
          if (is.character(digi_png) && file.exists(digi_png)) file.remove(digi_png)
        }
      }

      if (verbose && (r %% max(1L, floor(n_rep / 5L)) == 0L)) {
        message("  done: setting ", i, "/", nrow(design), " | rep ", r, "/", n_rep)
      }
    }
  }

  # finalize table
  results <- bind_rows_safe(rows)

  # a compact "error view" for quick debugging
  error_rows <- results[
    (!is.na(results$dgm_error) & nzchar(results$dgm_error)) |
      (!is.na(results$digitize_error) & nzchar(results$digitize_error)) |
      (!is.na(results$reconstruct_error) & nzchar(results$reconstruct_error)) |
      ("metrics_error" %in% names(results) & !is.na(results$metrics_error) & nzchar(results$metrics_error)),
    , drop = FALSE
  ]

  list(
    results = results,
    errors = error_rows
  )
}


# --------------------------
# Minimal test example
# --------------------------

# scenario library (build once)
lib <- make_scenario_library(mean_time = 10)

# Example 1 ----

design <- data.frame(
  scenario_id = c("W_shape0p6_HR085", "W_shape0p6_HR067", "W_shape0p6_HR085", "W_shape0p6_HR067"),
  censoring = c("random","random","back","back"),
  target_censoring = c(0.30, 0.30, 0.30, 0.30),
  n_control = c(250, 250, 250, 250),
  n_treatment = c(250, 250, 250, 250),
  stringsAsFactors = FALSE
)

out_dir <- here::here("test", "0608567_2cen_100_test")


bench <- run_single_trial_benchmark(
  design = design,
  n_rep = 10,
  scenario_lib = lib,
  base_seed = 2026,
  out_dir = out_dir,
  keep_files = FALSE,
  verbose = TRUE,
  digitize_mode = "keep_censor_marks",
  digitize_make_plot = FALSE
)
write.csv(bench2$results, file = file.path(here::here("test"), "results.csv"), row.names = FALSE)





# bench4$results %>%
#   group_by(censoring) %>%
#   summarise(
#     swap_rate = mean(curve_swapped, na.rm = TRUE),
#     rmse_mean = mean(rmse_mean, na.rm = TRUE),
#     km_rmse_mean = mean(km_rmse_mean, na.rm = TRUE),
#     cox_logHR_diff = mean(cox_logHR_diff, na.rm = TRUE)
#   ) %>%
#   arrange(censoring)


# Example 2 ----

design2 <- data.frame(
  scenario_id = c("W_shape0p6_HR067", "W_shape0p6_HR067", "W_shape0p6_HR067", "W_shape0p6_HR067", "W_shape0p6_HR067"),
  censoring = c("random","informative","exp","front","back"),
  target_censoring = c(0.30, 0.30, 0.30, 0.30, 0.30),
  n_control = c(250, 250, 250, 250, 250),
  n_treatment = c(250, 250, 250, 250, 250),
  stringsAsFactors = FALSE
)

out_dir2 <- here::here("test", "numric_rep10_HR067_5cen")


bench2 <- run_single_trial_benchmark(
  design = design2,
  n_rep = 100,
  scenario_lib = lib,
  base_seed = 2,
  out_dir = out_dir2,
  keep_files = FALSE,
  verbose = TRUE,
  digitize_mode = "keep_censor_marks",
  digitize_make_plot = FALSE
)

write.csv(bench2$results, file = file.path(here::here("test"), "single_trial_06067_5cen_300.csv"), row.names = FALSE)


# Example 3 ----

design3 <- data.frame(
  scenario_id = c("W_shape0p6_HR085", "W_shape0p6_HR085", "W_shape0p6_HR085", "W_shape0p6_HR085", "W_shape0p6_HR085"),
  censoring = c("random","informative","exp","front","back"),
  target_censoring = c(0.30, 0.30, 0.30, 0.30, 0.30),
  n_control = c(250, 250, 250, 250, 250),
  n_treatment = c(250, 250, 250, 250, 250),
  stringsAsFactors = FALSE
)

out_dir3 <- here::here("test", "numric_rep100_HR085_5cen")


bench3 <- run_single_trial_benchmark(
  design = design3,
  n_rep = 100,
  scenario_lib = lib,
  base_seed = 2,
  out_dir = out_dir3,
  keep_files = FALSE,
  verbose = TRUE,
  digitize_mode = "keep_censor_marks",
  digitize_make_plot = FALSE
)

write.csv(bench3$results, file = file.path(here::here("test"), "single_trial_06085_5cen_100.csv"), row.names = FALSE)


# Example 4 ----

design4 <- data.frame(
  scenario_id = c("W_shape2p0_HR085", "W_shape2p0_HR085", "W_shape2p0_HR085", "W_shape2p0_HR085", "W_shape2p0_HR085"),
  censoring = c("random","informative","exp","front","back"),
  target_censoring = c(0.30, 0.30, 0.30, 0.30, 0.30),
  n_control = c(250, 250, 250, 250, 250),
  n_treatment = c(250, 250, 250, 250, 250),
  stringsAsFactors = FALSE
)

out_dir4 <- here::here("test", "numric_rep100_20HR085_5cen")

bench4 <- run_single_trial_benchmark(
  design = design4,
  n_rep = 100,
  scenario_lib = lib,
  base_seed = 2,
  out_dir = out_dir4,
  keep_files = FALSE,
  verbose = TRUE,
  digitize_mode = "keep_censor_marks",
  digitize_make_plot = FALSE
)

write.csv(bench4$results, file = file.path(here::here("test"), "single_trial_20085_5cen_100.csv"), row.names = FALSE)


# Example 5 ----

design5 <- data.frame(
  scenario_id = c("W_shape0p6_HR067", "W_shape0p6_HR067", "W_shape0p6_HR067", "W_shape0p6_HR067", "W_shape0p6_HR067"),
  censoring = c("random","informative","exp","front","back"),
  target_censoring = c(0.30, 0.30, 0.30, 0.30, 0.30),
  n_control = c(250, 250, 250, 250, 250),
  n_treatment = c(250, 250, 250, 250, 250),
  stringsAsFactors = FALSE
)

out_dir5 <- here::here("test", "numric_rep500_06HR067_5cen")


# 5) run small benchmark
bench5 <- run_single_trial_benchmark(
  design = design5,
  n_rep = 500,
  scenario_lib = lib,
  base_seed = 2,
  out_dir = out_dir5,
  keep_files = FALSE,
  verbose = TRUE,
  digitize_mode = "keep_censor_marks",
  digitize_make_plot = FALSE
)

write.csv(bench5$results, file = file.path(here::here("test"), "single_trial_06067_5cen_500.csv"), row.names = FALSE)




# inspect outputs ----


res <- bench5$results

# mark failures
res$ok <- is.na(res$reconstruct_error) | !nzchar(res$reconstruct_error)

# summarize by censoring x method
summary_tbl <- aggregate(
  cbind(km_rmse_mean, cox_logHR_diff) ~ censoring + method,
  data = res[res$ok & is.finite(res$km_rmse_mean) & is.finite(res$cox_logHR_diff), ],
  FUN = mean
)

fail_rate <- aggregate(
  ok ~ censoring + method,
  data = res,
  FUN = function(x) 1 - mean(x, na.rm = TRUE)
)
names(fail_rate)[names(fail_rate) == "ok"] <- "fail_rate"

summary_tbl <- merge(summary_tbl, fail_rate, by = c("censoring"), all = TRUE)
summary_tbl
summary_tbl[order(summary_tbl$censoring), ]


bench5$results %>%
  group_by(method) %>%
  summarise(
    rmse_mean = mean(rmse_mean, na.rm = TRUE),
    km_rmse_mean = mean(km_rmse_mean, na.rm = TRUE),
    cox_logHR_diff = mean(cox_logHR_diff, na.rm = TRUE)
  ) %>%
  arrange(method)

bench5$results %>%
  group_by(censoring) %>%
  summarise(
    rmse_mean = mean(rmse_mean, na.rm = TRUE),
    km_rmse_mean = mean(km_rmse_mean, na.rm = TRUE),
    cox_logHR_diff = mean(cox_logHR_diff, na.rm = TRUE)
  ) %>%
  arrange(censoring)


bench5$results %>%
  dplyr::group_by(censoring, method) %>%
  dplyr::summarise(
    n = dplyr::n(),
    err_rate = mean(!is.na(reconstruct_error)),
    rmse_mean = mean(rmse_mean, na.rm = TRUE),
    km_rmse_mean = mean(km_rmse_mean, na.rm = TRUE),
    cox_logHR_diff = mean(cox_logHR_diff, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(censoring, method)


plot(sort(bench2$results$cox_logHR_diff))
plot(sort(bench2$results$rmse_mean))
plot(sort(bench2$results$km_rmse_mean))
plot(sort(bench2$results$n_points_used))
plot(sort(bench2$results$nrisk_maxabs_treatment))
summary(bench2$results$rmse_mean)
