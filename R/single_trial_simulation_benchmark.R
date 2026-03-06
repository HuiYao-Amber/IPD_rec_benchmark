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
    keep_failed_files = TRUE,    # kept for compatibility (numeric digitization doesn't save extra files)
    verbose = TRUE,

    # ---- digitization knobs (kept for compatibility; numeric digitization uses sd_S + risk_table only) ----
    digitize_mode = "keep_censor_marks",
    digitize_bg_lightness = 0.4,
    digitize_nr_neighbors = 20,
    digitize_sd_S = 0.04,
    digitize_make_plot = FALSE,

    # ---- reconstruction knobs ----
    kmdata_interpolate = FALSE,

    # ---- metrics knobs ----
    metrics_do_rmst = TRUE,
    metrics_grid_n = 300,

    # NEW (2026-02):
    max_tries_per_rep = 50,
    timeout_sec = 30
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
    x <- lapply(x, function(z) as.data.frame(z, stringsAsFactors = FALSE))
    all_cols <- unique(unlist(lapply(x, names)))
    x <- lapply(x, function(z){
      for (cc in setdiff(all_cols, names(z))) z[[cc]] <- NA
      z[, all_cols, drop = FALSE]
    })
    do.call(rbind, x)
  }

  # helper: safe timeout wrapper (no hard dependency on R.utils)
  run_with_timeout <- function(expr_fun, timeout_sec){
    if (requireNamespace("R.utils", quietly = TRUE)) {
      return(R.utils::withTimeout(expr_fun(), timeout = timeout_sec, onTimeout = "error"))
    }
    expr_fun()
  }

  # NEW (2026-02): expand KMtoIPD into two scenarios for benchmark
  expand_methods <- function(methods){
    jobs <- list()
    for (m in methods) {
      if (m == "KMtoIPD") {
        jobs[[length(jobs) + 1L]] <- list(method_base = "KMtoIPD", kmtoipd_scenario = "no_marks",  method = "KMtoIPD_no_marks")
        jobs[[length(jobs) + 1L]] <- list(method_base = "KMtoIPD", kmtoipd_scenario = "with_marks", method = "KMtoIPD_with_marks")
      } else {
        jobs[[length(jobs) + 1L]] <- list(method_base = m, kmtoipd_scenario = NA_character_, method = m)
      }
    }
    jobs
  }
  jobs <- expand_methods(methods)

  # NEW: fail-fast order to reduce runtime in hard scenarios
  priority <- c("kmdata", "IPDfromKM", "KMtoIPD_no_marks", "KMtoIPD_with_marks")
  job_names <- vapply(jobs, function(z) z$method, character(1))
  ord_ff <- order(match(job_names, priority, nomatch = 999L))
  jobs_ff <- jobs[ord_ff]

  if (verbose) {
    message("Benchmark starts: ", nrow(design), " settings x ", n_rep, " successful reps; methods: ",
            paste(job_names, collapse = ", "))
  }

  rows <- list()
  errors <- list()

  # Main loops
  for (i in seq_len(nrow(design))) {

    scenario_id      <- as.character(design$scenario_id[i])
    censoring        <- as.character(design$censoring[i])
    target_censoring <- as.numeric(design$target_censoring[i])
    n_control        <- as.integer(design$n_control[i])
    n_treatment      <- as.integer(design$n_treatment[i])

    if (verbose) {
      message("Setting ", i, "/", nrow(design), ": ", scenario_id, " | censoring=", censoring,
              " | target=", target_censoring, " | n=", n_control, ",", n_treatment)
    }

    for (r in seq_len(n_rep)) {

      success <- FALSE

      # NEW (2026-02): retry within each rep until all reconstructions succeed
      for (try_idx in seq_len(max_tries_per_rep)) {

        seed_run <- base_seed + i * 100000L + r * 1000L + try_idx

        run_id <- paste0(
          scenario_id, "_", censoring, "_c", sprintf("%03d", round(100 * target_censoring)),
          "_n", n_control, "v", n_treatment, "_rep", sprintf("%03d", r),
          "_try", sprintf("%02d", try_idx)
        )

        run_dir <- if (keep_files) file.path(tmp_root, scenario_id, censoring, paste0("rep_", r), paste0("try_", try_idx)) else tmp_root
        if (!dir.exists(run_dir)) dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

        km_png <- NA_character_

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
            prefix = run_id,
            plot = "km",
            add_censor_marks = FALSE
          )
        }, error = function(e){
          dgm_err <<- e$message
          NULL
        })
        if (!is.null(res) && !is.null(res$km_png)) km_png <- res$km_png

        if (is.null(res)) {
          errors[[length(errors) + 1L]] <- data.frame(
            run_id = run_id, stage = "DGM", method = NA_character_, error = dgm_err,
            stringsAsFactors = FALSE
          )
          # cleanup for failed try
          if (!keep_files && is.character(km_png) && file.exists(km_png)) file.remove(km_png)
          next
        }

        # 2) Digitization (numeric mode; includes risk-table times + censor marks)
        digi_err <- NA_character_
        digi <- tryCatch({
          digitize_km_event_times(
            ipd_true = res$ipd,
            sd_S = digitize_sd_S,
            seed = seed_run + 1L,
            risk_table = res$risk_table # NEW (2026-02): ensures risk-table time has St
          )
        }, error = function(e){
          digi_err <<- e$message
          NULL
        })

        if (is.null(digi)) {
          errors[[length(errors) + 1L]] <- data.frame(
            run_id = run_id, stage = "Digitize", method = NA_character_, error = digi_err,
            stringsAsFactors = FALSE
          )
          if (!keep_files && is.character(km_png) && file.exists(km_png)) file.remove(km_png)
          next
        }

        # 3) Reconstructions (fail-fast order). If any fails -> discard this try.
        rec_list <- list()
        rec_fail <- FALSE
        rec_fail_msg <- NA_character_
        rec_fail_method <- NA_character_

        for (job in jobs_ff) {

          rec_err <- NA_character_
          rec <- tryCatch({
            run_with_timeout(function(){
              reconstruct_ipd_twoarm(
                digi = digi,
                risk_table = res$risk_table,
                method = job$method_base,
                curve_map = curve_map,
                total_events = list(
                  Control = res$totals$total_events_ctrl,
                  Treatment = res$totals$total_events_trt
                ),
                kmtoipd_scenario = if (is.na(job$kmtoipd_scenario)) "no_marks" else job$kmtoipd_scenario,
                kmdata_interpolate = kmdata_interpolate
              )
            }, timeout_sec = timeout_sec)
          }, error = function(e){
            rec_err <<- e$message
            NULL
          })

          if (is.null(rec) || is.null(rec$ipd)) {
            rec_fail <- TRUE
            rec_fail_method <- job$method
            rec_fail_msg <- rec_err
            errors[[length(errors) + 1L]] <- data.frame(
              run_id = run_id, stage = "Reconstruct", method = job$method, error = rec_err,
              stringsAsFactors = FALSE
            )
            break
          }

          rec_list[[job$method]] <- rec
        }

        if (rec_fail) {
          # cleanup for failed try
          if (!keep_files && is.character(km_png) && file.exists(km_png)) file.remove(km_png)
          next
        }

        # If we got here, ALL reconstructions succeeded -> compute metrics + store rows
        tau <- res$totals$tau_study_end

        for (job in jobs) {
          method_label <- job$method
          rec <- rec_list[[method_label]]

          base_row <- list(
            run_id = run_id,
            scenario_id = scenario_id,
            censoring = censoring,
            target_censoring = target_censoring,
            n_control = n_control,
            n_treatment = n_treatment,
            rep = r,
            try = try_idx,
            seed = seed_run,
            method = method_label,
            method_base = job$method_base,
            kmtoipd_scenario = job$kmtoipd_scenario,
            km_png = if (keep_files) km_png else NA_character_,
            dgm_error = NA_character_,
            digitize_error = NA_character_,
            reconstruct_error = NA_character_
          )

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
            list(metrics_error = e$message)
          })

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

          rows[[length(rows) + 1L]] <- out_row
        }

        # 4) Cleanup (unless keep_files)
        if (!keep_files) {
          if (is.character(km_png) && file.exists(km_png)) file.remove(km_png)
        }

        success <- TRUE
        if (verbose) message("  rep ", r, " success (try ", try_idx, ")")
        break
      }

      if (!success && verbose) {
        message("  rep ", r, " failed after ", max_tries_per_rep, " tries; discarded")
      }
    }
  }

  results <- bind_rows_safe(rows)
  error_rows <- bind_rows_safe(errors)

  list(
    results = results,
    errors = error_rows
  )
}


# --------------------------
# Minimal test example (will NOT run unless you set if(FALSE) -> if(TRUE))
# --------------------------
lib <- make_scenario_library(mean_time = 10)
out_dir_bench <- "benchmark_results_0226"

design <- data.frame(
  scenario_id = c("W_shape0p6_HR067","W_shape0p6_HR067","W_shape0p6_HR067","W_shape0p6_HR067","W_shape0p6_HR067",
                  "W_shape2p0_HR067","W_shape2p0_HR067","W_shape2p0_HR067","W_shape2p0_HR067","W_shape2p0_HR067"),
  censoring = c("front","random","exp","informative","back",
                "front","random","exp","informative","back"),
  target_censoring = c(0.30),
  n_control = c(200),
  n_treatment = c(200),
  stringsAsFactors = FALSE
)

bench <- run_single_trial_benchmark(
  design = design,
  n_rep = 300,
  scenario_lib = lib,
  base_seed = 1,
  out_dir = out_dir_bench,
  keep_files = FALSE,
  verbose = TRUE,
  digitize_sd_S = 0.04,
  kmdata_interpolate = TRUE,
  max_tries_per_rep = 30
)

write.csv(bench$results, file.path(out_dir_bench, "benchmark_results.csv"), row.names = FALSE)

colnames(bench$results)
table(bench$results$scenario_id)
table(bench$results$method)
head(bench$errors)


res <- bench$results %>% filter(scenario_id == "W_shape0p6_HR067")


# 期望的方法集合（以结果中实际出现的为准）
methods_in_res <- sort(unique(res$method))
print(methods_in_res)
sort(unique(res$scenario_id))
sort(unique(res$target_censoring))

# 每个 (scenario + censoring + target + n + rep) 的方法数
key <- with(res, paste(scenario_id, censoring, target_censoring, n_control, n_treatment, rep, sep="|"))
n_methods_by_rep <- tapply(res$method, key, function(x) length(unique(x)))

summary(n_methods_by_rep)

# 找出“不完整”的 rep（方法数不足）
expected_k <- length(methods_in_res)
bad_rep <- names(n_methods_by_rep)[n_methods_by_rep != expected_k]
length(bad_rep)
if (length(bad_rep)) {
  message("Incomplete reps found: ", length(bad_rep))
  print(head(bad_rep, 10))
}


try_unique_by_rep <- tapply(res$try, key, function(x) length(unique(x)))
table(try_unique_by_rep)

# 如果出现 >1，说明同一rep混入了不同try（不应该）
bad_try <- names(try_unique_by_rep)[try_unique_by_rep != 1]
if (length(bad_try)) {
  message("Rep with multiple try values: ", length(bad_try))
  print(head(bad_try, 10))
}


dup_flag <- duplicated(res[, c("run_id","method")])
sum(dup_flag)

if (any(dup_flag)) {
  message("Found duplicated (run_id, method) rows:")
  print(res[dup_flag, c("run_id","method","scenario_id","rep","try")][1:10, ])
}

