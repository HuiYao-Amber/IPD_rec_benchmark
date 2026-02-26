plot_censor_density_profile_DGM <- function(ipd, R, scenario_id, censoring, tau_star,
                                            target_censoring,
                                            out_dir, prefix, dpi = 300) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  censor_png <- file.path(out_dir, paste0(prefix, "_", censoring, "_CENSOR_density.png"))
  x <- seq(0, tau_star, length.out = 512)

  # --- theoretical pdf (when available/meaningful) ---
  pdf_df <- NULL
  if (censoring == "random") {
    pdf_df <- data.frame(
      x = x, y = rep(1 / tau_star, length(x)),
      label = sprintf("Theoretical pdf: Uniform(0, τ=%.2f)", tau_star)
    )
  }
  if (censoring == "exp") {
    # DGM exp: lambda = -log(1-p0)/tau with p0=target_censoring
    p0 <- max(0, min(0.95, target_censoring))
    lambda <- -log(1 - p0) / max(1e-9, tau_star)
    pdf_df <- data.frame(
      x = x, y = lambda * exp(-lambda * x),
      label = sprintf("Theoretical pdf: Exp(λ=%.3g), p0=%.2f", lambda, p0)
    )
  }
  if (censoring == "front") {
    a <- 0.7; b <- 2.0
    pdf_df <- data.frame(
      x = x,
      y = dbeta(x / tau_star, a, b) / tau_star,
      label = sprintf("Theoretical pdf: τ*Beta(%.1f,%.1f)", a, b)
    )
  }
  if (censoring == "back") {
    a <- 2.0; b <- 0.7
    pdf_df <- data.frame(
      x = x,
      y = dbeta(x / tau_star, a, b) / tau_star,
      label = sprintf("Theoretical pdf: τ*Beta(%.1f,%.1f)", a, b)
    )
  }
  # informative: no simple closed-form -> skip

  # --- generator clocks empirical density (truncate to [0, τ]) ---
  C_all <- c(as.numeric(R$C_control), as.numeric(R$C_treat))
  C_all <- C_all[is.finite(C_all) & C_all >= 0 & C_all <= tau_star]
  emp_gen <- if (length(C_all) >= 5) {
    d <- density(C_all, from = 0, to = tau_star, n = 512)
    data.frame(x = d$x, y = d$y, src = "Empirical of generator C (truncated)")
  } else NULL

  # --- observed censor times (status==0), by arm ---
  cens_obs <- ipd[ipd$status == 0, ]
  emp_obs <- if (nrow(cens_obs) >= 5) {
    rbindlist(lapply(split(cens_obs, cens_obs$arm), function(dd) {
      d <- density(dd$time, from = 0, to = tau_star, n = 512)
      data.frame(x = d$x, y = d$y, arm = unique(dd$arm))
    }), use.names = TRUE)
  } else NULL

  # --- subtitle stats (same spirit as your old code) ---
  sub_txt <- sprintf(
    "%s | %s | τ=%.2f | obs censor=%.2f (Ctrl=%d, Trt=%d) | events(Ctrl=%d, Trt=%d) | medTime(Ctrl=%.1f, Trt=%.1f)",
    scenario_id, censoring, tau_star,
    mean(ipd$status == 0),
    sum(ipd$arm == "Control"   & ipd$status == 0),
    sum(ipd$arm == "Treatment" & ipd$status == 0),
    sum(ipd$arm == "Control"   & ipd$status == 1),
    sum(ipd$arm == "Treatment" & ipd$status == 1),
    median(ipd$time[ipd$arm == "Control"]),
    median(ipd$time[ipd$arm == "Treatment"])
  )

  cap <- switch(
    censoring,
    random      = "Black: Uniform(0,τ) theoretical pdf; grey dashed: empirical generator clocks; colored: observed censor times by arm.",
    exp         = "Black: Exp pdf; grey dashed: empirical generator clocks truncated to [0,τ]; colored: observed censor times by arm (includes admin at τ).",
    front       = "Black: τ*Beta(0.7,2.0) theoretical pdf; grey dashed: empirical generator clocks; colored: observed censor times by arm.",
    back        = "Black: τ*Beta(2.0,0.7) theoretical pdf; grey dashed: empirical generator clocks; colored: observed censor times by arm.",
    informative = "Grey dashed: empirical generator clocks; colored: observed censor times by arm. Informative censoring is a mixture depending on T; no simple closed-form pdf."
  )

  p <- ggplot() +
    { if (!is.null(pdf_df)) geom_line(data = pdf_df, aes(x, y), linewidth = 1.1, color = "black") } +
    { if (!is.null(pdf_df)) geom_text(data = pdf_df[256, ],
                                      aes(x, y, label = unique(pdf_df$label)),
                                      vjust = -1.0, size = 3.2, color = "black") } +
    { if (!is.null(emp_gen)) geom_line(data = emp_gen, aes(x, y),
                                       linetype = "dashed", color = "grey40") } +
    { if (!is.null(emp_obs)) geom_line(data = emp_obs, aes(x, y, color = arm),
                                       linewidth = 1.1) } +
    scale_color_manual(values = c(Control = "#1f77b4", Treatment = "#d62728"),
                       name = "Observed censor density") +
    coord_cartesian(xlim = c(0, tau_star), ylim = c(0, NA)) +
    labs(x = "time t", y = "density f(t)",
         title = "Censoring profile: theoretical vs empirical",
         subtitle = sub_txt,
         caption = cap) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          plot.caption = element_text(size = 9))

  ggsave(censor_png, p, width = 7, height = 4.5, dpi = dpi)
  invisible(censor_png)
}


# -------------------------
# Run Task 2: one fixed truth, compare all censoring mechanisms
# -------------------------
run_task2_censor_plots_DGM <- function(
    scenario_id = "W_shape0p6_HR067",
    n_control = 250,
    n_treatment = 250,
    target_censoring = 0.30,
    seed = 1,
    out_dir = "task2_censor_plots",
    prefix = "task2",
    censor_set = c("random", "exp", "front", "back", "informative"),
    dpi = 300
){
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  set.seed(seed)
  scenario_lib <- make_scenario_library(mean_time = 10)
  scn <- scenario_lib[[scenario_id]]
  if (is.null(scn)) stop("scenario_id not found: ", scenario_id)

  # Fix truth ONCE
  latent <- generate_latent_event_times(scn, n_control, n_treatment)
  U <- make_censor_uniforms(n_control, n_treatment)

  # Loop censoring mechanisms (same latent T and same U streams)
  out_paths <- character(0)

  for (censoring in censor_set) {
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

    out_paths <- c(out_paths, plot_censor_density_profile_DGM(
      ipd = ipd, R = R,
      scenario_id = scenario_id,
      censoring = censoring,
      tau_star = tau_star,
      target_censoring = target_censoring,
      out_dir = out_dir,
      prefix = prefix,
      dpi = dpi
    ))
  }

  invisible(out_paths)
}

# ---- run it ----
run_task2_censor_plots_DGM(
  scenario_id = "W_shape0p6_HR067",
  n_control = 250,
  n_treatment = 250,
  target_censoring = 0.30,
  seed = 20260220,
  out_dir = "task2_censor_plots",
  prefix = "singletrial_truthfixed",
  dpi = 300
)










# =========================
# Task 2 (revised): Density plots for 7 censoring mechanisms
#   Drop: Non-informative exponential-shaped (TruncExp on [0,t_max])
#   Keep: Uniform (NI), Exp+Admin (NI), Front (NI), Back (NI),
#         Front (I), Back (I), Exp+Admin (I, frailty-dependent)
# =========================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

set.seed(1)

# ---- knobs you should align with single_trial_simulation_DGM ----
n_calib   <- 200000
n_plot    <- 200000
t_max     <- 30
c_target  <- 0.30

# event-time DGM (placeholder; swap to match your DGM)
weib_shape <- 0.6
weib_scale <- 12   # <-- adjust if needed

# informative severity (stress-test)
sigma_u <- 1.0
gamma_u <- 1.5     # higher U => earlier censoring via higher rate / time-scaling

# base Beta shapes (extreme)
beta_front_ab <- c(0.3, 4.0)
beta_back_ab  <- c(4.0, 0.3)

# ---- stable draws for calibration/root finding ----
make_draws <- function(n) {
  list(
    uT = runif(n),
    uC = runif(n),
    bF = rbeta(n, beta_front_ab[1], beta_front_ab[2]),
    bB = rbeta(n, beta_back_ab[1], beta_back_ab[2]),
    zU = rnorm(n)
  )
}

# ---- event generator: Weibull with shared frailty exp(U) ----
gen_T_weibull_frailty <- function(draws, shape, scale, U = NULL) {
  u <- draws$uT
  if (is.null(U)) {
    return(scale * (-log(u))^(1 / shape))
  } else {
    return(scale * ((-log(u)) / exp(U))^(1 / shape))
  }
}

# ---- censor generators ----
C_uniform <- function(draws, tau) {
  draws$uC * tau
}

C_exp_admin <- function(draws, lambda, t_max) {
  C <- -log(1 - draws$uC) / lambda
  pmin(C, t_max)
}

# front/back via Beta + power (p>=1 makes it more extreme)
C_beta_front <- function(draws, p, t_max) {
  t_max * (draws$bF^p)
}
C_beta_back <- function(draws, p, t_max) {
  t_max * (1 - (1 - draws$bB)^p)
}

# ---- informative versions ----
C_beta_front_inf <- function(draws, p, t_max, U, gamma_u) {
  C0 <- C_beta_front(draws, p, t_max)
  pmax(0, pmin(t_max, C0 * exp(-gamma_u * U)))
}
C_beta_back_inf <- function(draws, p, t_max, U, gamma_u) {
  C0 <- C_beta_back(draws, p, t_max)
  pmax(0, pmin(t_max, C0 * exp(-gamma_u * U)))
}

# Informative exponential censoring (Admin-truncated), rate depends on frailty
C_exp_admin_inf <- function(draws, lambda, t_max, U, gamma_u) {
  rate <- lambda * exp(gamma_u * U)
  C <- -log(1 - draws$uC) / rate
  pmin(C, t_max)
}

# ---- 1D calibration via bracketing + uniroot ----
calibrate_theta <- function(obj_fun, lower, upper, expand = 2, max_expand = 12) {
  fL <- obj_fun(lower)
  fU <- obj_fun(upper)
  k <- 0
  while (sign(fL) == sign(fU) && k < max_expand) {
    lower <- lower / expand
    upper <- upper * expand
    fL <- obj_fun(lower)
    fU <- obj_fun(upper)
    k <- k + 1
  }
  if (sign(fL) == sign(fU)) {
    stop("Calibration failed: could not bracket root. Try different initial bounds.")
  }
  uniroot(function(x) obj_fun(x), lower = lower, upper = upper, tol = 1e-4)$root
}

# ---- calibration draws ----
draws_cal <- make_draws(n_calib)

# Non-informative event times
T_noninf <- gen_T_weibull_frailty(draws_cal, weib_shape, weib_scale, U = NULL)

# Informative frailty + event times
U_inf <- sigma_u * draws_cal$zU
T_inf <- gen_T_weibull_frailty(draws_cal, weib_shape, weib_scale, U = U_inf)

# ---- calibrate each mechanism to achieve P(C < T) ~ c_target ----
# 1) Uniform (NI): tau
tau_hat <- calibrate_theta(
  obj_fun = function(tau) mean(C_uniform(draws_cal, tau) < T_noninf) - c_target,
  lower = 1e-3, upper = t_max
)

# 2) Exponential + admin truncation (NI): lambda
lambda_admin_hat <- calibrate_theta(
  obj_fun = function(lam) mean(C_exp_admin(draws_cal, lam, t_max) < T_noninf) - c_target,
  lower = 1e-6, upper = 10
)

# 3) Front (NI): power p
p_front_hat <- calibrate_theta(
  obj_fun = function(p) mean(C_beta_front(draws_cal, p, t_max) < T_noninf) - c_target,
  lower = 0.5, upper = 20
)
p_front_hat <- max(p_front_hat, 1)

# 4) Back (NI): power p
p_back_hat <- calibrate_theta(
  obj_fun = function(p) mean(C_beta_back(draws_cal, p, t_max) < T_noninf) - c_target,
  lower = 0.5, upper = 20
)
p_back_hat <- max(p_back_hat, 1)

# 5) Front (I): power p
p_front_inf_hat <- calibrate_theta(
  obj_fun = function(p) mean(C_beta_front_inf(draws_cal, p, t_max, U_inf, gamma_u) < T_inf) - c_target,
  lower = 0.5, upper = 20
)
p_front_inf_hat <- max(p_front_inf_hat, 1)

# 6) Back (I): power p
p_back_inf_hat <- calibrate_theta(
  obj_fun = function(p) mean(C_beta_back_inf(draws_cal, p, t_max, U_inf, gamma_u) < T_inf) - c_target,
  lower = 0.5, upper = 20
)
p_back_inf_hat <- max(p_back_inf_hat, 1)

# 7) Exponential + admin truncation (I): lambda (frailty-dependent rate)
lambda_admin_inf_hat <- calibrate_theta(
  obj_fun = function(lam) mean(C_exp_admin_inf(draws_cal, lam, t_max, U_inf, gamma_u) < T_inf) - c_target,
  lower = 1e-6, upper = 10
)

params <- list(
  "Uniform (NI)"                       = list(type="uniform", tau=tau_hat),
  "Exponential+Admin (NI)"             = list(type="exp_admin", lambda=lambda_admin_hat),
  "Front (NI): Beta+Power"             = list(type="front", p=p_front_hat),
  "Back (NI): Beta+Power"              = list(type="back",  p=p_back_hat),
  "Front (I): Beta+Power+Frailty"      = list(type="front_inf", p=p_front_inf_hat),
  "Back (I): Beta+Power+Frailty"       = list(type="back_inf",  p=p_back_inf_hat),
  "Exponential+Admin (I): FrailtyRate" = list(type="exp_admin_inf", lambda=lambda_admin_inf_hat)
)

# ---- generate plot data ----
draws_plot <- make_draws(n_plot)
U_plot <- sigma_u * draws_plot$zU
T_plot_noninf <- gen_T_weibull_frailty(draws_plot, weib_shape, weib_scale, U = NULL)
T_plot_inf    <- gen_T_weibull_frailty(draws_plot, weib_shape, weib_scale, U = U_plot)

gen_C_by_param <- function(nm, par) {
  switch(par$type,
         uniform       = C_uniform(draws_plot, par$tau),
         exp_admin     = C_exp_admin(draws_plot, par$lambda, t_max),
         front         = C_beta_front(draws_plot, par$p, t_max),
         back          = C_beta_back(draws_plot, par$p, t_max),
         front_inf     = C_beta_front_inf(draws_plot, par$p, t_max, U_plot, gamma_u),
         back_inf      = C_beta_back_inf(draws_plot, par$p, t_max, U_plot, gamma_u),
         exp_admin_inf = C_exp_admin_inf(draws_plot, par$lambda, t_max, U_plot, gamma_u),
         stop("Unknown type")
  )
}

dt <- rbindlist(lapply(names(params), function(nm) {
  C <- gen_C_by_param(nm, params[[nm]])
  T_use <- if (grepl("\\(I\\)", nm)) T_plot_inf else T_plot_noninf
  data.table(
    mechanism = nm,
    C = C,
    censored = as.integer(C < T_use),
    is_admin = as.integer(C >= t_max - 1e-12)
  )
}))

sanity <- dt[, .(
  achieved_cens = mean(censored),
  p_admin = mean(is_admin)
), by = mechanism][order(mechanism)]
print(sanity)

# ---- density plot: show continuous part; annotate P(C=t_max) ----
mass <- dt[, .(p_admin = mean(is_admin)), by = mechanism]

p <- ggplot(dt[C < t_max - 1e-10], aes(x = C)) +
  geom_density() +
  facet_wrap(~ mechanism, ncol = 2, scales = "free_y") +
  geom_text(
    data = mass,
    aes(x = 0.65 * t_max, y = Inf, label = sprintf("P(C=t_max)=%.2f", p_admin)),
    inherit.aes = FALSE, vjust = 1.2, size = 3
  ) +
  labs(
    title = sprintf("Censoring-time density (7 mechanisms; target P(C<T)=%.2f, t_max=%.1f)", c_target, t_max),
    x = "Censoring time C",
    y = "Density"
  ) +
  theme_bw()
print(p)

# Optional overlay (continuous parts)
p2 <- ggplot(dt[C < t_max - 1e-10], aes(x = C, color = mechanism)) +
  geom_density() +
  labs(title = "Overlay of continuous censoring-time densities (continuous part)", x = "C", y = "Density") +
  theme_bw() +
  theme(legend.position = "bottom")
print(p2)
