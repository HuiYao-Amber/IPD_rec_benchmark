# single_trial_censor_density_plot.R
# Purpose: stand-alone censor density plotting for a single simulated case
# Requires: source("single_trial_simulation_DGM_informSplit.R") (or your DGM file)

plot_censor_density_case <- function(
    scenario_id,
    scenario_lib,
    n_control = 250,
    n_treatment = 250,
    censoring = c("random","exp","front","back","informative_front","informative_back","informative"),
    target_censoring = 0.30,
    seed = 1,
    out_dir = ".",
    prefix = "demo",
    dpi = 300
){
  censoring <- match.arg(censoring)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Need ggplot2 for this plotting function. Please install.packages('ggplot2').")
  }

  # ----------------- Reproduce internal generation (to access tau_star + raw C clocks) -----------------
  if (!is.null(seed)) set.seed(seed)

  scn <- scenario_lib[[scenario_id]]
  if (is.null(scn)) stop("scenario_id not found in scenario_lib: ", scenario_id)

  latent <- generate_latent_event_times(scn, n_control, n_treatment)
  U <- make_censor_uniforms(n_control, n_treatment)

  tau_star <- calibrate_tau(
    censoring = censoring,
    target_censoring = target_censoring,
    scn = scn,
    T_control = latent$T_control,
    T_treat   = latent$T_treat,
    U = U
  )

  # raw censor clocks (generator C) and realized ipd
  CC <- make_censor_times(
    tau = tau_star,
    censoring = censoring,
    target_censoring = target_censoring,
    T_control = latent$T_control,
    T_treat   = latent$T_treat,
    U = U
  )
  R <- realize_ipd(
    tau = tau_star,
    censoring = censoring,
    target_censoring = target_censoring,
    T_control = latent$T_control,
    T_treat   = latent$T_treat,
    U = U
  )
  ipd <- R$ipd

  # ----------------- Build theoretical pdf (when available) -----------------
  x <- seq(0, tau_star, length.out = 512)

  pdf_df <- NULL
  if (censoring == "random") {
    pdf_df <- data.frame(
      x = x, y = rep(1 / tau_star, length(x)),
      label = "Theoretical pdf: Uniform(0, \u03C4)"
    )
  } else if (censoring == "exp") {
    # NOTE: must match make_censor_times()
    p0 <- max(0, min(0.95, target_censoring))
    lambda <- -log(1 - p0) / max(1e-9, tau_star)
    # truncated Exp on [0, tau_star]
    denom <- (1 - exp(-lambda * tau_star))
    y <- (lambda * exp(-lambda * x)) / denom
    pdf_df <- data.frame(
      x = x, y = y,
      label = "Theoretical pdf: Exp(\u03BB) truncated to [0, \u03C4]"
    )
  } else if (censoring == "front") {
    a <- 0.7; b <- 2.0
    y <- dbeta(x / tau_star, a, b) / tau_star
    pdf_df <- data.frame(
      x = x, y = y,
      label = sprintf("Theoretical pdf: \u03C4*Beta(%.1f, %.1f)", a, b)
    )
  } else if (censoring == "back") {
    a <- 2.0; b <- 0.7
    y <- dbeta(x / tau_star, a, b) / tau_star
    pdf_df <- data.frame(
      x = x, y = y,
      label = sprintf("Theoretical pdf: \u03C4*Beta(%.1f, %.1f)", a, b)
    )
  } else {
    # informative_* : no simple closed-form marginal pdf because depends on T
    pdf_df <- NULL
  }

  # ----------------- Empirical densities -----------------
  # generator clocks C (truncate to [0, tau_star])
  C_all <- c(as.numeric(CC$C_control), as.numeric(CC$C_treat))
  C_all <- C_all[is.finite(C_all) & C_all >= 0 & C_all <= tau_star]

  emp_gen <- NULL
  if (length(C_all) >= 5) {
    d <- stats::density(C_all, from = 0, to = tau_star, n = 512)
    emp_gen <- data.frame(x = d$x, y = d$y, src = "Empirical of generator C (truncated)")
  }

  # observed censor times (status==0), by arm
  cens_obs <- ipd[ipd$status == 0, , drop = FALSE]
  emp_obs <- NULL
  if (nrow(cens_obs) >= 5) {
    split_list <- split(cens_obs, cens_obs$arm)
    emp_obs <- do.call(rbind, lapply(split_list, function(dd){
      if (nrow(dd) < 5) return(NULL)
      d <- stats::density(dd$time, from = 0, to = tau_star, n = 512)
      data.frame(x = d$x, y = d$y, arm = unique(dd$arm))
    }))
  }

  # subtitle summary
  sub_txt <- sprintf(
    "%s | %s | \u03C4=%.2f | obs censor=%.2f (Ctrl=%d, Trt=%d)",
    scenario_id, censoring, tau_star,
    mean(ipd$status == 0),
    sum(ipd$arm == "Control"   & ipd$status == 0),
    sum(ipd$arm == "Treatment" & ipd$status == 0)
  )

  # ----------------- Plot -----------------
  p <- ggplot2::ggplot() +
    { if (!is.null(pdf_df)) ggplot2::geom_line(data = pdf_df, ggplot2::aes(x, y),
                                               linewidth = 1.1, color = "black") } +
    { if (!is.null(pdf_df)) ggplot2::geom_text(
      data = pdf_df[256, , drop = FALSE],
      ggplot2::aes(x, y, label = unique(pdf_df$label)),
      vjust = -1.0, size = 3.2, color = "black"
    ) } +
    { if (!is.null(emp_gen)) ggplot2::geom_line(data = emp_gen, ggplot2::aes(x, y),
                                                linetype = "dashed", color = "grey40") } +
    { if (!is.null(emp_obs)) ggplot2::geom_line(data = emp_obs, ggplot2::aes(x, y, color = arm),
                                                linewidth = 1.1) } +
    ggplot2::coord_cartesian(xlim = c(0, tau_star), ylim = c(0, NA)) +
    ggplot2::labs(
      x = "time t", y = "density f(t)",
      title = "Censoring profile: theoretical vs empirical",
      subtitle = sub_txt,
      caption = if (censoring %in% c("informative_front","informative_back")) {
        "Informative censoring depends on Ti=min(T,\u03C4); black theoretical pdf omitted; grey dashed=raw clocks; colored=observed censor times."
      } else {
        "Black=theoretical pdf; grey dashed=raw clocks; colored=observed censor times (status==0)."
      }
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(size = 9))

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  censor_png <- file.path(out_dir, paste0(prefix, "_", censoring, "_CENSOR_density.png"))
  ggplot2::ggsave(censor_png, p, width = 7, height = 4.5, dpi = dpi)

  invisible(list(
    plot_file = censor_png,
    tau_star = tau_star,
    ipd = ipd,
    C_raw = CC
  ))
}


# Convenience: draw all 6 mechanisms for the SAME scenario + seed
plot_censor_density_all <- function(
    scenario_id,
    scenario_lib,
    n_control = 250,
    n_treatment = 250,
    target_censoring = 0.30,
    seed = 1,
    out_dir = ".",
    prefix = "demo",
    dpi = 300
){
  mechs <- c("random","exp","front","back","informative_front","informative_back","informative")
  out <- lapply(mechs, function(cc){
    plot_censor_density_case(
      scenario_id = scenario_id,
      scenario_lib = scenario_lib,
      n_control = n_control,
      n_treatment = n_treatment,
      censoring = cc,
      target_censoring = target_censoring,
      seed = seed,
      out_dir = out_dir,
      prefix = prefix,
      dpi = dpi
    )
  })
  names(out) <- mechs
  invisible(out)
}
