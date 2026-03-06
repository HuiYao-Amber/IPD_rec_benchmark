
out_dir_demo <- "R/single_trial_case_demo0226"
lib <- make_scenario_library(mean_time = 10)
seed_run <- 1

res1 <- simulate_trial(
  scenario_id = "W_shape0p6_HR067",
  scenario_lib = lib,
  n_control = 250,
  n_treatment = 250,
  censoring = "random",
  target_censoring = 0.30,
  # censor_par = list(a = 0.7, b = 2.0, p = 1.0),
  seed = seed_run,
  out_dir = out_dir_demo,
  prefix = "W_shape0p6_HR067_c030_n250_seed1",
  plot = "km",
  add_censor_marks = FALSE
)

plot_censor_density_case(
  scenario_id = "W_shape0p6_HR067",
  scenario_lib = lib,
  n_control = 250,
  n_treatment = 250,
  censoring = "informative_front",
  target_censoring = 0.30,
  seed = seed_run,
  out_dir = out_dir_demo,
  prefix = "censor"
)

plot_censor_density_all(
  scenario_id = "W_shape0p6_HR067",
  scenario_lib = lib,
  n_control = 250,
  n_treatment = 250,
  target_censoring = 0.30,
  seed = seed_run,
  out_dir = out_dir_demo,
  prefix = "all_censor"
)



# digi <- digitize_km_image(
#   img_path = res1$km_png,
#   x_end = res1$axis$x_end_plot,
#   x_increment = res1$axis$x_increment,
#   y_increment = 0.2,
#   mode = "keep_censor_marks",
#   num_curves = 2,
#   bg_lightness = 0.3,
#   nr_neighbors = 20,
#   seed = seed_run + 1L,
#   make_plot = TRUE,
#   plot_out = file.path(out_dir_demo, "digitized_points.png")
# )


digi <- digitize_km_event_times(
  ipd_true = res1$ipd,
  sd_S = 0.04,
  seed = seed_run + 1L,
  risk_table = res1$risk_table,
  mark_merge_eps = 0.01
)


length(digi$data$time)
print(res1$risk_table)

rec1 <- reconstruct_ipd_twoarm(digi, res1$risk_table, method = "IPDfromKM",
                               curve_map = c(Control = 1, Treatment = 2))
rec2 <- reconstruct_ipd_twoarm(digi, res1$risk_table, method = "kmdata",
                               total_events = c(Control = res1$totals$total_events_ctrl,
                                                Treatment = res1$totals$total_events_trt),
                               curve_map = c(Control = 1, Treatment = 2))

# Scenario A: no marks
rec3a <- reconstruct_ipd_twoarm(digi, res1$risk_table,
                               method="KMtoIPD",
                               kmtoipd_scenario="no_marks")

# Scenario B: with marks
rec3b <- reconstruct_ipd_twoarm(digi, res1$risk_table,
                               method="KMtoIPD",
                               kmtoipd_scenario="with_marks")


sanity_plot_overlay(
  truth_ipd = res1$ipd,
  digi_dt = digi$data,
  rec_list = list(rec1 = rec1, rec2 = rec2, rec3a = rec3a, rec3b = rec3b),
  rec_names = c("IPDfromKM", "kmdata::ipd", "KMtoIPD_useMarks", "KMtoIPD_noMarks"),
  curve_map = c(Control = 1, Treatment = 2),  # 如果你发现反了就 swap
  out_file = file.path(out_dir_demo, "overlay.png")
)


(tau <- res1$totals$tau_study_end)
curve_map <- c(Control = 1, Treatment = 2)

m1 <- evaluate_one_method(res1$ipd, digi$data, rec1$ipd, res1$risk_table, tau = tau, curve_map = curve_map, do_rmst = FALSE)
m2 <- evaluate_one_method(res1$ipd, digi$data, rec2$ipd, res1$risk_table, tau = tau, curve_map = curve_map, do_rmst = FALSE)
m3a <- evaluate_one_method(res1$ipd, digi$data, rec3a$ipd, res1$risk_table, tau = tau, curve_map = curve_map, do_rmst = FALSE)
m3b <- evaluate_one_method(res1$ipd, digi$data, rec3b$ipd, res1$risk_table, tau = tau, curve_map = curve_map, do_rmst = FALSE)

as.data.frame(
  method = c("IPDfromKM", "kmdata::ipd", "KMtoIPD_useMarks", "KMtoIPD_noMarks"),
  rbind(m1, m2, m3a, m3b)
)







# ---- 画图 ----
# 固定随机种子以保证可重复性
base_seed <- 1

# 场景库（Weibull 形状 0.6，平均生存时间 10）
lib <- make_scenario_library(mean_time = 10)

# 感兴趣的场景 ID（来自您的 design2）
scenario_ids <- c("W_shape0p6_HR067", "W_shape0p6_HR085")

# 五种删失类型
censoring_types <- c("random", "informative", "exp", "front", "back")

# 样本量
n_control <- 250
n_treatment <- 250
target_censoring <- 0.30

# 存储所有模拟数据的列表
all_data <- list()

for (hr_id in scenario_ids) {
  for (cens in censoring_types) {
    # 模拟一次试验（每次使用相同种子，使得潜在事件时间相同，仅删失机制不同）
    set.seed(base_seed)
    sim <- simulate_trial(
      scenario_id = hr_id,
      scenario_lib = lib,
      n_control = n_control,
      n_treatment = n_treatment,
      censoring = cens,
      target_censoring = target_censoring,
      seed = base_seed,
      out_dir = NULL,
      plot = "none"
    )

    # 添加标识列
    dat <- sim$ipd
    dat$HR <- hr_id
    dat$censoring <- cens
    all_data[[paste(hr_id, cens, sep = "_")]] <- dat
  }
}

# 合并所有数据
plot_df <- bind_rows(all_data)


# 生成时间网格（0 到 30，步长 0.1）
time_grid <- seq(0, 30, length.out = 300)

theoretical_curves <- data.frame()
for (hr_id in scenario_ids) {
  scn <- lib[[hr_id]]
  # Weibull 生存函数 S(t) = exp(-(t/scale)^shape)
  S_control <- exp(-(time_grid / scn$control$scale)^scn$control$shape)
  S_treat   <- exp(-(time_grid / scn$treat$scale)^scn$treat$shape)

  theo <- rbind(
    data.frame(HR = hr_id, arm = "Control", time = time_grid, surv = S_control),
    data.frame(HR = hr_id, arm = "Treatment", time = time_grid, surv = S_treat)
  )
  theoretical_curves <- rbind(theoretical_curves, theo)
}


# 对每个组合、每个臂分别拟合 KM
km_points <- plot_df %>%
  group_by(HR, censoring, arm) %>%
  summarise(
    fit = list(survfit(Surv(time, status) ~ 1)),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    times = list(c(0, fit$time)),          # 添加起始点 time=0
    surv = list(c(1, fit$surv))            # 起始生存=1
  ) %>%
  unnest(cols = c(times, surv)) %>%
  ungroup()


# 将 HR 列转换为更友好的标签
km_points$HR_label <- factor(km_points$HR,
                             levels = scenario_ids,
                             labels = c("HR = 0.67", "HR = 0.85"))
theoretical_curves$HR_label <- factor(theoretical_curves$HR,
                                      levels = scenario_ids,
                                      labels = c("HR = 0.67", "HR = 0.85"))

# 删失类型顺序和标签
cens_labels <- c("random" = "Random", "informative" = "Informative",
                 "exp" = "Exponential", "front" = "Front", "back" = "Back")

p <- ggplot() +
  geom_step(data = km_points,
            aes(x = times, y = surv, color = arm),
            direction = "hv", size = 1) +
  geom_line(data = theoretical_curves,
            aes(x = time, y = surv, linetype = arm),
            color = "black", alpha = 0.5) +
  scale_color_manual(values = c("Control" = "#1f77b4", "Treatment" = "#d62728")) +
  scale_linetype_manual(values = c("Control" = "dashed", "Treatment" = "dotted")) +
  labs(x = "Time", y = "Survival Probability",
       color = "Observed arm", linetype = "Theoretical arm") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "grey90")) +
  facet_grid(censoring ~ HR_label,
             labeller = labeller(censoring = cens_labels))

print(p)
