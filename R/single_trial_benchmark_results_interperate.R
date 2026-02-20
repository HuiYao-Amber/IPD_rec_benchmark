library(dplyr)
library(tidyr)

# write.csv(bench5$results, file = file.path(here::here("test"), "single_trial_06067_5cen_500.csv"), row.names = FALSE)
# res <- read.csv(file.path(here::here("test"), "single_trial_06067_5cen_500.csv"))

res <- bench5$results

dim(res)


#检查稳定性
res %>%
  mutate(
    ok_dgm = is.na(dgm_error),
    ok_digi = is.na(digitize_error)
  ) %>%
  group_by(censoring) %>%
  summarise(
    n = n(),
    ok_dgm_rate = mean(ok_dgm),
    ok_digi_rate = mean(ok_digi),
    .groups = "drop"
  ) %>%
  arrange(censoring)

res %>%
  mutate(
    ok_rec = is.na(reconstruct_error)
  ) %>%
  group_by(method) %>%
  summarise(
    n = n(),
    ok_rec_rate = mean(ok_rec),
    .groups = "drop"
  ) %>%
  arrange(method)



res %>%
  mutate(
    ok_dgm = is.na(dgm_error),
    ok_digi = is.na(digitize_error),
    ok_rec = is.na(reconstruct_error),
    ok_all = ok_dgm & ok_digi & ok_rec
  ) %>%
  group_by(censoring, method) %>%
  summarise(
    n = n(),
    ok_rate = mean(ok_all),
    rec_fail = mean(!ok_rec),
    .groups = "drop"
  ) %>%
  arrange(censoring, method,rec_fail)


# 检查删失数量
res %>%
  group_by(censoring) %>%
  summarise(
    censor_prop = mean(dgm_censor_prop, na.rm = TRUE),
    ev_ctrl = mean(dgm_total_events_ctrl, na.rm = TRUE),
    ev_trt  = mean(dgm_total_events_trt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(censoring)


# km_rmse_mean 的分布（只看成功重建的）
res_ok <- res %>%
  filter(is.na(reconstruct_error), is.finite(km_rmse_mean))

res_ok %>%
  group_by(censoring, method) %>%
  summarise(
    mean_km_rmse = mean(km_rmse_mean),
    sd_km_rmse   = sd(km_rmse_mean),
    q50 = quantile(km_rmse_mean, 0.50),
    q90 = quantile(km_rmse_mean, 0.90),
    .groups = "drop"
  ) %>%
  arrange(censoring, method)

ggplot(res_ok, aes(
  x = factor(censoring, levels = c("random", "exp", "front", "back","informative")),
  y = km_rmse_mean,
  fill = factor(method, levels = c( "KMtoIPD","IPDfromKM", "kmdata"))
)) +
  geom_boxplot(outlier.alpha = 0.2) +
  coord_cartesian(ylim = c(0, 0.1)) +
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(
    title = "KM RMSE (reconstructed vs true) across censoring types",
    x = "Censoring mechanism",
    y = "KM RMSE (mean of arms)"
  )


# 可选保存
# ggsave("km_rmse_by_censoring.png", p, width = 8, height = 5, dpi = 200)


# cox_logHR_diff 的分布（只看成功重建的）
res_ok %>%
  group_by(censoring, method) %>%
  summarise(
    mean_beta_diff = mean(cox_logHR_diff, na.rm = TRUE),
    hr_ratio = exp(mean_beta_diff),  # HR_rec / HR_true on average
    .groups = "drop"
  ) %>%
  arrange(censoring, method)

ggplot(res_ok, aes(
  x = factor(censoring, levels = c("random", "exp", "front", "back","informative")),
  y = cox_logHR_diff,
  fill = factor(method, levels = c( "KMtoIPD","IPDfromKM", "kmdata"))
)) +
  geom_boxplot(outlier.alpha = 0.2) +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(
    title = "Cox logHR difference (reconstructed - true)",
    x = "Censoring mechanism",
    y = "Δ logHR"
  )



# rmst_diff 的分布（只看成功重建的）
res_ok %>%
  group_by(censoring, method) %>%
  summarise(
    mean_rmst_diff = mean(rmst_diff, na.rm = TRUE),
    sd_rmst   = sd(rmst_diff, na.rm = TRUE),
    q50 = quantile(rmst_diff, 0.50),
    q90 = quantile(rmst_diff, 0.90),
    .groups = "drop"
  ) %>%
  arrange(censoring, method)

ggplot(res_ok, aes(
  x = factor(censoring, levels = c("random", "exp", "front", "back","informative")),
  y = rmst_diff,
  fill = factor(method, levels = c( "KMtoIPD","IPDfromKM", "kmdata"))
)) +
  geom_boxplot(outlier.alpha = 0.2) +
  coord_cartesian(ylim = c(-3, 3)) +
  theme_bw() +
  theme(legend.position = "bottom")+
  labs(
    title = "Restricted Mean Survival Time",
    x = "Censoring mechanism",
    y = "Δ rmst"
  )




# 每个 trial(rep) 只保留一行
trial_level <- res %>%
  group_by(run_id) %>%
  slice(1) %>%          # or distinct(run_id, .keep_all = TRUE)
  ungroup()
dim(trial_level)


# reconstruction-level: three methods
method_level <- res



dgm_summary <- trial_level %>%
  group_by(censoring) %>%
  summarise(
    n_trials = n(),
    censor_prop_mean = mean(dgm_censor_prop, na.rm = TRUE),
    censor_prop_sd   = sd(dgm_censor_prop, na.rm = TRUE),
    events_ctrl_mean = mean(dgm_total_events_ctrl, na.rm = TRUE),
    events_trt_mean  = mean(dgm_total_events_trt, na.rm = TRUE),
    tau_mean         = mean(tau, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(censoring)

dgm_summary



recon_summary <- method_level %>%
  group_by(censoring, method) %>%
  summarise(
    n = n(),
    err_rate = mean(!is.na(reconstruct_error)),
    km_rmse_mean = mean(km_rmse_mean, na.rm = TRUE),
    km_rmse_p95  = quantile(km_rmse_mean, 0.95, na.rm = TRUE),
    # Cox: report both signed diff and absolute diff
    cox_diff_mean = mean(cox_logHR_diff, na.rm = TRUE),
    cox_diff_sd   = sd(cox_logHR_diff, na.rm = TRUE),
    cox_abs_mean  = mean(abs(cox_logHR_diff), na.rm = TRUE),
    cox_ratio_mean = mean(exp(cox_logHR_diff), na.rm = TRUE),  # approx HR_rec / HR_true
    # Monte Carlo standard error for the mean (how stable your mean estimate is)
    cox_diff_mcse = cox_diff_sd / sqrt(sum(is.finite(cox_logHR_diff))),
    .groups = "drop"
  ) %>%
  arrange(censoring, method)

recon_summary



paired <- method_level %>%
  select(run_id, censoring, method, km_rmse_mean, cox_logHR_diff) %>%
  mutate(cox_abs = abs(cox_logHR_diff)) %>%
  pivot_wider(names_from = method,
              values_from = c(km_rmse_mean, cox_abs))

# 以 IPDfromKM 为基准，看看另外两种方法相对提升/变差
paired_delta <- paired %>%
  group_by(censoring) %>%
  summarise(
    n_trials = n(),
    d_km_rmse_kmdata  = mean(km_rmse_mean_kmdata  - km_rmse_mean_IPDfromKM, na.rm = TRUE),
    d_km_rmse_KMtoIPD = mean(km_rmse_mean_KMtoIPD - km_rmse_mean_IPDfromKM, na.rm = TRUE),
    p_kmdata_better   = mean((km_rmse_mean_kmdata  - km_rmse_mean_IPDfromKM) < 0, na.rm = TRUE),
    p_kmtoipd_better  = mean((km_rmse_mean_KMtoIPD - km_rmse_mean_IPDfromKM) < 0, na.rm = TRUE),

    d_coxabs_kmdata   = mean(cox_abs_kmdata  - cox_abs_IPDfromKM, na.rm = TRUE),
    d_coxabs_KMtoIPD  = mean(cox_abs_KMtoIPD - cox_abs_IPDfromKM, na.rm = TRUE),
    p_kmdata_cox_better  = mean((cox_abs_kmdata  - cox_abs_IPDfromKM) < 0, na.rm = TRUE),
    p_kmtoipd_cox_better = mean((cox_abs_KMtoIPD - cox_abs_IPDfromKM) < 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(censoring)

paired_delta
