source("R_meta/meta_sim_v3.R")


grid <- build_scenario_grid(
  p_miss_vec = c(0, 0.2, 0.4, 0.6, 0.8),
  tau_levels = list(low = 0.05, middle = 0.1, high = 0.25),
  n_levels   = list(standard = 250),
  K = 30
)

sim <- run_meta_simulation(grid, n_rep = 300)
sim2 <- run_meta_simulation_paired(grid, n_rep=300, digitize_sd_S = 0.04)
sim3 <- run_meta_simulation_paired(grid, n_rep=300, digitize_sd_S = 0.04)
str(sim3)

write.csv(sim2$res_methods_long, file = "R_meta/meta_sim_k30_n300_outputs_sd004_more_points/sim_reps300_results_long.csv")
write.csv(sim2$res_studies_long, file = "R_meta/meta_sim_k30_n300_outputs_sd004_more_points/sim_reps300_studies_long.csv")

write.csv(sim3$res_methods_long, file = "R_meta/meta_sim_k30_n300_outputs_sd004/sim_reps300_results_long.csv")
write.csv(sim3$res_studies_long, file = "R_meta/meta_sim_k30_n300_outputs_sd004/sim_reps300_studies_long.csv")

perf <- summarise_performance(sim$res_methods_long, mu_HR = 0.75)
plot_bias_rmse(perf, out_dir = "R_meta/meta_sim_k30_n100_outputs")

perf2 <- summarise_performance(sim2$res_methods_long, mu_HR=0.75)
plot_bias_rmse(perf2, out_dir="R_meta")

perf3 <- summarise_performance(sim3$res_methods_long, mu_HR=0.75)
plot_bias_rmse(perf3, out_dir="R_meta/meta_sim_k30_n300_outputs_sd004")

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

mu_HR <- 0.75
mu_logHR <- log(mu_HR)

meth_levels <- c("AD_only","Hybrid","RM","Pooled_IPD")
meth_labels <- c(
  AD_only    = "AD-only (drop missing HR)",
  Hybrid     = "Hybrid (HR + recon for missing)",
  RM         = "RM (all recon)",
  Pooled_IPD = "Pooled IPD (oracle)"
)



methods <- as.data.table(sim2$res_methods_long)
methods <- methods[method %in% meth_levels]
methods[, method := factor(method, levels = meth_levels, labels = meth_labels[meth_levels])]

# 过滤掉失败的rep（如果你想保守）
methods_ok <- methods[ok == TRUE & is.finite(logHR_hat) & is.finite(se_hat)]

perf2 <- methods_ok[, .(
  n_rep = .N,
  mean_logHR = mean(logHR_hat),
  bias_logHR = mean(logHR_hat - mu_logHR),
  abs_bias_logHR = mean(abs(logHR_hat - mu_logHR)),
  rmse_logHR = sqrt(mean((logHR_hat - mu_logHR)^2)),
  mean_se = mean(se_hat),
  mean_tau2_hat = mean(tau2_hat, na.rm=TRUE),
  mean_I2_hat = mean(I2_hat, na.rm=TRUE),
  mean_k_used = mean(k_used)
), by = .(scenario_id, p_miss, tau_label, tau_logHR, n_label, n_per_arm, K, method)]



# bias 的第二种画法
p_bias <- ggplot(perf2, aes(x = p_miss, y = bias_logHR,
                            color = method,    # 新增：按 method 分配颜色
                            group = method)) +  # 保持分组
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~tau_label, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Missing HR proportion",
    y = "Bias in log(HR) (estimate - log(0.75))",
    color = "Method",      # 将图例标题改为 "Method"
    title = "Meta-level bias vs missingness"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
p_bias

# 检查bias的汇总表
perf2 %>%
  select(p_miss, tau_label, method, bias_logHR) %>%
  filter(tau_label == "low") %>%
  arrange(tau_label, method, p_miss)


# RMSE 的第二种画法
p_rmse <- ggplot(perf2, aes(x = p_miss, y = rmse_logHR, group = method, color = method)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~tau_label, nrow = 1) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9), labels = percent_format(accuracy = 1)) +
  labs(
    x = "Missing HR proportion",
    y = "RMSE of log(HR)",
    color = "Method",
    title = "Meta-level RMSE vs missingness"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
p_rmse

# 检查rmse的汇总表
perf2 %>%
  select(p_miss, tau_label, method, rmse_logHR) %>%
  filter(tau_label == "low") %>%
  filter(method == "RM (all recon)") %>%
  arrange(tau_label, method, p_miss)


p_kused <- ggplot(methods_ok, aes(x = factor(p_miss), y = k_used, color = method)) +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_grid(tau_label ~ method) +
  labs(
    x = "Missing HR proportion",
    y = "k used in meta",
    title = "How many studies each method uses",
    color = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "right")

p_kused

# 计算 k_used 的统计摘要（按 p_miss、tau_label、method 分组）
methods_ok %>%
  group_by(p_miss, tau_label, method) %>%
  summarise(
    min = min(k_used, na.rm = TRUE),
    q1 = quantile(k_used, 0.25, na.rm = TRUE),
    median = median(k_used, na.rm = TRUE),
    q3 = quantile(k_used, 0.75, na.rm = TRUE),
    max = max(k_used, na.rm = TRUE),
    mean = mean(k_used, na.rm = TRUE),
    sd = sd(k_used, na.rm = TRUE),
    n = n(),  # 每个组合的观测数
    .groups = "drop"
  ) %>%
  arrange(tau_label, method, p_miss)




ad <- methods_ok[method == meth_labels["AD_only"]]

p_ad_scatter <- ggplot(ad, aes(x = k_used, y = abs(logHR_hat - mu_logHR))) +
  geom_point(alpha = 0.35) +
  facet_wrap(~tau_label) +
  labs(
    x = "k_used (AD-only)",
    y = "Absolute error |logHR_hat - log(0.75)|",
    title = "AD-only: error vs number of included studies"
  ) +
  theme_bw()

p_ad_scatter

# 数值相关性（每个tau、每个p_miss分别看）
ad_cor <- ad[, .(
  cor_abs_err_k = cor(abs(logHR_hat - mu_logHR), k_used, use="complete.obs"),
  mean_k = mean(k_used),
  mean_abs_err = mean(abs(logHR_hat - mu_logHR))
), by = .(tau_label, p_miss)]
ad_cor



studies <- as.data.table(sim3$res_studies_long)

studies_ok <- studies[true_ok == TRUE & recon_ok == TRUE &
                        is.finite(logHR_true_hat) & is.finite(logHR_recon_hat)]

studies_ok[, recon_err := logHR_recon_hat - logHR_true_hat]

# 总体分布（按tau）
p_recon_err <- ggplot(studies_ok, aes(x = recon_err)) +
  geom_histogram(bins = 50) +
  facet_wrap(~tau_label, nrow = 1) +
  labs(
    x = "logHR(recon) - logHR(true)",
    y = "Count",
    title = "Trial-level reconstruction error distribution"
  ) +
  theme_bw()

p_recon_err



p_recon_err_box <- ggplot(studies_ok, aes(x = factor(p_miss), y = recon_err)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~tau_label, nrow = 1) +
  labs(
    x = "Missing HR proportion (scenario label)",
    y = "logHR(recon) - logHR(true)",
    title = "Reconstruction error by scenario"
  ) +
  theme_bw()

p_recon_err_box


p_scatter_true_recon <- ggplot(studies_ok, aes(x = logHR_true_hat, y = logHR_recon_hat)) +
  geom_point(alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.6) +
  facet_wrap(~tau_label, nrow = 1) +
  labs(
    x = "Trial logHR from TRUE IPD",
    y = "Trial logHR from RECON IPD",
    title = "Trial-level: reconstructed HR tracks true HR"
  ) +
  theme_bw()

p_scatter_true_recon


recon_summary <- studies_ok[, .(
  n_trials = .N,
  mean_err = mean(recon_err),
  sd_err = sd(recon_err),
  rmse_err = sqrt(mean(recon_err^2)),
  cor_true_recon = cor(logHR_true_hat, logHR_recon_hat)
), by = .(tau_label, p_miss)]
recon_summary


methods_ok[, ci_lo := logHR_hat - qnorm(0.975) * se_hat]
methods_ok[, ci_hi := logHR_hat + qnorm(0.975) * se_hat]
methods_ok[, cover := (ci_lo <= mu_logHR & mu_logHR <= ci_hi)]

cover_dt <- methods_ok[, .(
  coverage = mean(cover, na.rm = TRUE),
  mean_len = mean(ci_hi - ci_lo, na.rm = TRUE),
  n_rep = .N
), by = .(p_miss, tau_label, method)]

p_cov <- ggplot(cover_dt, aes(x = p_miss, y = coverage, group = method, color = method)) +
  geom_hline(yintercept = 0.95, linewidth = 0.4) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~tau_label, nrow = 1) +
  scale_x_continuous(breaks = c(0.3,0.6,0.9), labels = percent_format(accuracy=1)) +
  scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy=1)) +
  labs(
    x = "Missing HR proportion",
    y = "Coverage of 95% CI",
    color = "Method",
    title = "Coverage vs missingness"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

p_cov
