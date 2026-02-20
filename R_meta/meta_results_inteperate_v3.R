# After running 'sim3 <- run_meta_simulation_paired(grid, n_rep=300, digitize_sd_S = 0.04)'
# I would like to analysis and interperate the results in 'sim3$res_methods_long' and 'sim3$res_studies_long'. I will focus on the performance of the different methods (AD-only, Hybrid, RM, Pooled IPD) in terms of bias, RMSE, and coverage probability. I will also look at how the digitization error (sd_S = 0.04) affects the results compared to a smaller digitization error (sd_S = 0.01).

library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)

# 将 data.table 转为 data.frame 方便 dplyr 操作（也可直接用 data.table）
methods_df <- as.data.frame(sim2$res_methods_long)

# 按场景和方法分组汇总
perf_summary <- methods_df %>%
  group_by(scenario_id, p_miss, tau_label, tau_logHR, n_per_arm, K, method) %>%
  summarise(
    n_rep = n(),                              # 重复次数
    fail_rate = mean(!ok, na.rm = TRUE),      # meta分析失败的比例（如收敛失败）
    mean_k_used = mean(k_used, na.rm = TRUE), # 平均使用试验数
    bias = mean(err, na.rm = TRUE),            # 偏差
    abs_bias = mean(abs(err), na.rm = TRUE),   # 绝对偏差
    rmse = sqrt(mean(err^2, na.rm = TRUE)),    # RMSE
    coverage = mean(cover, na.rm = TRUE),      # 覆盖概率
    .groups = "drop"
  )

# 查看汇总结果: 60 = 3(tau) x 5(p_miss) x 4(method)
print(perf_summary, n=60)

perf_summary %>%
  group_by(method) %>%
  summarise(mean_bias = mean(bias), mean_rmse = mean(rmse), mean_coverage = mean(coverage), mean_k_used = mean(mean_k_used))

perf_summary %>%
  group_by(method, p_miss) %>%
  summarise(mean_bias = mean(bias), mean_rmse = mean(rmse), mean_coverage = mean(coverage), mean_k_used = mean(mean_k_used))

perf_summary %>%
  group_by( tau_logHR) %>%
  summarise(mean_bias = mean(bias), mean_rmse = mean(rmse), mean_coverage = mean(coverage), mean_k_used = mean(mean_k_used))


# 重建成功的概率
studies_df <- as.data.frame(sim2$res_studies_long)

# 计算每个场景下，缺失试验的重建成功率
recon_success <- studies_df %>%
  filter(is_missing == TRUE) %>%               # 只考虑缺失的试验
  group_by(p_miss, tau_label, n_label, K, n_per_arm) %>%
  summarise(
    n_missing = n(),                            # 缺失试验总数
    n_recon_ok = sum(recon_ok, na.rm = TRUE),   # 重建成功的数量
    recon_success_rate = n_recon_ok / n_missing, # 成功率
    .groups = "drop"
  )

print(recon_success)




# 为了方便绘图，将 p_miss 转为数值因子，tau_label 作为分面
perf_plot <- perf_summary %>%
  filter(method != "Pooled_IPD") %>%   # 金标准作为参考线，可以不画或单独画
  mutate(p_miss = as.numeric(p_miss))

# 偏差图
ggplot(perf_plot, aes(x = p_miss, y = bias, color = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ tau_logHR, labeller = labeller(tau_label = c(low = "Low", middle = "Middle", high = "High heterogeneity"))) +
  labs(
    x = "Proportion of missing HR",
    y = "Bias (log HR scale)",
    color = "Method",
    title = "Bias of meta-analysis methods under different missing proportions"
  ) +
  theme_bw()

# RMSE图
ggplot(perf_plot, aes(x = p_miss, y = rmse, color = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ tau_logHR, labeller = labeller(tau_label = c(low = "Low", middle = "Middle", high = "High heterogeneity"))) +
  labs(
    x = "Proportion of missing HR",
    y = "RMSE (log HR scale)",
    color = "Method",
    title = "RMSE of meta-analysis methods"
  ) +
  theme_bw()

# 方差图
# 从 res_methods_long 重新计算每个场景各方法的平均标准误
perf_summary_var <- sim3$res_methods_long %>%
  as.data.frame() %>%
  group_by(scenario_id, p_miss, tau_label, n_label, method) %>%
  summarise(
    mean_se = mean(se_hat, na.rm = TRUE),
    mean_var = mean(se_hat^2, na.rm = TRUE),
    .groups = "drop"
  )
# 准备绘图数据（可包含 Pooled_IPD 作为参考）
plot_data_var <- perf_summary_var %>%
  mutate(p_miss = as.numeric(p_miss))


# 平均标准误图
ggplot(plot_data_var, aes(x = p_miss, y = mean_var, color = method, group = method)) +
  geom_line() +
  geom_point(size = 2) +
  facet_wrap(~ tau_label, labeller = labeller(tau_label = c(low="Low",  high="High",middle="Middle"))) +
  labs(
    x = "Proportion of missing HR",
    y = "Mean variance (log HR scale)",
    color = "Method",
    title = "Precision of meta-analysis methods (smaller SE = better precision)"
  ) +
  theme_bw()
