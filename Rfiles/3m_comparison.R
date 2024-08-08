# Here is the function to output the result of 3 methods comparison
# sample size n (per arm), interim period, PET0 and overall alpha should be set in advanced
# PET1 as small as possible and compare the overall power.

# We take Minimax design in Jung 2017 table 1 as the standard


m3_compare <- function(n, sim_size, acc_time, cen_time, interim, lambda_H0, lambda_H1, H1_type, 
                    HR1, HR2, search_times, tar_a1, tar_pow1_low, tar_alpha) 
{

data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                            cen_time = cen_time, arm = 0, interim = interim)    
data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                            cen_time = cen_time, arm = 1, interim = interim) 

if (H1_type == 'PH'){
    data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H1, dist = 'exp', 
                            cen_time = cen_time, arm = 1, interim = interim)
}
else if (H1_type == 'NPH'){
    data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, 
                          dist = 'pcw_exp', cen_time = cen_time,HR1 = HR1, HR2 = HR2, 
                          change_time = change_time, arm = 1, interim = interim)
}

rmst_h0_int <- RMST_sim_cal(n = n,data_E = data_E_H0[ , c(2,3,1)], 
                                data_C = data_C[ , c(2,3,1)],tau = interim, sim_size = sim_size)
rmst_h1_int <- RMST_sim_cal(n = n,data_E = data_E_H1[ , c(2,3,1)], 
                                data_C = data_C[ , c(2,3,1)],tau = interim, sim_size = sim_size)
rmst_h0_fin <- RMST_sim_cal(n = n,data_E = data_E_H0[ , c(4,5,1)], 
                                data_C = data_C[ , c(4,5,1)],tau = acc_time + cen_time,sim_size = sim_size)
rmst_h1_fin <- RMST_sim_cal(n = n,data_E = data_E_H1[ , c(4,5,1)], 
                                data_C = data_C[ , c(4,5,1)],tau = acc_time + cen_time,sim_size = sim_size)
rmst_data <- rbind(rmst_h0_int, rmst_h1_int, rmst_h0_fin, rmst_h1_fin)
#Log rank
z_stats_h0_int <- log_rank_sim(data_C = data_C[ , c(2,3,1)], 
                                data_E = data_E_H0[ , c(2,3,1)], sim_size =  sim_size,
                            n = n, alpha = 0.05, sided = 'greater')$z_stats
z_stats_h1_int <- log_rank_sim(data_C = data_C[ , c(2,3,1)], 
                                data_E = data_E_H1[ , c(2,3,1)], sim_size =  sim_size,
                            n = n, alpha = 0.05, sided = 'greater')$z_stats
z_stats_h0_fin <- log_rank_sim(data_C = data_C[ , c(4,5,1)], 
                                data_E = data_E_H0[ , c(4,5,1)], sim_size =  sim_size,
                            n = n, alpha = 0.05, sided = 'greater')$z_stats
z_stats_h1_fin <- log_rank_sim(data_C = data_C[ , c(4,5,1)], 
                                data_E = data_E_H1[ , c(4,5,1)], sim_size =  sim_size,
                            n = n, alpha = 0.05, sided = 'greater')$z_stats
logrank_data <- rbind(z_stats_h0_int, z_stats_h1_int, z_stats_h0_fin, z_stats_h1_fin)

#Grid search critical value---------------------------
#RMST
m_low <- quantile((rmst_data[2,] - rmst_data[1,]), 0.1)
# Smallest RMST difference: interim under H0
m_up <- quantile((rmst_data[8,] - rmst_data[7,]), 0.9)
# Largest RMST difference: final under H1
search_step_rmst <- (m_up - m_low) / search_times
t_low <- quantile((rmst_data[2,]), 0.1)
# Smallest experiment RMST: interim under H0
t_up <- quantile((rmst_data[8,]), 0.9)
#Log rank
c_low <- quantile(logrank_data, 0.2)
c_up <- quantile(logrank_data, 0.8)
search_step_lr <- (c_up - c_low) / search_times

# Control PET0 
best_our_rmst <- find_m_t_RMST(m_low = m_low, t_low = t_low, t_up = t_up, rmst_data = rmst_data, 
                            search_times = search_times, search_step = search_step_rmst, tar_a1 = tar_a1, 
                            tar_pow1_low = tar_pow1_low, tar_a2 = tar_alpha, sim_size = sim_size)

best_simple_rmst <- find_m_t_RMST(m_low = m_low, t_low = -Inf, t_up = t_up, rmst_data = rmst_data, 
                            search_times = search_times, search_step = search_step_rmst, tar_a1 = tar_a1, 
                            tar_pow1_low = tar_pow1_low, tar_a2 = tar_alpha, sim_size = sim_size)

best_log_rank <- find_m_logrank(m_low = c_low, logrank_data = logrank_data, 
                            search_times = search_times, search_step = search_step_lr, tar_a1 = tar_a1, 
                            tar_pow1_low = tar_pow1_low, tar_a2 = tar_alpha, sim_size = sim_size)

return(list(best_our_rmst = best_our_rmst, 
            best_simple_rmst = best_simple_rmst,
            best_log_rank = best_log_rank))

}
