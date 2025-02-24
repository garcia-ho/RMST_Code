# Here is the function to output the result of 3 methods comparison
# sample size n (per arm), interim period, PET0 and overall alpha should be set in advanced
# PET1 as small as possible and compare the overall power.

# We take Minimax design in Jung 2017 table 1 as the standard

m3_compare <- function(n, sim_size, acc_time, cen_time, interim, lambda_H0, tau = NULL,
                        lambda_H1, H1_type, HR1, HR2, alpha, change_time) 
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
                            dist = 'pcw_exp', cen_time = cen_time, HR1 = HR1, HR2 = HR2, 
                            change_time = change_time, arm = 1, interim = interim)
}

if (is.null(tau)){ # the preset cut off tau for final stage
        tau_f <- acc_time + cen_time
    }
    else{
        tau_f <- tau
    }


result_list <-foreach(i = 1:length(interim), .combine = rbind) %do% {

    interim_i <- interim[i]

    rmst_h0_int <- RMST_sim_cal(n = n, data_E = data_E_H0[, c(2, 3, 1), i], 
                                data_C = data_C[, c(2, 3, 1), i], tau = interim_i, sim_size = sim_size)
    rmst_h1_int <- RMST_sim_cal(n = n, data_E = data_E_H1[, c(2, 3, 1), i], 
                                data_C = data_C[, c(2, 3, 1), i], tau = interim_i, sim_size = sim_size)
    rmst_h0_fin <- RMST_sim_cal(n = n, data_E = data_E_H0[, c(4, 5, 1), i], 
                                data_C = data_C[, c(4, 5, 1), i], tau = tau_f, sim_size = sim_size)
    rmst_h1_fin <- RMST_sim_cal(n = n, data_E = data_E_H1[, c(4, 5, 1), i], 
                                data_C = data_C[, c(4, 5, 1), i], tau = tau_f, sim_size = sim_size)
    rmst_data <- rbind(rmst_h0_int, rmst_h1_int, rmst_h0_fin, rmst_h1_fin)

    mu_cov_h0 <- mu_cov_mc(rmst_int = rmst_h0_int, rmst_fin = rmst_h0_fin, sim_size = sim_size)
    mu_cov_h1 <- mu_cov_mc(rmst_int = rmst_h1_int, rmst_fin = rmst_h1_fin, sim_size = sim_size)

    # log rank data
    try_lr <- tryCatch({
        lr_h0_int <- log_rank_sim(data_C = data_C[, c(2, 3, 1), i], 
                                  data_E = data_E_H0[, c(2, 3, 1), i], 
                                  sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
        lr_h1_int <- log_rank_sim(data_C = data_C[, c(2, 3, 1), i], 
                                  data_E = data_E_H1[, c(2, 3, 1), i], 
                                  sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
        lr_h0_fin <- log_rank_sim(data_C = data_C[, c(4, 5, 1), i], 
                                  data_E = data_E_H0[, c(4, 5, 1), i], 
                                  sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
        lr_h1_fin <- log_rank_sim(data_C = data_C[, c(4, 5, 1), i], 
                                  data_E = data_E_H1[, c(4, 5, 1), i], 
                                  sim_size = sim_size, n = n, alpha = alpha, sided = 'greater')
        list(lr_h0_int = lr_h0_int, lr_h1_int = lr_h1_int, lr_h0_fin = lr_h0_fin, lr_h1_fin = lr_h1_fin)
    }, error = function(e) {
        return(NULL)
    })

    if (is.null(try_lr)) {
        return(NULL)
    }
    lr_h0_int <- try_lr$lr_h0_int
    lr_h1_int <- try_lr$lr_h1_int
    lr_h0_fin <- try_lr$lr_h0_fin
    lr_h1_fin <- try_lr$lr_h1_fin

    if (is.null(lr_h0_int) || is.null(lr_h1_int) || is.null(lr_h0_fin) || is.null(lr_h1_fin)) {
        return(NULL)
    }
    # Get W/sigma
    z_stats_h1_int <- lr_h1_int$z_stats
    z_stats_h1_fin <- lr_h1_fin$z_stats
    z_stats_h0_int <- lr_h0_int$z_stats
    z_stats_h0_fin <- lr_h0_fin$z_stats
    logrank_data <- rbind(z_stats_h0_int, z_stats_h1_int, z_stats_h0_fin, z_stats_h1_fin)

    # corr(W1, W | H0)
    corr_h0 <- sqrt(mean(lr_h0_int$var_w) / mean(lr_h0_fin$var_w))

    best_our <- adp_grid_src(rmst_data = rmst_data, mu_cov_h0 = mu_cov_h0, mu_cov_h1 = mu_cov_h1, 
                             int_n = interim_i * 2 * n * acc_time, fin_n = 2 * n, alpha = alpha, 
                             sim_size = sim_size, method = 'Complex')
    best_lr <- find_m_logrank(logrank_data = logrank_data, search_times = 100, corr_h0 = corr_h0,
                              alpha = alpha, sim_size = sim_size)
    best_rmst <- adp_grid_src(rmst_data = rmst_data, mu_cov_h0 = mu_cov_h0, mu_cov_h1 = mu_cov_h1, 
                              int_n = interim_i * 2 * n * acc_time, fin_n = 2 * n, alpha = alpha, 
                              sim_size = sim_size, method = 'Simple')

    c(best_lr$alpha[1], best_rmst$alpha[1], best_our$alpha[1], best_lr$power[1], best_rmst$power[1], best_our$power[1],
      best_lr$PET0[1], best_rmst$PET0[1], best_our$PET0[1], best_lr$PET1[1], best_rmst$PET1[1], best_our$PET1[1])

    }
}
