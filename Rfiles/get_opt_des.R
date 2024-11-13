#*********************************
# Here is the function to get optimal design. 
# It can return the interim sample size n and other measurement of the design with min(E(N))
# The overall sample size N is fixed.
#_________________________________________________________________
# Parameter 'method' : 'Simple' means simple RMST difference test, 'Complex' means our method.
# Parameter 'int_step' : difference of interim sample size between two searching step.
#                        e.g int_step = 4 means increasing int_n by 4(2 for each group) each time
# Parameter 'n': To maintaining consistency the input n should be the overall sample size of each group.
# Parameter 'logrank': If logrank = TRUE, search for the optimal n for log_rank test
#_________________________________________________________________

get_opt_des <- function(n, sim_size, acc_time, cen_time, int_step, method, lambda_H0, 
                        lambda_H1, H1_type, HR1, HR2, change_time, alpha, power) 
{
    N <- 2 * n #overall sample size of two groups
    r <- N / acc_time
    int_factor <- seq(0.4, 0.7, by = int_step / N)  # Each time interim sample size increase by 6
    interim_list <- int_factor * acc_time

    data_C <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                            cen_time = cen_time, arm = 0, interim = interim_list)    
    data_E_H0 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, dist = 'exp', 
                            cen_time = cen_time, arm = 1, interim = interim_list) 
    if (H1_type == 'PH'){
    data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H1, dist = 'exp', 
                            cen_time = cen_time, arm = 1, interim = interim_list)
    }
    else if (H1_type == 'NPH'){
    data_E_H1 <- expo_gen_2stages(N = n * sim_size, acc_time = acc_time, lambda = lambda_H0, 
                            dist = 'pcw_exp', cen_time = cen_time, HR1 = HR1, HR2 = HR2, 
                            change_time = change_time, arm = 1, interim = interim_list)
    }

    all_result <- data.frame()
    for (i in 1 : length(interim_list))
    {  
        interim <- interim_list[i]  # interim is the length of interim period
        if (method != 'logrank')
        {
            rmst_h0_int <- RMST_sim_cal(n = n,data_E = data_E_H0[ , c(2,3,1), i], 
                                data_C = data_C[ , c(2,3,1), i],tau = interim, sim_size = sim_size)
            rmst_h1_int <- RMST_sim_cal(n = n,data_E = data_E_H1[ , c(2,3,1), i], 
                                data_C = data_C[ , c(2,3,1), i],tau = interim, sim_size = sim_size)
            rmst_h0_fin <- RMST_sim_cal(n = n,data_E = data_E_H0[ , c(4,5,1), i], 
                                data_C = data_C[ , c(4,5,1), i],
                                tau = acc_time + cen_time, sim_size = sim_size)
            rmst_h1_fin <- RMST_sim_cal(n = n,data_E = data_E_H1[ , c(4,5,1), i], 
                                data_C = data_C[ , c(4,5,1), i],
                                tau = acc_time + cen_time, sim_size = sim_size)
            rmst_data <- rbind(rmst_h0_int, rmst_h1_int, rmst_h0_fin, rmst_h1_fin)

            mu_cov_h0 <- mu_cov_mc(rmst_int = rmst_h0_int, rmst_fin = rmst_h0_fin, sim_size = sim_size)
            mu_cov_h1 <- mu_cov_mc(rmst_int = rmst_h1_int, rmst_fin = rmst_h1_fin, sim_size = sim_size)

            best_our <- adp_grid_src(rmst_data = rmst_data, mu_cov_h0 = mu_cov_h0, mu_cov_h1 = mu_cov_h1, 
                    int_n = interim * r, fin_n = N, sim_size = sim_size, method = method,
                    alpha = alpha, power = power)
        }

        else if (method == 'logrank')   # search the min(E(N)) using log rank test
        {
            lr_h0_int <- log_rank_sim(data_C = data_C[ , c(2,3,1), i], data_E = data_E_H0[ , c(2,3,1), i], 
                        sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
            lr_h1_int <- log_rank_sim(data_C = data_C[ , c(2,3,1), i], data_E = data_E_H1[ , c(2,3,1), i], 
                        sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
            lr_h0_fin <- log_rank_sim(data_C = data_C[ , c(4,5,1), i], data_E = data_E_H0[ , c(4,5,1), i], 
                        sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
            lr_h1_fin <- log_rank_sim(data_C = data_C[ , c(4,5,1), i], data_E = data_E_H1[ , c(4,5,1), i], 
                        sim_size =  sim_size, n = n, alpha = alpha, sided = 'greater')
            
            # Get W/sigma
            z_stats_h1_int <- lr_h1_int$z_stats
            z_stats_h1_fin <- lr_h1_fin$z_stats
            z_stats_h0_int <- lr_h0_int$z_stats
            z_stats_h0_fin <- lr_h0_fin$z_stats
            logrank_data <- rbind(z_stats_h0_int, z_stats_h1_int, z_stats_h0_fin, z_stats_h1_fin) 

            # corr(W1, W | H0)
            corr_h0 <- sqrt(mean(lr_h0_int$var_w) / mean(lr_h0_fin$var_w))         
            best_our <- find_m_logrank(logrank_data = logrank_data, sim_size = sim_size, corr_h0 = corr_h0,
                            search_times = 100, alpha = alpha, power = power, int_n = interim * r, fin_n = N)
        }

        best_our$interim_n <- ceiling(interim * r)
        all_result <- rbind(all_result, best_our)
    }

    all_result <- na.omit(all_result)    #valid result 
    if (dim(all_result)[1] == 0) {
        if(method == 'logrank'){
            return(data.frame(m1 = 0, m2 = 0, PET0 = 0, PET1 = 0, alpha = 0, 
                        power = 0, EN0 = NA, EN1 = NA, EN = NA, interim_n = NA))
        }
        else{
             return(data.frame(m1 = 0, t1 = 0, m2 = 0, t2 = 0, lambda = 0,
                        gamma = 0, PET0 = 0, PET1 = 0, alpha = 0, 
                        power = 0, EN0 = NA, EN1 = NA, EN = NA, interim_n = NA))
        }
    }
    else{
        return(all_result[which(all_result$EN == min(all_result$EN, na.rm = TRUE)), ])
        # return the result with minimal EN
    }
   
}