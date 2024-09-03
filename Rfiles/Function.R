#—————————————————————————————————————————————————————
# All functions used in this project are stored here!! 
#—————————————————————————————————————————————————————

#-------------- 1. expo_gen_2stages------------
# Generate exponential dist survival data for 2 stages test.
# The censoring distribution in interim period is different from the whole trail 
# It returns an (N * 3 * 2) array : obs_survival_int, event_int, obs_survival_fin, event_fin, arm  
# 3 columns are: obs_survival, event, arm  
# If interim is a c(), result[, , i] is the result of ith interim (N * 5)

  # N: Number of patients
  # acc_time: Accrual time period with constant rate
  # lambda: for exponential distribution
  # cen_time: extra censoring period after accrual period
  # arm: group label()

expo_gen_2stages <- function(N,dist,acc_time,cen_time,lambda,HR1,HR2,arm,interim,change_time)
{

if (dist == 'exp') {
  survival_time_all <- rexp(N, rate = lambda)
} else if (dist == 'pcw_exp') { #piecewise exponential with one change time point
  survival_time_all <- rpwexp(n = N,fail_rate = data.frame(rate = c(lambda * HR1, lambda * HR2), 
                                            duration = c(change_time)))
}

arm_all <- rep(arm, N)
entry_time_all <- runif(N, min = 0, max = acc_time)
censor_time_fin <- runif(N, min = cen_time, max = acc_time + cen_time)
all_data <- cbind(survival_time_all, entry_time_all, censor_time_fin, arm_all)

# _______________________Define function for event counting________________________________
# filter the sample die before interim
cal_event_int <- function(row) { 
  survival_time_all <- row[1]
  censor_time_fin <- row[3]
  entry_time_all <- row[2]
  if (entry_time_all >= interim_val) { # Not in the interim
    return(c(0,0))
  } else if (entry_time_all + min(survival_time_all, censor_time_fin) < interim_val) {
    return(c(min(survival_time_all, censor_time_fin), 1))
  } else if (entry_time_all + min(survival_time_all, censor_time_fin) > interim_val &&
           entry_time_all < interim_val) {
    return(c(interim_val - entry_time_all, 0))
  }
}
# calculate the event in final study
cal_event_fin <- function(row) { 
  survival_time_all <- row[1]
  entry_time_all <- row[2]
  censor_time_fin <- row[3]
  obs_time_int <- row[5]
  event_int <- row[6]
if(entry_time_all >= interim_val) {  # censoring in stage II
  return(c(min(survival_time_all,censor_time_fin), 
          as.integer(survival_time_all <= censor_time_fin)))
}
else if (event_int == 0 && obs_time_int != 0) { # censoring in interim
  return(c(min(survival_time_all,censor_time_fin), 
          as.integer(survival_time_all <= censor_time_fin)))
}
else {
  return(c(obs_time_int,event_int))
}
}
#_________________________________________________________________

if (length(interim) == 1) {
    interim_val <- interim
   obs_event_int <- apply(all_data, 1, cal_event_int)
   sur_data_int <- cbind(all_data,t(obs_event_int))
   obs_event_fin <- apply(sur_data_int, 1, cal_event_fin) 
   all_data <- cbind(sur_data_int,t(obs_event_fin)) 
   return(all_data[ ,c(4,5,6,7,8)])
# In all_data 4th column is arm, 5th obs_time_int, 6th event_int, 7th obs_time_fin, 8th event_fin
}

else 
  {  # interim is a c() it return different interim event 
    result <- array(NA, dim = c(N , 5, length(interim)))
    for (i in 1 : length(interim))
    {
      interim_val <- interim[i]
      obs_event_int <- apply(all_data, 1, cal_event_int)
      sur_data_int <- cbind(all_data,t(obs_event_int))
      obs_event_fin <- apply(sur_data_int, 1, cal_event_fin) 
      result[, , i] <- cbind(sur_data_int,t(obs_event_fin))[ ,c(4,5,6,7,8)]
    }
   return(result)
  }
}





#------ 2. RMST_sim_cal-------
# Generate Survival function and calculate the RMST of each arm.
# It returns 2 * sim_size matrix. First row control, second one experiment. sim_size is the times of simulation
# The out put of this function is array to accelerate the grid search code.

RMST_sim_cal <- function(n,data_E,data_C,tau,sim_size)
{
  # n is the sample size of each group
  # sim_size is the times of simulation
  sim_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'survRM2') %dopar% {
    pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
    # the simulation data is not guaranteed to be larger than tau
    if (tau < min(max(data_E[((k-1)*n+1):(k*n),1]),  max(data_C[((k-1)*n+1):(k*n),1]))) {  
        rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], tau = tau) # No need for adjustment
      }  else {
        rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3])
        # The rmst2 function will automatically adjust tau for us if it's not specified
      }
    c(rmst_result$RMST.arm0$rmst[1],rmst_result$RMST.arm1$rmst[1])
    }
    
    return(sim_result)
}





#---------------- 3. RMST_sim_test-------
# Different from RMST_sim_cal, It return a dataframe of two_sided RMST test rejection times over simulation times.
# It also counts the times of tau adjustment
# This function can be used to compare our rejection method with classical RMST diff test 
# It can return the test result over sim_size, tau adjustment proportion and p value of each test

RMST_sim_test <- function(n, data_E, data_C, tau, sim_size, alpha, sided)
{
  # n is the sample size of each group
  # sim_size is the times of simulation
    tau_adj_count <- 0
    sim_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'survRM2') %dopar% {
          test_res <- 0
          p <- 0
          pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
          
        if (tau < min(max(data_E[((k-1)*n+1):(k*n),1]),  max(data_C[((k-1)*n+1):(k*n),1]))) {
            rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], tau = tau, alpha = alpha)
            # the simulation data is not guaranteed to be larger than tau
          }  else {
            rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], alpha = alpha)
            tau_adj_count <- tau_adj_count + 1   # tau is adjusted automatically
          }
        
        if (sided == 'two_sided') {
          p <- rmst_result$unadjusted.result[1,4]  
          # The p value of RMST difference test. It's a two sided test in the package
          } else if (sided == 'greater') {
          diff <- rmst_result$unadjusted.result[1,1]
          std <- (rmst_result$unadjusted.result[1,1] - 
                  rmst_result$unadjusted.result[1,2]) / qnorm(1 - alpha/2)
          p <- 1 - pnorm(diff / std)
          }
        if ( p <= alpha ) {
            test_res <- 1
          } else {
            test_res <- 0
          }
        c(test_res,tau_adj_count, p )  # The last element of the 2nd row is tau_adj_count
        }

    return(list(test_result = data.frame('rejection' = sum(sim_result[1, ]) / sim_size, 
                                      'tau adjustment' = sim_result[2, sim_size] / sim_size),
                  p_value = sim_result[3, ])) 
}




#--------- 4. emp_est -------
# Estimate the empirical RMST and τ
emp_est <- function(acc_time,cen_time,lambda_H0,lambda_H1)
# estimate the empirical RMST under H1 (hazard ratio), and the tau
  {
  sur_0 <- expo_gen(N = 10000,acc_time = acc_time, lambda = lambda_H0, cen_time = cen_time,arm = 0)
  sur_1 <- expo_gen(N = 10000,acc_time = acc_time, lambda = lambda_H1, cen_time = cen_time,arm = 1)
  est_tau <- mean(sur_0[1,]) + 3 * sd(sur_0[1,])
  pre_data <- rbind(sur_0,sur_1)
  rmst_result <- rmst2(pre_data[,1], pre_data[,2], pre_data[,3], tau = est_tau)
  emp_rmst <- c(rmst_result$RMST.arm0$result['RMST','Est.'],rmst_result$RMST.arm1$result['RMST','Est.'])

  return (c(emp_rmst[1],emp_rmst[2],est_tau))
  }





#----------5. log_rank_sim --------
# For one-sided log rank test simulation. Return the simulated alpha

log_rank_sim <- function(data_C, data_E, sim_size, n, alpha, sided)
{
  # n is the sample size of each group
  # sim_size is the times of simulation
  logrank_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'nph') %dopar% {
    pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
    if (sided == 'greater') {
        result <- logrank.test(pre_data[,1],pre_data[,2],pre_data[,3],alternative = "greater")
      } 
    else if (sided == 'two_sided') {
        result <- logrank.test(pre_data[,1],pre_data[,2],pre_data[,3],alternative = "two.sided")
      }
    p <- result$test$p
    z_stats <- result$test$z  # return the z statistics for two stages trials
    c(p, z_stats)
  }
  return(list(rejection = sum(logrank_result[1, ] <= alpha) / sim_size, # rejection times
              z_stats = logrank_result[2, ]))  # the z statistics W/sigma of every simulation
}





#----------6. mdir_sim --------
# For one-sided multiple direction log rank test simulation. Return the simulated alpha

mdir_sim <- function(data_C, data_E, sim_size, n, alpha, iter)
{
  # This function only accept data.frame input and need correct colnames
  data_C <- as.data.frame(data_C)
  data_E <- as.data.frame(data_E)
  colnames(data_C) <- c('time', 'event', 'group')
  colnames(data_E) <- c('time', 'event', 'group')
  mdir_result <- foreach(k = 1:sim_size, .combine = 'cbind', .packages = 'mdir.logrank') %dopar% {
    pre_data <- rbind(data_C[((k-1)*n+1):(k*n),],data_E[((k-1)*n+1):(k*n),])
    result <- mdir.onesided(data = pre_data, group1 = 0, iter = iter, wild = "norm")
    p <- result$p_value
    test_stats <- result$stat # return the z statistics for two stages trials
    c(p, test_stats)
  }
  return(list(rejection = sum(mdir_result[1, ] <= alpha) / sim_size, # rejection times
              test_stats = mdir_result[2, ]))  # the z statistics W/sigma of every simulation
}





#--------------7. theo_RMST---------
# Calculate the theorecital RMST and RSDST of a explicit survival function

theo_RMST <- function(lambda, dist, tau) 
{
  if (dist == 'exp') 
    {
    surv_fun <- function(t) {
             exp(-lambda * t)
        }
    RMST <- integrate(surv_fun, lower = 0, upper = tau)$value
    
    }

    return(RMST)
}





#---------------8. PET_norm----------------
# It returns c(PET0, PET1)
# Input the RMST value and estimated variance of each group. and the critical value

PET_norm <- function(mu_c,var_c,mu_e,var_e,m1,t1)
{
    mu_h0 <- c(0, mu_c)
    sigma_h0 <- matrix(c(2*var_c, var_c, var_c, var_c), nrow = 2)
    mu_h1 <- c(mu_e-mu_c, mu_e)
    sigma_h1 <- matrix(c(var_e+var_c, var_e, var_e, var_c), nrow = 2)
    upper <- c(Inf, Inf)
    lower <- c(m1,t1)
    p_rj_h0 <- 1 - pmvnorm(lower, upper, mean = mu_h0, sigma = sigma_h0) # p(E-C>m1 & E>t1|H0)
    p_rj_h1 <- 1 - pmvnorm(lower, upper, mean = mu_h1, sigma = sigma_h1)
    
    return(c(p_rj_h0,p_rj_h1))
    
}





#---------------------- 9. mu_cov_mc -----------------------
# This function use Monte Carlo simulation to calculate the var_cov matrix of 
#   [E1-C1, E1, E2-C2, E2]
# The calculation formula refers to the Lu(2021) sequential trials.

#*************
# The input rmst_int should be the interim rmst data of two groups (RMST_sim_cal output)
# rmst_fin is the final rmst data of two groups. sim_size is the MC simulation times B
# true_rmst_int is the theoretical result of interim rmst of two groups c(RMST_C, RMST_E)

mu_cov_mc <- function(rmst_int, rmst_fin, sim_size){

    diff_C <- rbind(rmst_int[1, ] - mean(rmst_int[1, ]), rmst_fin[1, ] - mean(rmst_fin[1, ]))
    cov_C <- matrix(0, nrow = 2, ncol = 2)
    for (i in 1:sim_size) 
      {
        product <- diff_C[, i] %*% t(diff_C[, i])  
        cov_C <- cov_C + product  
      }
    cov_C <- cov_C / sim_size   # [ Var(C1)  Cov(C1, C2)
                                #  Cov(C1, C2)  Var(C2) ]

    diff_E <- rbind(rmst_int[2, ] - mean(rmst_int[2, ]), rmst_fin[2, ] - mean(rmst_fin[2, ]))
    cov_E <- matrix(0, nrow = 2, ncol = 2)
    for (i in 1:sim_size) 
      {
        product <- diff_E[,i] %*% t(diff_E[,i])  
        cov_E <- cov_E + product  
      }
    cov_E <- cov_E / sim_size   # [ Var(E1)  Cov(E1, E2)
                                #  Cov(E1, E2)  Var(E2) ]

    var_cov_all <- matrix(
            c(cov_E[1,1]+cov_C[1,1], cov_E[1,1], cov_C[1,2]+cov_E[1,2], cov_E[1,2],
              cov_E[1,1], cov_E[1,1], cov_E[1,2], cov_E[1,2],
              cov_C[1,2]+cov_E[1,2], cov_E[1,2], cov_E[2,2]+cov_C[2,2], cov_E[2,2],
              cov_E[1,2], cov_E[1,2], cov_E[2,2], cov_E[2,2]), nrow = 4, ncol = 4)

    mu <- c(mean(rmst_int[2,] - rmst_int[1,]), mean(rmst_int[2,]),
            mean(rmst_fin[2,] - rmst_fin[1,]), mean(rmst_fin[2,]))

    return(list(mu = mu,
                sigma = var_cov_all))
  }







#------------------ 10. find_m_t_RMST ------------------
# Grid searching loop Works for simple RMST difference test (E1 - C1 > m1 & E2 - C2 > m2)
# Given rmst data of interim and all, find the best m1, m2 for 
# n need to be given in this loop


find_m_t_RMST <- function(rmst_data, search_times, alpha, sim_size) {
  rmst_h0_int <- rmst_data[c(1,2) , ]
  rmst_h1_int <- rmst_data[c(3,4) , ]
  rmst_h0_fin <- rmst_data[c(5,6) , ]
  rmst_h1_fin <- rmst_data[c(7,8) , ]

  c_low <- 0  #superiority test
  c_up <- quantile((rmst_data[8,] - rmst_data[7,]), 0.9)

  result_m1 <- foreach(m1 = seq(from = c_low, to = c_up, by = (c_up - c_low) / search_times), 
                      .combine = 'cbind') %dopar% { 
      best_m2 <- c()
      opt_power <- 0
      for( m2 in seq(from = c_low, to = c_up, by = (c_up - c_low) / search_times))
        {
        proc_h0 <- sum((rmst_h0_int[2, ] - rmst_h0_int[1, ] > m1) & 
                       (rmst_h0_fin[2, ] - rmst_h0_fin[1, ] > m2))
        proc_h1 <- sum((rmst_h1_int[2, ] - rmst_h1_int[1, ] > m1) & 
                       (rmst_h1_fin[2, ] - rmst_h1_fin[1, ] > m2))
        if (proc_h0/sim_size > 0 & proc_h0/sim_size < alpha & proc_h1/sim_size >= opt_power) 
          {
            opt_power <- proc_h1/sim_size
            PET0 <- sum((rmst_h0_int[2, ] - rmst_h0_int[1, ] <= m1)) / sim_size
            PET1 <- sum((rmst_h1_int[2, ] - rmst_h1_int[1, ]  <= m1)) / sim_size
            best_m2 <- c(m1, m2, PET0, PET1, proc_h0/sim_size, opt_power)
          }
        }
      return(best_m2)
      }
  powerful_m1 <- result_m1[, which(result_m1[6, ] == max(result_m1[6, ]))]

  if (is.null(powerful_m1)) 
    {   # Return NULL when something goes wrong
      return(NULL)
    }

  if(is.null(dim(powerful_m1)))
    {
      opt_pet0_m1 <- powerful_m1
    }
  else
    {  # find the largest PET0 if multiply solution exist
      opt_pet0_m1 <- powerful_m1[, which(powerful_m1[3, ] == max(powerful_m1[3, ]))]
    }
  
    opt_pet0_m1 <- data.frame(t(opt_pet0_m1))
    colnames(opt_pet0_m1) <- c('m1', 'm2', 'PET0', 'PET1', 'alpha', 'power')
    return(opt_pet0_m1)
}






#------------------11. find_m_logrank------------------
# Similar to find_m_t_RMST, this one is for log rank test 2 stages design
# The z statistics (W/sigma) in logrank test of each simulation is required as an input
# The z statistics can be obtained by log_rank_sim function $ z_stats
# This function help you find the Z1 > m1 & Z > m2 to control the overall power.
# You need to tune tar_a1 and tar_pow1_low to control PET0 and PET1

find_m_logrank <- function( logrank_data, search_times, alpha, sim_size) 
{
  z_stats_h0_int <-  logrank_data[1, ]
  z_stats_h1_int <-  logrank_data[2, ]
  z_stats_h0_fin <-  logrank_data[3, ]
  z_stats_h1_fin <-  logrank_data[4, ]

  c_low <- quantile(logrank_data, 0.1)
  c_up <- quantile(logrank_data, 0.9) 

  result_m1 <- foreach(m1 = seq(from = c_low, to = c_up, by = (c_up - c_low) / search_times), 
                      .combine = 'cbind') %dopar% { 
      best_m2 <- c()
      opt_power <- 0
      for( m2 in seq(from = c_low, to = c_up, by = (c_up - c_low) / search_times))
        {
        proc_h0 <- sum((z_stats_h0_int > m1) & (z_stats_h0_fin > m2))
        proc_h1 <- sum((z_stats_h1_int > m1) & (z_stats_h1_fin  > m2))
        if (proc_h0/sim_size > 0 & proc_h0/sim_size < alpha & proc_h1/sim_size >= opt_power) 
          {
            opt_power <- proc_h1/sim_size
            PET0 <- sum((z_stats_h0_int <= m1)) / sim_size
            PET1 <- sum((z_stats_h1_int <= m1)) / sim_size
            best_m2 <- c(m1, m2, PET0, PET1, proc_h0/sim_size, opt_power)
          }
        }
      return(best_m2)
      }
    powerful_m1 <- result_m1[, which(result_m1[6, ] == max(result_m1[6, ]))]
    if (is.null(powerful_m1)) 
      {   # Return NULL when something goes wrong
        return(NULL)
      }

    if(is.null(dim(powerful_m1)))
      {
        opt_pet0_m1 <- powerful_m1
      }
    else
      {  # find the largest PET0 if multiply solution exist
        opt_pet0_m1 <- powerful_m1[, which(powerful_m1[3, ] == max(powerful_m1[3, ]))]
      }
      opt_pet0_m1 <- data.frame(t(opt_pet0_m1))
      colnames(opt_pet0_m1) <- c('m1', 'm2', 'PET0', 'PET1', 'alpha', 'power')
      return(opt_pet0_m1)
}




#-----------------------12. adp_grid_src----------------------
# This is the function for adaptive grid search of parameters lambda and gamma.
# In order to solve the critical values m & t, we made extra assumption based on exponential probability.
# Input mu_cov_h0 is the mu and cov of (E1-C1, E1, E2-C2, E2) under H0, output of func: mu_cov_mc
# interim sample size n.
# overall sample size N. 
# stated alpha.

adp_grid_src <- function(rmst_data, mu_cov_h0, mu_cov_h1, int_n, fin_n, alpha, sim_size) 
  {
      # Interim
      mu1 <- mu_cov_h1$mu[c(1,2)]
      sigma1 <- mu_cov_h1$sigma[1:2, 1:2]
      # Final
      mu2 <- mu_cov_h1$mu[c(3,4)]
      sigma2 <- mu_cov_h1$sigma[3:4, 3:4]

      # Function to minimize solve t in P(E-C > m & E > t) = tar_prob given m
      norm_2d <- function(t, m, mean, sigma, tar_prob) 
        {
          prob <- pmvnorm(lower = c(m, t), 
                          upper = rep(Inf, 2), 
                          mean = mean, 
                          sigma = sigma)
          return (prob - tar_prob)
        }

      rmst_h0_int <- rmst_data[c(1,2) , ]
      rmst_h1_int <- rmst_data[c(3,4) , ]
      rmst_h0_fin <- rmst_data[c(5,6) , ]
      rmst_h1_fin <- rmst_data[c(7,8) , ]
      #Grid search
    crit_val_res <- foreach(lambda = seq(0.01, 0.99, 0.01), .combine = 'cbind') %dopar%
      {   
        best_gamma <- c()
        best_power <- 0
        for (gamma in seq(0, 1, by = 0.01))
          {
            p1_tar <- exp(-gamma * (int_n / fin_n))            # P(E1-C1 > m1)
            p2_tar <- lambda * exp(-gamma * (int_n / fin_n))   # P(E1-C1 > m1, E1 > t1)
            p3_tar <- exp(-gamma * ( fin_n / fin_n))            # P(E2-C2 > m2)
            p4_tar <- lambda * exp(-gamma * (fin_n / fin_n))   # P(E2-C2 > m2, E2 > t2)

            # First equation P(E1-C1 > m1) = p1_tar
            m1 <- qnorm(1 - p1_tar, mean = mu1[1], sd = sqrt(sigma1[1, 1]))
            # Second equation P(E1-C1 > m1 & E1 > t1) = p2_tar
            t1 <- uniroot(norm_2d, interval = c(0, 100), m = m1, 
                          mean = mu1, sigma = sigma1, tar_prob = p2_tar)$root
            # Third equation
            m2 <- qnorm(1 - p3_tar, mean = mu2[1], sd = sqrt(sigma2[1, 1])) 
            # Forth equation
            t2 <- uniroot(norm_2d, interval = c(0, 100), m = m2, 
                        mean = mu2, sigma = sigma2, tar_prob = p4_tar)$root

            proc_h0 <- sum((rmst_h0_int[2, ] - rmst_h0_int[1, ] > m1) & (rmst_h0_int[2, ] > t1) &
                      (rmst_h0_fin[2, ] - rmst_h0_fin[1, ] > m2) & (rmst_h0_fin[2, ] > t2))
            proc_h1 <- sum((rmst_h1_int[2, ] - rmst_h1_int[1, ] > m1) & (rmst_h1_int[2, ] > t1) &
                      (rmst_h1_fin[2, ] - rmst_h1_fin[1, ] > m2) & (rmst_h1_fin[2, ] > t2))

            if (m1 > 0 & m2 > 0 & proc_h0 / sim_size <= alpha 
                & proc_h1 / sim_size > best_power)  #control alpha, find the most powerful set
              {
                best_power <- proc_h1 / sim_size
                best_gamma <- c(m1, t1, m2, t2, lambda, gamma, proc_h0/sim_size, proc_h1/sim_size)
              }
          }
        best_gamma 
        }

      if (is.null(crit_val_res)) 
      {   # Return NULL when something goes wrong
        return(NULL)
      }

      best_res <- crit_val_res[, which(crit_val_res[8, ] == max(crit_val_res[8, ]))]
      threshold <- best_res[1:4] # critical values
      PET0 <- sum((rmst_h0_int[2, ] - rmst_h0_int[1, ] < best_res[1]) | 
                  (rmst_h0_int[2, ] < best_res[2])) / sim_size
      PET1 <- sum((rmst_h1_int[2, ] - rmst_h1_int[1, ] < best_res[1]) | 
                  (rmst_h1_int[2, ] < best_res[2])) / sim_size

      return(data.frame(m1 = threshold[1],
                        t1 = threshold[2],
                        m2 = threshold[3],
                        t2 = threshold[4],
                        lambda = best_res[5],
                        gamma = best_res[6],
                        PET0 = PET0,
                        PET1 = PET1,
                        alpha = best_res[7],
                        power = best_res[8]))

    }



#-----------------------13. compare_line_plot----------------------
# Used to draw line plot for 3m_comparison under different scenario
# Input is the dataframe of 3m_comparison output. ** The order of variable is fixed **


compare_line_plot <- function(data, var_name) 
  {
    options(repr.plot.width = 20, repr.plot.height = 8)

    color_palette <- c("Our_power" = "darkred", "Our_alpha" = "darkred", 
                      "LR_power" = "lightgreen", "LR_alpha" = "lightgreen",
                      "Rmst_power" = "blue", "Rmst_alpha" = "blue",
                      "Our_PET0" = "darkred", "Our_PET1" = "darkred", 
                      "LR_PET0" = "lightgreen", "LR_PET1" = "lightgreen",
                      "Rmst_PET0" = "blue", "Rmst_PET1" = "blue")

    a_power_delta <- data.frame(data[, c(1,2,3,4,5,6,7)])
    colnames(a_power_delta) <- c(var_name,'LR_alpha','Rmst_alpha','Our_alpha',
                          'LR_power', 'Rmst_power','Our_power')
    a_power_long <- a_power_delta %>%
        pivot_longer(cols = -!!sym(var_name), names_to = "variable", values_to = "value")%>%
        mutate(linetype_group = ifelse(variable %in% 
            c("LR_alpha", "Rmst_alpha", "Our_alpha"), "Alpha", "Power"))

    plot1 <- ggplot(a_power_long, aes(x = !!sym(var_name), y = value, 
        color = variable, linetype = linetype_group)) +
    geom_point(size = 3) +
    geom_line(linewidth = 1) +
    scale_linetype_manual(values = c("Alpha" = "solid", "Power" = "dotted")) +
    labs( linetype = "Line Type", color = "Variable",
          title = "Line Plot with Different Line Types") +

    scale_y_continuous(breaks = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), limits = c(0, 1)) +
    scale_color_manual(values = color_palette) +
    labs(x = var_name, y = "Value", color = "Variable",
      title = 'Type I error and Power') +
    theme_minimal(base_size = 18) + 
    theme(plot.title = element_text(hjust = 0.5),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    guides(linetype = guide_legend(override.aes = list(color = "black")),
        color = guide_legend(override.aes = list(linetype = "solid")))

#_________next plot_________
    pet_delta <- data.frame(data[, c(1,8,9,10,11,12,13)])
    colnames(pet_delta) <- c(var_name,'LR_PET0', 'Rmst_PET0', 'Our_PET0',
                      'LR_PET1', 'Rmst_PET1', 'Our_PET1')
    pet_long <- pet_delta %>%
        pivot_longer(cols = -!!sym(var_name), names_to = "variable", values_to = "value")%>%
        mutate(linetype_group = ifelse(variable %in% 
              c("LR_PET0", "Rmst_PET0", "Our_PET0"), "PET0", "PET1"))

    plot2 <- ggplot(pet_long, aes(x = !!sym(var_name), y = value, 
        color = variable, linetype = linetype_group)) +
        geom_point(size = 3) +
        geom_line(linewidth = 1) +
        scale_linetype_manual(values = c("PET0" = "solid", "PET1" = "dotted")) +
        labs( linetype = "Line Type", color = "Variable",
          title = "Line Plot with Different Line Types") +

        scale_y_continuous(limits = c(0, 1)) +
        scale_color_manual(values = color_palette) +
        labs(x = var_name, y = "Value", color = "Variable",
        title = 'PET0 and PET1') +
        theme_minimal(base_size = 18) + 
        theme(plot.title = element_text(hjust = 0.5),
              plot.background = element_rect(fill = "white", color = NA),
              plot.margin = unit(c(1, 1, 1, 1), "cm")) + 
        guides(linetype = guide_legend(override.aes = list(color = "black")),
              color = guide_legend(override.aes = list(linetype = "solid")))

    plot_grid(plot1, plot2, ncol = 2)

  }