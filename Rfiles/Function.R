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





#------ 3. RMST_sim_test-------
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





#---------------------- 9. cov_mc -----------------------
# This function use Monte Carlo simulation to calculate the var_cov matrix of 
#   [E1-C1, E1, E2-C2, E2]
# The calculation formula refers to the Lu(2021) sequential trials.

#*************
# The input rmst_int should be the interim rmst data of two groups (RMST_sim_cal output)
# rmst_fin is the final rmst data of two groups. sim_size is the MC simulation times B
# true_rmst_int is the theoretical result of interim rmst of two groups c(RMST_C, RMST_E)

cov_mc <- function(rmst_int, rmst_fin, true_rmst_int, true_rmst_fin, sim_size){

    diff_C <- rbind(rmst_int[1, ] - true_rmst_int[1], rmst_fin[1, ] - true_rmst_fin[1])
    cov_C <- matrix(0, nrow = 2, ncol = 2)
    for (i in 1:sim_size) 
      {
        product <- diff_C[, i] %*% t(diff_C[, i])  
        cov_C <- cov_C + product  
      }
    cov_C <- cov_C / sim_size   # [ Var(C1)  Cov(C1, C2)
                                #  Cov(C1, C2)  Var(C2) ]

    diff_E <- rbind(rmst_int[2,] - true_rmst_int[2], rmst_fin[2,] - true_rmst_fin[2])
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









#------------------9. find_m_t_RMST------------------
# Works for both Our RMST test and simple RMST difference test
# Given rmst data of interim and all, find the best m1,t1, m2,t2
# n need to be given in this loop
# if t1, t2 is abandon, t_low need to be set as -Inf, It will only return m1,m2

find_m_t_RMST <- function(m_low, t_low, t_up, rmst_data, search_times, search_step,
                     tar_a1, tar_pow1_low, tar_a2, sim_size) {
  rmst_h0_int <- rmst_data[c(1,2) , ]
  rmst_h1_int <- rmst_data[c(3,4) , ]
  rmst_h0_all <- rmst_data[c(1,2,5,6) , ]
  rmst_h1_all <- rmst_data[c(3,4,7,8) , ]
  result_m1_t1 <- c()
  result_m1_t1 <- foreach(i = 1:search_times, .combine = 'cbind') %dopar% { 
      m1 = m_low + i * search_step
      result_t1 <- c()
      if (t_low == -Inf) 
        {   # for simple RMST difference test, No second condition
          proc_h0 <- sum((rmst_h0_int[2, ] - rmst_h0_int[1, ] > m1))
          proc_h1 <- sum((rmst_h1_int[2, ] - rmst_h1_int[1, ] > m1))
          if (proc_h0/sim_size < tar_a1 & proc_h1/sim_size >= tar_pow1_low ) 
            {
              result_t1 <- c(m1,-Inf,proc_h0/sim_size,proc_h1/sim_size)         
            }
        }
      else 
        {
          for (t1 in seq(from = t_low, to = t_up, by = (t_up - t_low) /  search_times))
          {
            proc_h0 <- sum((rmst_h0_int[2, ] - rmst_h0_int[1, ] > m1) & (rmst_h0_int[2, ] > t1))
            proc_h1 <- sum((rmst_h1_int[2, ] - rmst_h1_int[1, ] > m1) & (rmst_h1_int[2, ] > t1))
            if (proc_h0/sim_size < tar_a1 & proc_h1/sim_size >= tar_pow1_low )
              {
                result_t1 <- cbind(result_t1, c(m1, t1, proc_h0/sim_size, proc_h1/sim_size))
              }    
          }
        }
      result_t1
      }

    if (is.null(result_m1_t1)) {
      # Return NULL when something goes wrong
        return(NULL)
      }

    powerful_m1_t1 <- result_m1_t1[, which(result_m1_t1[4,] == max(result_m1_t1[4,]))]
    #Find the most powerful m1,t1

    if (t_low == -Inf | class(powerful_m1_t1)[1] == 'numeric') {  #powerful_m1_t1 is a 1*4 vector
      bestmt <- powerful_m1_t1
    }
    else {
      bestmt <- powerful_m1_t1[, which(powerful_m1_t1[2,] == max(powerful_m1_t1[2,])) ]
    } # return the result with the larges t1
    # Among these most powerful, find the m1 with smallest absolute value 
    m1 <- bestmt[1]
    t1 <- bestmt[2]  # t1 = -Inf when simple RMST difference test

    result_fin <- c()
    result_fin <- foreach(i = 1:search_times, .combine = 'cbind') %dopar% 
      { 
        m2 = m_low + i * search_step
        opt_power <- 0
        opt_m2 <- 0
        opt_t2 <- 0
        opt_mt <- c()
        if (t1 == -Inf) 
          {
            proc_h0 <- sum((rmst_h0_all[2, ] - rmst_h0_all[1, ] > m1) &
                          (rmst_h0_all[4, ] - rmst_h0_all[3, ] > m2))
            proc_h1 <- sum((rmst_h1_all[2, ] - rmst_h1_all[1, ] > m1) & 
                          (rmst_h1_all[4, ] - rmst_h1_all[3, ] > m2))
            if ( proc_h0/sim_size < tar_a2 & proc_h1/sim_size >= opt_power) 
                {
                  opt_power <- proc_h1/sim_size
                  opt_m2 <- m2
                  opt_t2 <- -Inf
                  opt_mt <- c(opt_m2,opt_t2,proc_h0/sim_size,opt_power)
                }
            }
      else 
        {
          for (t2 in seq(from = t_low, to = t_up, by = (t_up - t_low) /  search_times)) 
            {
              proc_h0 <- sum((rmst_h0_all[2, ] - rmst_h0_all[1, ] > m1) & (rmst_h0_all[2, ] > t1) &
                            (rmst_h0_all[4, ] - rmst_h0_all[3, ] > m2) & (rmst_h0_all[4, ] > t2))     
              proc_h1 <- sum((rmst_h1_all[2, ] - rmst_h1_all[1, ] > m1) & (rmst_h1_all[2, ] > t1) &
                            (rmst_h1_all[4, ] - rmst_h1_all[3, ] > m2) & (rmst_h1_all[4, ] > t2))
                  # return the best 
              if ( proc_h0/sim_size < tar_a2 & proc_h1/sim_size >= opt_power) 
                {
                  opt_power <- proc_h1/sim_size
                  opt_m2 <- m2
                  opt_t2 <- t2
                  opt_mt <- c(opt_m2, opt_t2, proc_h0/sim_size, opt_power)
                }
            }
        }
      opt_mt
      }

  if (is.null(result_fin)) 
    {
      return(NULL)    # Return NULL when something goes wrong
    }

  powerful_fin <- result_fin[, which(result_fin[4,] == max(result_fin[4,]))]
  if (t1 == -Inf | class(powerful_fin)[1] == 'numeric')
    {
      best_result <- powerful_fin
    }
  else 
    {
      best_result <- powerful_fin[ , which(powerful_fin[2,] == max(powerful_fin[2,])) ]
    }    
    return(data.frame(m1 = m1,
                  t1 = t1,
                  PET0 = 1 - bestmt[3],
                  PET1 = 1 - bestmt[4],
                  m2 = best_result[1],
                  t2 = best_result[2],
                  alpha = best_result[3],
                  Power = best_result[4]
                  ))
}





#------------------10. find_m_logrank------------------
# Similar to find_m_t_RMST, this one is for log rank test 2 stages design
# The z statistics (W/sigma) in logrank test of each simulation is required as an input
# The z statistics can be obtained by log_rank_sim function $ z_stats
# This function help you find the Z1 > m1 & Z > m2 to control the overall power.
# You need to tune tar_a1 and tar_pow1_low to control PET0 and PET1

find_m_logrank <- function(m_low, logrank_data, search_times, search_step,
                          tar_a1, tar_pow1_low, tar_a2, sim_size) {
  z_stats_h0_int <-  logrank_data[1, ]
  z_stats_h1_int <-  logrank_data[2, ]
  z_stats_h0_all <-  logrank_data[c(1, 3), ]
  z_stats_h1_all <-  logrank_data[c(2, 4), ]          
  result_m1 <- foreach(i = 1:search_times, .combine = 'cbind') %dopar% { 
      m1 = m_low + i * search_step
      result_t1 <- c()
      proc_h0 <- sum(z_stats_h0_int > m1)
      proc_h1 <- sum(z_stats_h1_int > m1)
      if (proc_h0/sim_size > 0 & proc_h0/sim_size < tar_a1 &
          proc_h1/sim_size >= tar_pow1_low ) 
        {
          result_t1 <- c(m1,proc_h0/sim_size,proc_h1/sim_size)         
        }
      result_t1
      }
  powerful_m1 <- result_m1[, which(result_m1[3,] == max(result_m1[3,]))]
    #Find the most powerful m1,t1
  if (is.null(powerful_m1)) 
    {   # Return NULL when something goes wrong
      return(NULL)
    }
  bestmt <- powerful_m1
  m1 <- bestmt[1]
  result_fin <- c()
  result_fin <- foreach(i = 1:search_times, .combine = 'cbind') %dopar% { 
      m2 = m_low + i * search_step
      opt_alpha <- 1
      opt_power <- 0
      opt_m2 <- 0
      opt_mt <- c()
      proc_h0 <- sum((z_stats_h0_all[1, ] > m1) & (z_stats_h0_all[2, ]  > m2))
      proc_h1 <- sum((z_stats_h1_all[1, ] > m1) & (z_stats_h1_all[2, ]  > m2))
      if (proc_h0/sim_size > 0 & proc_h0/sim_size < tar_a2 & proc_h1/sim_size >= opt_power) {
          opt_alpha <- proc_h0/sim_size
          opt_power <- proc_h1/sim_size
          opt_m2 <- m2
          opt_mt <- c(opt_m2,opt_alpha,opt_power)
        }
      opt_mt
  }
  if (is.null(result_fin)) 
    {
      return(NULL)
    }
  else if (class(result_fin)[1] == 'numeric') 
    {
      powerful_fin <- result_fin
    }
  else 
    {
      powerful_fin <- result_fin[, which(result_fin[3,] == max(result_fin[3,]))]
    }
  return(data.frame(m1 = m1,
                  PET0 = 1 - bestmt[2],
                  PET1 = 1 - bestmt[3],
                  m2 = powerful_fin[1],
                  alpha = powerful_fin[2],
                  Power = powerful_fin[3]
                  ))
}



