# RMST Two Stages Clinical Trial Design
This is the storage of the start up project.  
Hope I could publish a paper based on this for my Ph.D. application.  

## Compare 4 methods:
- 2 stages log rank test design (Jung 2017)   
    -  [Minjung Kwak & Sin-Ho Jung (2017) Optimal two-stage log-rank test for randomized phase II clinical trials, Journal of Biopharmaceutical Statistics, 27:4, 639-658, DOI: 10.1080/10543406.2016.1167073](https://www.tandfonline.com/doi/abs/10.1080/10543406.2016.1167073)  
    

-  Simple RMST difference (Eaton 2020, Shan 2021)    
    -  [Eaton A, Therneau T, Le-Rademacher J. Designing clinical trials with (restricted) mean survival time endpoint: Practical considerations. Clinical Trials. 2020;17(3):285-294. doi:10.1177/1740774520905563](https://journals.sagepub.com/doi/abs/10.1177/1740774520905563)

-  Multiple Direction Log rank test (mdir:  Ditzhaus 2020)  
    -  [Ditzhaus, M., Pauly, M. (2018). Wild bootstrap logrank tests with broader power functions for testing superiority. arXiv preprint arXiv:arXiv:1808.05627.](https://www.sciencedirect.com/science/article/abs/pii/S0167947319300362)

-  Special RMST rejection region (inspired by Litwin 2017)  
    -  [Litwin S, Basickes S, Ross E A. Two‐sample binary phase 2 trials with low type I error and low sample size[J]. Statistics in medicine, 2017, 36(9): 1383-1394.](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7226)



Basically, our method is testing both RMST difference and the RMST of experiment group.  
## Hypothesis
$H_0:\ \lambda_E(t) > \lambda_C(t)$, where $\lambda$ is the hazard function  
$H_1:\ \lambda_E(t) < \lambda_C(t)$  

## Two stages trial rejection region (one-sided):
- 1. Log rank test
    - $\frac{W_1}{\sigma_1}>c_1$ & $\frac{W}{\sigma}>c$, where $\frac{W}{\sigma}>c$ is the normal log rank test statistics    
According to Jung(2017), two stages log-ranktest can reach a similar power as corresponding single stage test if the interim period, $c_1,\ c$ are well selected. 

$E, C$ are the RMST value of experiment group and the control group respectively.   
$\tau$ is the cutoff time point of interim period(1) and overall trial(2)  

- 2. Simple RMST Difference  
    - $E(\tau_1) - C(\tau_1) > m_1$ & $E(\tau_2)-C(\tau_2)>m_2$

- 3. Our RMST Rejection method
    - $E(\tau_1)-C(\tau_1)>m_1$ & $E(\tau_1)>t_1$ & $E(\tau_2)-C(\tau_2)>m_2\$ & $E(\tau_2)>t_2$  

The critical values of $m_1,t_1,m_2,t_2$ are calculated by 10000 times Monte Carlo simulation
  

All codes are in R. Presented in notebook R kernel.  
All functions that used for simulation is stored at [Function.R](Rfiles/Function.R).  
The simulation processes of single stage and two stages are stored at different ipynb files.
No Bayesian methods would be compared with. 

--------------
Vamos García

