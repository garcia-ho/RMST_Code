# RMST Two Stages Clinical Trial Design
This is the storage of the start up project.  
Hope I could publish a paper based on this for my Ph.D. application.  

## Compare 3 methods:
- 2 stages log rank test design (Jung 2017)   
    -  [Minjung Kwak & Sin-Ho Jung (2017) Optimal two-stage log-rank test for randomized phase II clinical trials, Journal of Biopharmaceutical Statistics, 27:4, 639-658, DOI: 10.1080/10543406.2016.1167073](https://www.tandfonline.com/doi/abs/10.1080/10543406.2016.1167073)  
    

-  Simple RMST difference (Eaton 2020, Shan 2021)    
    -  [Eaton A, Therneau T, Le-Rademacher J. Designing clinical trials with (restricted) mean survival time endpoint: Practical considerations. Clinical Trials. 2020;17(3):285-294. doi:10.1177/1740774520905563](https://journals.sagepub.com/doi/abs/10.1177/1740774520905563)


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

****
## Asymptotic Normality of RMST:
```math
[E_1-C_1,\ E_1,\ E_2-C_2,\ E_2] \sim \mathcal{N}\ (\ \hat{\mu},\ \hat{\Sigma}\ )
```

```math
\hat{\Sigma} = 
 \left[
 \begin{matrix}
   Var(E_1)+Var(C_1) & Var(E_1) & Cov(E_1,E_2)+Cov(C_1,C_2) & Cov(E_1,E_2) \\
                   & Var(E_1) & Cov(E_1,E_2) & Cov(E_1,E_2) \\
                   &         & Var(E_2)+Var(C_2) & Var(E_2)\\
                   &         &                 & Var(E_2)\\
  \end{matrix}
\right]
```

### Estimated by Monte Carlo Simulation
Similar to Lu(2021):  
[Lu Y, Tian L. Statistical considerations for sequential analysis of the restricted mean survival time for randomized clinical trials[J]. Statistics in biopharmaceutical research, 2021, 13(2): 210-218.](https://www.tandfonline.com/doi/abs/10.1080/19466315.2020.1816491)

```math
 \left[
 \begin{matrix}
   Var(\hat{E_1}) & Cov(\hat{E_1}, \hat{E_2}) \\
   Cov(\hat{E_1}, \hat{E_2}) & Var{\hat{E_2}} \\
  \end{matrix}
  \right] = 

\frac{1}{B}\sum\limits_{i = 1}^{B} 
 \left[
 \begin{matrix}
   \hat{RMST_{E1}} - \bar{RMST_{E1}} \\
   \hat{RMST_{E2}} - \bar{RMST_{E2}} \\
  \end{matrix}
  \right] ^{\bigotimes2}_i  
 ```  
  
```math
  \bar{RMST_{E1}} = \frac{1}{B}\sum^{B}_{i = 1}RMST_{E1,i}
```  
The theoretical RMST value did not consider accrual and censoring distribution.

****

## Parameter Optimization (Grid Search)
Reference: Zhou(2017) BOP2 Bayesian design:  
[Zhou H, Lee J J, Yuan Y. BOP2: Bayesian optimal design for phase II clinical trials with simple and complex endpoints[J]. Statistics in medicine, 2017, 36(21): 3302-3314.](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.7338)  


In order to solve the critical values ($m_1,\ t_1,\ m_2,\ t_2$), a function that can control the normal probability is required. The following $\mathcal{f}(n)$ is what we proposed. 

```math
\large \mathcal{f}(\Tilde{N}) = \mathcal{e}^{-\ \gamma · \frac{\Tilde{N}}{N}}
```

$\Tilde{N}$ is the sample size(2 arms) of interim period. $N$ is the final total sample size of 2 arms. $\hat{D}_i = \hat{R}_{Ei} - \hat{R}_{Ci}$
Then we set the following constraints:   

```math
\begin{aligned}
P(\hat{D}_1(\tau_1)> m_1) &= f(\Tilde{N}) \\
P(\hat{D}_1(\tau_1) > m_1 \cap \hat{R}_{E1}(\tau_1)> q_1) &= \lambda \cdot f(\Tilde{N}) \\
P(\hat{D}_2(\tau_2) > m_2) &= f(N) \\
P(\hat{D}_2(\tau_2) > m_2 \cap \hat{R}_{E2}(\tau_2)> q_2) &= \lambda \cdot f(N) \\
\lambda \in (0,1).\\ 
\end{aligned}
```

$\mathcal{f}(\tilde{N})$ is a monotonously decereasing funciton of n, which means that two probability constraints in interim period will go up when the interim sample size n decrease.  
#### It leads to a small early stop probability with an insufficient interim sample size.  
Then we grid search $(\lambda, \gamma)$ . Each pair of $(\lambda, \gamma)$ determines a set of ($m_1,\ t_1,\ m_2,\ t_2$) by normal calculation. Record critical values sets that yield the desirable overall type I error $\alpha$:  

```math  
\large \alpha = P(\hat{D}_1(\tau_1) > m_1\ \cap\ \hat{R}_{E1}(\tau_1)> q_1\ \cap\ \hat{D}_2(\tau_2) > m_2\ \cap\ \hat{R}_{E2}(\tau_2)> q_2\ |\ H_0),\\
```

2 grid searching strategies can be conducted:
#### 1. Two-stage Most Powerful Design with Fixed Sample Size
#### 2. Two-stage Optimal Design Minimizing $\overline{EN}$ 


****
## Result
### Under Proportional Hazard Scenario  
We used $\frac{\lambda{E}}{\lambda{C}}= 0.67$ as base setting. The Type I error, Power, Probability of Early Termination under H0 and H1 are shown in the graph. We apply constant accrual rate $r = 90$(patients per year) with a fixed overall sample size N = 200 (100 per arm) . The longer the interim period, the larger the interim sample size n.  
The performances of three methods are compared with different Hazard Ratio(under H1) and different interim period. Our methods demonstrateed competitive performance compared to log-rank test and simple RMST test without extra rejection restriction.

<p align="center"> 
<img src="figures/Comparison/PH.png" alt="Result Under PH" width="900" height="600">
</p>


### Under Early Difference Scenario
The survival function of two groups are shown below: 
<p align="center"> 
<img src="figures/Comparison/Early_diff_sur.png" alt="Survival of Early Difference" width="420" height="300">
</p>

The comparison of three methods under different setting is shown below:  
![Result Under Early Difference Scenario](figures/Comparison/Early_diff.png)  

The graph shows us that our method have higher power while preserving the type I error under different early difference settings. Higher PET0 and relatively acceptable sacrifice of PET1. 


****
### Visualization of Sculpted Rejection Region 

![The sculpted critical region](figures/rejection_plot.jpg)  

The data points in overall rejection region of H0 in 10000 simulated trials. Each blue or orange dot in subplot A representing a pair of $(R_{E1}, R_{C1})$, which is one simulated trial in interim stage, has a corresponding dot with the same color in subplot B indicating the final stage of the same trial. A pair of corresponding dots that both fall in shaded region will be dyed green. Lines are most powerful decision boundaries of Simple and Sculpted RMST in two stages.  
185 blue dots and 99 orange dots move out of the shaded region from interim to final, which is the main source of difference in overall type one error. More precisely, there are more simulated trials with small experimental RMST in final stage $R_{E2}$ fail to reject $H_0$ due to the $R_{E2} > 0.8983$ condition in Sculpted RMST (blue dots to the left of shaded area) than that with low final RMST difference $D_2$ (orange dots below shaded area). These additional dots moved outside of the sculpted region compared to Simple RMST represent the saving type I error.   
Move on to the discussion of parameters $(\gamma, \lambda)$ in grid search procedure. $\gamma$ influences the relative position of the orange lines in subplot A compared to those in subplot B. Meanwhile, $\lambda$ governs the relative position of the vertical and horizontal orange lines within each subplot. 






****
### Robustness
Within the hypothesis of a trial: 
```math
\begin{aligned}
H_0&:\ \lambda_E = \lambda_C = \lambda_0,\\
H_1&:\ \lambda_C = \lambda_0,\ \lambda_E = \lambda_1, \nonumber
\end{aligned}
```
The following graph shows the robustness of three methods when the real hazard drift away from the assumed value. The critical values calculated under this setting of orange dotted lines are saved then plugged into the simulated data under different $\lambda_0$(subplot A) or $\lambda_1$(subplot B). 10000 pairs of experimental and control groups survival data are generated under each $\lambda_0$ or $\lambda_1$ while preserving the hazard of the other group. The change of Power and PET1 is shown below(The PET0 and type I error should be consistant across different HR).  The actual $\alpha$ and power of the Sculpted RMST will deviate more from the expected stated level compared to the other two methods if the assumed hazard in control group is inaccurate. However if only $\lambda_1$ drifts from the stated number, our method can perform as robustly as the Simple RMST and log-rank test.

![Robustness under HR](figures/Robustness.jpg)  




****
All codes are in R. Presented in notebook R kernel.  
All functions that used for simulation is stored at [Function.R](Rfiles/Function.R).  
The simulation processes of single stage and two stages are stored at different ipynb files.


--------------
Vamos García

