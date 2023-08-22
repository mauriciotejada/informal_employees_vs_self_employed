# Replication files for: Are informal self-employment and informal employment as employee behaviorally distinct labor force states?

## Reference

Luca Flabbi and Mauricio M. Tejada, "Are informal self-employment and informal employment as employee behaviorally distinct labor force states?", Economics Letters, Volume 231, October 2023, 111278

[Link to the paper](https://doi.org/10.1016/j.econlet.2023.111278).

## Abstract

The paper performs both a parametric and non-parametric analysis to address a fundamental question in the growing literature using search models to study labor market informality: Should informal self-employment and informal employment as an employee be considered two different labor market
states? Both analyses strongly reject equality between the two states, cautioning against aggregating them in a common "informality state". The parametric model identifies that the variation in informal self-employment income and the short duration of informal employee jobs are the primary factors that contribute to the observed differences between these labor market states.

## Replication instructions

To replicate the results run `bash run.sh` in the terminal. Note that *Julia 1.8 with instantiated environment is required*.

All the results where generated using a *MacBook Pro M2 with 16 GB of RAM and runing macOS 13.4*.

Note: Bootstrapped results may sligtly differ from those in the paper because the use of parallel computing. Eventought a random seed was used, the order of the working processor in each iteration is still random. 

All the results are stored in the `/res/` folder. The naming convetion is as follows:

- `est_param_0.txt`: Estimated parameters of the unrestricted model.
- `est_param_r_0.tx`: Estimated parameters of the restricted model.
- `likfun_0.txt`: Value of the loglikelihood function of the unrestricted model.
- `likfun_r_0.txt`: Value of the loglikelihood function of the restricted model.
- `imp_values_0.txt`: Implied parameters of the unrestricted model (hazard rates and offered wages stats).
- `imp_values_r_0.txt`: Implied parameters of the restricted model (hazard rates and offered wages stats).
- `est_param_boot_0.txt`: Estimated parameters of the unrestricted model in each iteration of the bootstrap.
- `stderror_boot_0.txt`: Bootstrapped standard errores of the estimated parameters of the unrestricted model.
- `est_param_boot_r_0.txt`: Estimated parameters of the restricted model in each iteration of the bootstrap.
- `stderror_boot_r_0.txt`: Bootstrapped standard errores of the estimated parameters of the restricted model.
- `est_param_boot_b_0.txt`: Estimated b parameter of the unrestricted model in each iteration of the bootstrap.
- `stderror_boot_b_0.txt`: Bootstrapped standard errores of the estimated b parameter of the unrestricted model.
- `est_param_boot_b_r_0.txt`: Estimated b parameter of the restricted model in each iteration of the bootstrap.
- `stderror_boot_b_r_0.txt`: Bootstrapped standard errores of the estimated b parameter of the restricted model.