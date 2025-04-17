##load packages
library(dplyr)
library(rbi)
library(rbi.helpers)
library(ggplot2)
library(tinytex)
library(latex2exp)
library(ggpubr)
library(pander)
library(lubridate)

#load the data
v <- read.csv("~/Desktop/Project IV/andre_estimates_21_02.txt", sep  = "\t") %>%
  rowSums()
y <- data.frame(value = v) %>%
  mutate(time = seq(7, by = 7, length.out = n())) %>%
  dplyr::select(time, value)
ncores <- 8
minParticles <- max(ncores, 16)

#define the model
dureau_model_str <- "
model dureau {
  obs y

  state S
  state E
  state I
  state R
  state x

  state Z

  input N
  param k
  param gamma
  param sigma // Noise driver
  param E0
  param I0
  param R0
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(1.59, 0.02, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(1.08, 0.075, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    x0 ~ uniform(-5,2)
    I0 ~ uniform(-16, -9)
    E0 ~ uniform(-16, -9)
    R0 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    tau ~ uniform(0, 1)
  }

  sub initial {
    S <- N
    R <- R0*S
    S <- S - R

    E <- exp(E0 + log(S))
    S <- S - E
    I <- exp(I0 + log(S))
    S <- S - I
    x <- x0
    Z <- 0
  }

  sub transition(delta = 1) {
    Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {
    y ~ log_normal(log(max(Z/10.0, 0)), tau)
  }

  sub proposal_parameter {
    k ~ gaussian(k, 0.005)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    x0 ~ gaussian(x0, 0.05)
    E0 ~ gaussian(E0, 0.05)
    I0 ~ gaussian(I0, 0.05)
    R0 ~ gaussian(R0, 0.05)
    tau ~ gaussian(tau, 0.05)
  }
}"



#perform PMCMC 
dureau_model <- bi_model(lines = stringi::stri_split_lines(dureau_model_str)[[1]])
dureau_bi_model <- libbi(dureau_model)
input_lst <- list(N = 52196381)
end_time <- max(y$time)
obs_lst <- list(y = y %>% dplyr::filter(time <= end_time))

dureau_bi <- sample(dureau_model, end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, 
                    nparticles = minParticles, nthreads = ncores, proposal = 'prior', verbose = TRUE) %>% 
  adapt_particles(min = minParticles, max = minParticles*200) %>%
  adapt_proposal(min = 0.05, max = 0.4) %>%
  sample(nsamples = 5000, thin = 2) %>% # burn in 
  sample(nsamples = 5000, thin = 2)

#extract the samples
dureau_bi_lst <- bi_read(dureau_bi %>% sample_obs)

#plot incidence vs. time
fitY <- dureau_bi_lst$y %>% 
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(y %>% rename(Y = value))

g1 <- ggplot(data = fitY) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = Y), colour = "Red") +
  ylab("Incidence") +
  xlab("Time (days)")+
  theme(
    axis.title.x = element_text(size = 18),  # Change x-axis label size
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
  )

#plot transmission vs. time
plot_df <- dureau_bi_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()


g2 <- ggplot(data = plot_df) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75),fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Transmissibility ($\\beta(t)$)")) +
  xlab("Time (days)")+
  theme(
    axis.title.x = element_text(size = 18),  # Change x-axis label size
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
  )


#plot Rt vs. time

Rt <- c()
for(i in 1:length(dureau_bi_lst$x$value)){
  x = dureau_bi_lst$x$value[i]
  s = dureau_bi_lst$S$value[i]
  b = exp(x)
  N = 52196381
  g = dureau_bi_lst$gamma$value[dureau_bi_lst$x$time[i]]
  Rt[i] <- (b*s)/(g*N)
}

Rt_df <- data.frame(time = dureau_bi_lst$x$time, Rt = Rt)

plot_df <- Rt_df %>%
  group_by(Rt) %>%
  group_by(time)%>%
  mutate(
    q025 = quantile(Rt, 0.025),
    q25 = quantile(Rt, 0.25),
    q50 = quantile(Rt, 0.5),
    q75 = quantile(Rt, 0.75),
    q975 = quantile(Rt, 0.975)
  ) %>% ungroup()

g3 <- ggplot(data = plot_df) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75),fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_hline(yintercept=1, linetype = "dashed", color = "red")+
  ylab(expression(R[t])) +
  xlab("Time (days)")+ theme(
    axis.title.x = element_text(size = 18),  # Change x-axis label size
    axis.title.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
  )

ggarrange(g1,g2,g3, ncol = 3, nrow = 1, align = "v")


################################################################################
## FORECASTING ##
################################################################################

## start by forecasting around time step 50 to see how well it can account for the change in the dynamics for the second wave
n_forecast_samples <- 100

last_time <- max(dureau_bi_lst$S$time)

# Extract final state estimates using median (or mean)
# S_last <- mean(dureau_bi_lst$S$value[dureau_bi_lst$S$time == last_time])
# E_last <- mean(dureau_bi_lst$E$value[dureau_bi_lst$E$time == last_time])
# I_last <- mean(dureau_bi_lst$I$value[dureau_bi_lst$I$time == last_time])
# R_last <- mean(dureau_bi_lst$R$value[dureau_bi_lst$R$time == last_time])
# x_last <- mean(dureau_bi_lst$x$value[dureau_bi_lst$x$time == last_time])

posterior_states <- list("S", "E", "I", "R", "x") %>%
  set_names() %>%
  map(~ {
    dureau_bi_lst[[.x]] %>%
      filter(time == last_time) %>%
      slice_sample(n = n_forecast_samples)
  })

param_names <- c("k", "gamma", "sigma", "tau", "x0")
posterior_params <- param_names %>%
  set_names() %>%
  map(~ {
    bi_read(dureau_bi, type = "param")[[.x]] %>%
      slice_sample(n = n_forecast_samples)
  })




dureau_model_forecast_str <- "
model dureau {
  obs y

  state S
  state E
  state I
  state R
  state x

  state Z

  input N
  input S_init
  input E_init
  input I_init
  input R_init
  input x_init
  param k
  param gamma
  param sigma // Noise driver
  param x0
  param tau

  sub parameter {
    k ~ truncated_gaussian(1.59, 0.02, lower = 0) // k is the period here, not the rate, i.e. 1/k is the rate
    gamma ~ truncated_gaussian(1.08, 0.075, lower = 0) // gamma is the period, not the rate
    sigma ~ uniform(0,1)
    x0 ~ uniform(-2,5)
    tau ~ uniform(0, 1)
  }

  sub initial {
    S <- S_init
    E <- E_init
    I <- I_init
    R <- R_init
    Z <- 0
    x <- x_init
  }

  sub transition(delta = 1) {
    Z <- ((t_now) % 7 == 0 ? 0 : Z)
    noise e
    e ~ wiener()
    ode(alg = 'RK4(3)', h = 1.0, atoler = 1.0e-3, rtoler = 1.0e-8) {
      dx/dt = sigma*e
      dS/dt = -exp(x)*S*I/N
      dE/dt = exp(x)*S*I/N - E/k
      dI/dt = E/k-I/gamma
      dR/dt = I/gamma
      dZ/dt = E/k
    }
  }

  sub observation {} 

  sub proposal_parameter {
    k ~ gaussian(k, 0.005)
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    x0 ~ gaussian(x0, 0.05)
    tau ~ gaussian(tau, 0.05)
  }
}"

input_forecast <- list(
  N = 52196381,        # Total population
  S_init = S_last,  # Last estimated S
  E_init = E_last,
  I_init = I_last,  # Last estimated I
  R_init = R_last,  # Last estimated R
  x_init = x_last   # Last estimated x
)


dureau_forecast_model <- bi_model(lines = stringi::stri_split_lines(dureau_model_forecast_str)[[1]])
dureau_bi_forecast <- libbi(dureau_forecast_model)
forecast_horizon <- 28  # Forecast 3 weeks ahead
new_end_time <- last_time + forecast_horizon
dureau_bi_forecast <- sample(dureau_bi_forecast, start_time = last_time, end_time = new_end_time, input = input_forecast, obs = NULL,
                             nsamples = 1000, nparticles = dureau_bi$options$nparticles, nthreads = ncores, output_every = 1,
                             proposal = 'model', target = 'prediction', verbose = TRUE) %>%
                      sample(nsamples = 5000, thin = 5) %>% # burn in 
                      sample(nsamples = 5000, thin = 5)



dureau_bi_forecast_lst <- bi_read(dureau_bi_forecast %>% sample_obs)

plot_df_beta <- dureau_bi_lst$x[dureau_bi_lst$x$time <=new_end_time,] %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()

plot_df_beta_forecast <- dureau_bi_forecast_lst$x %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.05),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.95)
  ) %>% ungroup()



g_trans_comb <- ggplot() +
  geom_ribbon(data = plot_df_beta, aes(x = time, ymin = q25, ymax = q75), alpha = 0.3, fill = "royalblue3") +
  geom_ribbon(data = plot_df_beta, aes(x = time, ymin = q025, ymax = q975), alpha = 0.3, fill = "cornflowerblue") +
  geom_line(data = plot_df_beta,aes(x = time, y = q50)) +

  geom_ribbon(data = plot_df_beta_forecast, aes(x = time, ymin = q25, ymax = q75), alpha = 0.4, fill = "darkolivegreen3") +
  #geom_ribbon(data = plot_df_beta_forecast, aes(x = time, ymin = q025, ymax = q975), alpha = 0.3, fill = "darkolivegreen3") +
  geom_line(data = plot_df_beta_forecast,aes(x = time, y = q50), color = "red") +
  geom_vline(xintercept = last_time, color = "blue", linetype = "dashed")+
  ylab(TeX("Transmissibility")) +
  xlab("Time")
g_trans_comb 


Rt_df <- data.frame(time = dureau_bi_lst$x$time, Rt = Rt)

plot_df <- Rt_df[dureau_bi_lst$x$time <=new_end_time,] %>%
  group_by(Rt) %>%
  group_by(time)%>%
  mutate(
    q025 = quantile(Rt, 0.025),
    q25 = quantile(Rt, 0.25),
    q50 = quantile(Rt, 0.5),
    q75 = quantile(Rt, 0.75),
    q975 = quantile(Rt, 0.975)
  ) %>% ungroup()

Rt_f <- c()
for(i in 1:length(dureau_bi_forecast_lst$x$value)){
  x = dureau_bi_forecast_lst$x$value[i]
  s = dureau_bi_forecast_lst$S$value[i]
  b = exp(x)
  N = 52196381
  g = dureau_bi_forecast_lst$gamma$value[dureau_bi_forecast_lst$x$time[i]]
  Rt_f[i] <- (b*s)/(g*N)
}


Rt_forecast_df <- data.frame(time = dureau_bi_forecast_lst$x$time, Rt = Rt_f)

plot_df_f <- Rt_forecast_df %>%
  group_by(Rt) %>%
  group_by(time)%>%
  mutate(
    q025 = quantile(Rt, 0.05),
    q25 = quantile(Rt, 0.25),
    q50 = quantile(Rt, 0.5),
    q75 = quantile(Rt, 0.75),
    q975 = quantile(Rt, 0.95)
  ) %>% ungroup()


g3 <- ggplot(data = plot_df) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75),fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_hline(yintercept=1, linetype = "dashed", color = "red")+
  ylab(expression(R[0])) +
  xlab("Time")


g_rt_comb <- ggplot() +
  geom_ribbon(data = plot_df, aes(x = time, ymin = q25, ymax = q75), alpha = 0.3, fill = "royalblue3") +
  geom_ribbon(data = plot_df, aes(x = time, ymin = q025, ymax = q975), alpha = 0.3, fill = "cornflowerblue") +
  geom_line(data = plot_df,aes(x = time, y = q50)) +
  
  geom_ribbon(data = plot_df_f, aes(x = time, ymin = q25, ymax = q75), alpha = 0.4, fill = "darkolivegreen3") +
  #geom_ribbon(data = plot_df_beta_forecast, aes(x = time, ymin = q025, ymax = q975), alpha = 0.3, fill = "darkolivegreen3") +
  geom_line(data = plot_df_f,aes(x = time, y = q50), color = "red") +
  geom_vline(xintercept = last_time, color = "blue", linetype = "dashed")+
  geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
  ylab(TeX("$R_t$")) +
  xlab("Time")
g_rt_comb 

ggarrange(g_trans_comb,g_rt_comb, ncol = 2, align = "v")








