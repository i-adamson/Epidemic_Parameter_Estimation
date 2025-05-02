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
library(tidyverse)

################################################################################
##  Data Set and Initial Values ################################################
################################################################################
##load data set:

Dengue_cases <- read_csv("~/Desktop/Dengue_cases.csv")
dengue_data <- Dengue_cases %>%
  pivot_wider(names_from = type_dengue, values_from = number)

##Determine your fixed values:
N <- 5630000
alpha = 0.01
gamma = 1/14


dengue_dat <- dengue_data[1:100,]
dengue_dat
##pull incidence values for each strain:
z1 <- dengue_dat$Dengue
z2 <- dengue_dat$DHF

z1 <- na.omit(z1)
z2 <- na.omit(z2)

l = length(z1)
#Determine the values for number of susceptibles and recovered:
R1 <- rep(0, l)
R2 <- rep(0, l)
for(i in 2:l){
  R1[i] <- R1[i-1] + alpha*z1[i] - gamma*R1[i]
  R2[i] <- R2[i-1] + alpha*z2[i] - gamma*R1[i]
}

S0 <- N - (z1+z2+R1+R2)

################################################################################
## PMCMC RUN TO ESTIMATE BETA1 AND BETA2 ###############################
################################################################################
y1 <- data.frame(value = z1) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)
y2 <- data.frame(value = z2) %>%
  mutate(time = seq(1, by = 1, length.out = n())) %>%
  dplyr::select(time, value)

ncores <- 8
minParticles <- max(ncores, 16)


################################################################################
multi_model_str <- "
model multi {
  obs y1
  obs y2

  state S
  state I1
  state I2
  state R1
  state R2
  state x1
  state x2

  state Z1
  state Z2

  input N
  param gamma
  param alpha
  param sigma // Noise driver
  param I10
  param I20
  param R10
  param R20
  param x10
  param x20
  param tau1
  param tau2

  sub parameter {
    gamma ~ truncated_gaussian(1.08, 0.075, lower=0.1, upper=5)
    alpha ~ truncated_gaussian(1.08, 0.075, lower=0.1, upper=5)
    sigma ~ truncated_gaussian(0.5, 0.2, lower=1e-6, upper=2)
    tau1 ~ truncated_gaussian(0.1, 0.05, lower=1e-6, upper=1)
    tau2 ~ truncated_gaussian(0.1, 0.05, lower=1e-6, upper=1)
    x10 ~ uniform(-5,2)
    x20 ~ uniform(-5,2)
    I10 ~ uniform(-9, -5)
    I20 ~ uniform(-9, -5)
    R10 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
    R20 ~ truncated_gaussian(0.15, 0.15, lower = 0, upper = 1)
  }

  sub initial {
    S <- N
    R1 <- R10*S
    R2 <- R20*S
    S <- S - R1 - R2

    I1 <- exp(I10 + log(S))
    I2 <- exp(I20 + log(S))
    S <- S - I1 - I2
    x1 <- x10
    x2 <- x20
    Z1 <- 0
    Z2 <- 0
  }

  sub transition(delta = 1) {
    S <- max(S, 1)
    I1 <- max(I1, 0)
    I2 <- max(I2, 0)
    R1 <- max(R1, 0)
    R2 <- max(R2, 0)
    
    Z1 <- ((t_now) % 1 == 0 ? 0 : Z1)
    Z2 <- ((t_now) % 1 == 0 ? 0 : Z2)
    
    noise e
    e ~ wiener()
    
    ode(alg = 'RK4(3)', h = 0.1, atoler = 1.0e-8, rtoler = 1.0e-8) {
      dx1/dt = sigma*e
      dx2/dt = sigma*e
      dS/dt = -exp(x1)*S*I1/N -exp(x2)*S*I2/N + R1/alpha + R2/alpha
      dI1/dt = exp(x1)*S*I1/N - I1/gamma
      dI2/dt = exp(x2)*S*I2/N - I2/gamma
      dR1/dt = I1/gamma - R1/alpha
      dR2/dt = I2/gamma - R2/alpha
      dZ1/dt = exp(x1)*S*I1/N
      dZ2/dt = exp(x2)*S*I2/N
    }
  }

  sub observation {
    y1 ~ poisson(max(Z1, 1e-10))
    y2 ~ poisson(max(Z2, 1e-10))
}

  sub proposal_parameter {
    sigma ~ gaussian(sigma, 0.01)
    gamma ~ gaussian(gamma, 0.01)
    alpha ~ gaussian(alpha, 0.01)
    x10 ~ gaussian(x10, 0.05)
    x20 ~ gaussian(x20, 0.05)
    I10 ~ gaussian(I10, 0.05)
    I20 ~ gaussian(I20, 0.05)
    R10 ~ gaussian(R10, 0.05)
    R20 ~ gaussian(R20, 0.05)
    tau1 ~ gaussian(tau1, 0.05)
    tau2 ~ gaussian(tau2, 0.05)
  }
}"


#perform pmcmc
multi_model <- bi_model(lines = stringi::stri_split_lines(multi_model_str)[[1]])

multi_bi_model <- libbi(multi_model)
input_lst <- list(N = 5630000)
end_time <- max(y1$time)
obs_lst <- list(y1 = y1 %>% dplyr::filter(time <= end_time),
                y2 = y2 %>% dplyr::filter(time <= end_time))

multi_bi <- sample(multi_model, target = "joint", end_time = end_time, input = input_lst, obs = obs_lst, nsamples = 1000, 
                    nparticles = minParticles, nthreads = ncores, proposal = 'prior', verbose = TRUE) %>% 
  adapt_particles(min = minParticles, max = minParticles*200)  %>%
  adapt_proposal(min = 0.05, max = 0.4, adapt = "both") %>%
  sample(nsamples = 5000, thin = 5) %>% # burn in 
  sample(nsamples = 5000, thin = 5)


#extract the samples
multi_bi_lst <- bi_read(multi_bi %>% sample_obs)

##take some time here to plot things to see if it produces correctly.
fitY <- multi_bi_lst$y1 %>% 
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup() %>%
  left_join(y1 %>% rename(Y = value))

g1 <- ggplot(data = fitY) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75), fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  geom_point(aes(x = time, y = Y), colour = "Red") +
  ylab("Incidence") +
  xlab("Time")

plot_df_beta1 <- multi_bi_lst$x1[multi_bi_lst$x1$time < 101,] %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()


g_beta1 <- ggplot(data = plot_df_beta1) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75),fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("Transmissibility ($\\beta_1(t)$)")) +
  xlab("Time")
g_beta1
 1)  %>% # burn in 
  sample(nsamples = 5, thin = 1)

plot_df_beta2 <- multi_bi_lst$x2 %>% mutate(value = exp(value)) %>%
  group_by(time) %>%
  mutate(
    q025 = quantile(value, 0.025),
    q25 = quantile(value, 0.25),
    q50 = quantile(value, 0.5),
    q75 = quantile(value, 0.75),
    q975 = quantile(value, 0.975)
  ) %>% ungroup()


g_beta2 <- ggplot(data = plot_df_beta2) +
  geom_ribbon(aes(x = time, ymin = q25, ymax = q75),fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = time, y = q50)) +
  ylab(TeX("$\\tilde{beta}_2(t)$")) +
  xlab("Time")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),)
g_beta2

#combined plot
plot_df_beta_comb <- bind_rows(
  plot_df_beta1 %>% mutate(beta_type = "Beta1"),
  plot_df_beta2 %>% mutate(beta_type = "Beta2")
)

g_beta_comb <- ggplot(data = plot_df_beta_comb) +
  geom_ribbon(data = plot_df_beta_comb, aes(x = time, ymin = q25, ymax = q75, fill = beta_type), alpha = 0.3) +
  geom_ribbon(data = plot_df_beta_comb, aes(x = time, ymin = q025, ymax = q975, fill = beta_type), alpha = 0.3) +
  geom_line(aes(x = time, y = q50, color = beta_type)) +
  ylab(TeX("Transmissibility ($\\beta(t)$)")) +
  xlab("Time") +
  scale_fill_manual(values = c("Beta1" = "royalblue3", "Beta2" = "cornflowerblue")) +
  scale_color_manual(values = c("Beta1" = "blue", "Beta2" = "darkblue")) +
  theme_minimal()


##Calculate Rt given the truncated model and plot

beta_1 <- data.frame(sample = c(multi_bi_lst$x1$np), time = c(multi_bi_lst$x1$time), beta_1 = c(exp(multi_bi_lst$x1$value)))
beta_2 <- data.frame(sample = c(multi_bi_lst$x2$np), time = c(multi_bi_lst$x2$time), beta_2 = c(exp(multi_bi_lst$x2$value)))
gamma_samp <- data.frame(sample = c(multi_bi_lst$gamma$np), gamma = c(multi_bi_lst$gamma$value))
St <- multi_bi_lst$S


Rt <- matrix(ncol = 2, nrow = 0)
colnames(Rt) <- c("Time", "Rt")

nsamp <- max(beta_1$sample)
tmax <- max(beta_1$time)

for(i in 0:nsamp){
  for(j in 1:tmax){
    gamma <- gamma_samp$gamma[i+1]
    beta1 <- beta_1[beta_1$sample == i & beta_1$time == j, "beta_1"]
    beta2 <- beta_2[beta_2$sample == i & beta_2$time == j, "beta_2"]
    lambda1 <- beta1/gamma
    lambda2 <- beta2/gamma
    Ro <- max(lambda1, lambda2)
    S <- St[St$np == i & St$time == j, "value"]
    m <- Ro*S/N
    Rt <- rbind(Rt, c(j, m))
  }
}
Rt <- as.data.frame(Rt)

plot_df_Rt <- Rt %>%
  group_by(Time) %>%
  mutate(
    q025 = quantile(Rt, 0.025),
    q25 = quantile(Rt, 0.25),
    q50 = quantile(Rt, 0.5),
    q75 = quantile(Rt, 0.75),
    q975 = quantile(Rt, 0.975)
  ) %>% ungroup()

g_Rt <- ggplot(data = plot_df_Rt) +
  geom_ribbon(aes(x = Time, ymin = q25, ymax = q75),fill="royalblue3", alpha = 0.3) +
  geom_ribbon(aes(x = Time, ymin = q025, ymax = q975),fill="cornflowerblue", alpha = 0.3) +
  geom_line(aes(x = Time, y = q50)) +
  geom_hline(yintercept=1, linetype = "dashed", color = "red")+
  ylab(TeX("$R_t$")) +
  xlab("Time")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),)
g_Rt

ggarrange(g1, g_beta1, g_beta2, g_Rt, ncol = 2, nrow = 2, align = "v")




