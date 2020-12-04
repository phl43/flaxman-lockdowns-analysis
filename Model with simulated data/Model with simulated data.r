library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(dplyr)
library(tidyr)
library(EnvStats)
library(optparse)
library(stringr)
library(bayesplot)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(abind)

simulate_data <- function(Rt, population, interventions_timing) {
  # length of simulation
  length <- length(Rt)
  
  # read the data on the serial interval
  serial_interval <- readRDS("data/serial-interval.rds")
  
  # pads serial interval with 0 if N2 is greater than the length of the serial interval array
  if (length(serial_interval$fit) < length) {
    pad_serial_interval <- data.frame(
      "X" = (length(serial.interval$fit)+1):N2,
      "fit" = rep(1e-17, max(length - length(serial.interval$fit), 0))
    )
    serial_interval = rbind(serial_interval, pad_serial_interval)
  }
  
  # create the list that will be passed to Stan
  stan_simulated_data<- list(
    M = 1,
    P = 6,
    N0 = 6,
    N = as.array(length),
    N2 = length,
    Rt_adjusted = array(rep(Rt[1], length), c(length, 1)),
    cases = array(rep(0, length), c(length, 1)),
    deaths = array(rep(0, length), c(length, 1)),
    f = array(rep(0, length), c(length, 1)),
    pop = as.array(population),
    SI = serial_interval$fit[1:length]
  )
  
  # start of the first intervention
  first_start <- min(
    interventions_timing$schools_start,
    interventions_timing$isolation_start,
    interventions_timing$events_start,
    interventions_timing$lockdown_start,
    interventions_timing$distancing_start
    )
  
  # vectors that say when each intervention starts and how long it stays in effect
  schools <- c(rep(0, interventions_timing$schools_start - 1), rep(1, length - interventions_timing$schools_start + 1))
  isolation <- c(rep(0, interventions_timing$isolation_start - 1), rep(1, length - interventions_timing$isolation_start + 1))
  events <- c(rep(0, interventions_timing$events_start - 1), rep(1, length - interventions_timing$events_start + 1))
  first <- c(rep(0, first_start - 1), rep(1, length - first_start + 1))
  lockdown <- c(rep(0, interventions_timing$lockdown_start - 1), rep(1, length - interventions_timing$lockdown_start + 1))
  distancing <- c(rep(0, interventions_timing$distancing_start - 1), rep(1, length - interventions_timing$distancing_start + 1))
  
  # create features matrix
  stan_simulated_data$X <- array(
    cbind(
      schools,
      isolation,
      events,
      first,
      lockdown,
      distancing
    ),
    dim = c(1, length, 6)
  )
  
  # create the infection to death distribution
  
  mean1 <- 5.1;
  cv1 <- 0.86;
  mean2 <- 17.8;
  cv2 <- 0.45
  
  # infection-to-onset distribution
  x1 <- rgammaAlt(1e6, mean1, cv1)
  
  # onset-to-death distribution
  x2 <- rgammaAlt(1e6, mean2, cv2)
  
  ecdf.saved <- ecdf(x1 + x2)
  
  IFR <- 0.01
  convolution <- function(u) (IFR * ecdf.saved(u))
  
  stan_simulated_data$f[1, 1] = (convolution(1.5) - convolution(0))
  for(i in 2:stan_simulated_data$N2) {
    stan_simulated_data$f[i, 1] = (convolution(i + 0.5) - convolution(i - 0.5)) 
  }
  
  # I need to reverse f and SI for the computations below
  
  f_rev <- rev(stan_simulated_data$f)
  SI_rev <- rev(stan_simulated_data$SI)
  
  # seeding infections to get the epidemic started
  
  cumulated_infections <- rep(0, stan_simulated_data$N2)
  
  tau <- rexp(1, 0.03)
  seed <- rexp(stan_simulated_data$N0, 1 / tau)
  
  stan_simulated_data$cases[1, 1] <- seed[1]
  for (i in 2:stan_simulated_data$N0) {
    stan_simulated_data$cases[i, 1] <- seed[i]
    cumulated_infections[i] <- cumulated_infections[i - 1] + stan_simulated_data$cases[i - 1]
  }
  
  # simulate the infections
  
  for (i in (stan_simulated_data$N0 + 1):stan_simulated_data$N2) {
    convolution <- as.vector(head(stan_simulated_data$cases, i - 1)) %*% tail(SI_rev, i - 1)
    cumulated_infections[i] <- cumulated_infections[i - 1] + stan_simulated_data$cases[i - 1, 1]
    stan_simulated_data$Rt_adjusted[i, 1] <- Rt[i] * (stan_simulated_data$pop - cumulated_infections[i - 1]) / stan_simulated_data$pop
    stan_simulated_data$cases[i, 1] <- stan_simulated_data$Rt_adjusted[i] * convolution
  }
  
  # simulate the deaths
  
  ifr_noise <- rnorm(1, 1, 0.1)
  
  stan_simulated_data$deaths[1, 1] <- 1e-15 * stan_simulated_data$cases[1, 1]
  
  for (i in 2:stan_simulated_data$N2) {
    stan_simulated_data$deaths[i, 1] <- ifr_noise * as.vector(head(stan_simulated_data$cases, i - 1)) %*% as.vector(tail(f_rev, i - 1))
  }
  
  # introduce underreporting and add some noise to simulated data (only deaths are used to fit the model)
  cases_underreporting <- 0.1
  cases_noise <- rnorm(stan_simulated_data$N2, 1, 0.3)
  deaths_underreporting <- 0.9
  deaths_noise <- rnorm(stan_simulated_data$N2, 1, 0.1)
  stan_simulated_data$reported_cases <- stan_simulated_data$cases * cases_underreporting * cases_noise
  stan_simulated_data$reported_deaths <- stan_simulated_data$deaths * deaths_underreporting * deaths_noise
  
  # round the number of cases and deaths because they are declared as integers in Stan
  stan_simulated_data$cases <- round(stan_simulated_data$cases)
  stan_simulated_data$deaths <- round(stan_simulated_data$deaths)
  stan_simulated_data$reported_cases <- round(stan_simulated_data$reported_cases)
  stan_simulated_data$reported_deaths <- round(stan_simulated_data$reported_deaths)
  
  # set the time at which the model start being estimated
  stan_simulated_data$EpidemicStart <- as.array(which(cumsum(stan_simulated_data$deaths) > 10)[1])
  
  stan_simulated_data
}

source("utils/plot-3-panel.r")

run_simulation_and_plot_results <- function(model, Rt, interventions_timing) {
  stan_simulated_data <- simulate_data(
    Rt,
    50e6,
    interventions_timing
  )
  
  # fit the model on the simulated data
  fit <- sampling(
    model,
    data = stan_simulated_data,
    iter = 4000,
    warmup = 2000,
    chains = 4,
    thin = 1,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  
  out <- rstan::extract(fit)
  
  dates <- list(seq(ymd("2020-02-01"), ymd("2020-02-01") + stan_simulated_data$N2 - 1, by = "day"))
  
  model_fit <- tibble(
    date = reduce(dates, c),
    reported_cases = stan_simulated_data$reported_cases,
    reported_deaths = stan_simulated_data$reported_deaths,
    predicted_cases = colMeans(out$prediction[,1:stan_simulated_data$N,1]),
    predicted_cases_li = colQuantiles(out$prediction[,1:stan_simulated_data$N,1], probs=.025),
    predicted_cases_ui = colQuantiles(out$prediction[,1:stan_simulated_data$N,1], probs=.975),
    predicted_cases_li2 = colQuantiles(out$prediction[,1:stan_simulated_data$N,1], probs=.25),
    predicted_cases_ui2 = colQuantiles(out$prediction[,1:stan_simulated_data$N,1], probs=.75),
    predicted_deaths = colMeans(out$E_deaths[,1:stan_simulated_data$N,1]),
    predicted_deaths_li = colQuantiles(out$E_deaths[,1:stan_simulated_data$N,1], probs=.025),
    predicted_deaths_ui = colQuantiles(out$E_deaths[,1:stan_simulated_data$N,1], probs=.975),
    predicted_deaths_li2 = colQuantiles(out$E_deaths[,1:stan_simulated_data$N,1], probs=.25),
    predicted_deaths_ui2 = colQuantiles(out$E_deaths[,1:stan_simulated_data$N,1], probs=.75),
    Rt = colMeans(out$Rt_adj[,1:stan_simulated_data$N,1]),
    rt_li = colQuantiles(out$Rt_adj[,1:stan_simulated_data$N,1],probs=.025),
    rt_ui = colQuantiles(out$Rt_adj[,1:stan_simulated_data$N,1],probs=.975),
    rt_li2 = colQuantiles(out$Rt_adj[,1:stan_simulated_data$N,1],probs=.25),
    rt_ui2 = colQuantiles(out$Rt_adj[,1:stan_simulated_data$N,1],probs=.75)
  )
  
  simulated_data <- tibble(
    date = reduce(dates, c),
    Rt_adjusted = stan_simulated_data$Rt_adjusted,
    cases = stan_simulated_data$cases,
    deaths = stan_simulated_data$deaths,
    reported_cases = stan_simulated_data$reported_cases,
    reported_deaths = stan_simulated_data$reported_deaths
  )
  
  covariates <- tibble(
    schools_universities = ymd("2020-02-01") + interventions_timing$schools_start - 1,
    self_isolating_if_ill = ymd("2020-02-01") + interventions_timing$isolation_start - 1,
    public_events = ymd("2020-02-01") + interventions_timing$events_start - 1,
    lockdown = ymd("2020-02-01") + interventions_timing$lockdown_start - 1,
    social_distancing_encouraged = ymd("2020-02-01") + interventions_timing$distancing_start - 1
  )
  
  make_model_fit_three_panel_plot(model_fit, covariates, "Simulation")
  make_simulation_data_three_panel_plot(simulated_data, covariates, "Simulation")
  
  interventions <- c(
    "School Closure",
    "Self Isolation",
    "Public Events",
    "First Intervention",
    "Lockdown",
    "Social distancing"
  )
  
  # compute the mean of the posterior draws of the effect for each intervention depending on whether
  # the intervention was the first intervention or there was already other interventions in place
  # when it started and a 95% credible interval and prepare the data for display in a chart
  
  alpha <- as.matrix(out$alpha[,,1])
  colnames(alpha) <- interventions
  
  individual_effects_first <- mcmc_intervals_data(
    alpha[,c(1,2,3,5,6)] + alpha[,4],
    prob_outer = 0.95,
    transformation = function(x) 1 - exp(-x),
    point_est = "mean"
  ) %>%
    mutate(type = "First intervention")
  
  individual_effects_later <- mcmc_intervals_data(
    alpha[,c(1,2,3,5,6)],
    prob_outer = 0.95,
    transformation = function(x) 1 - exp(-x),
    point_est = "mean"
  ) %>%
    mutate(type = "Later intervention")
  
  individual_effects <- bind_rows(
    individual_effects_first %>% mutate(parameter = as.character(parameter)),
    individual_effects_later %>% mutate(parameter = as.character(parameter))
  )
  
  individual_effects$parameter <- gsub("t(", "", individual_effects$parameter, fixed = TRUE)
  individual_effects$parameter <- gsub(")", "", individual_effects$parameter, fixed = TRUE)
  individual_effects$parameter <- factor(
    as.character(individual_effects$parameter),
    levels = sort(interventions)
  )
  
  # plot the effects
  
  ggplot(individual_effects, aes(x = m, y = parameter)) +
    geom_point(size = 2, position = position_dodge(-.5), aes(color = type)) +
    geom_errorbarh(aes(xmin = ll, xmax = hh, height = .2, color = type), position = position_dodge(-.5)) +
    theme_bw() +
    ggtitle(expression(paste("Individual effects of interventions on ", R[t], " with 95% credibility interval"))) +
    xlab(expression(paste("Reduction in ", R[t]))) +
    ylab("") +
    labs(color = "Timing") +
    scale_x_continuous(
      labels = percent_format(accuracy = 1)
    ) +
    scale_y_discrete(
      limits = rev(tail(levels(individual_effects$parameter), 5))
    ) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    ggsave("figures/Individual effects of interventions.png", width = 12, height = 6)
}

# seed the random number generator for reproducibility
set.seed(79)

# make sure the dates in the plots are in English
Sys.setenv("LANGUAGE"="En")
Sys.setlocale("LC_ALL", "en_US")

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
model <- stan_model(paste0("stan-models/model.stan"))

interventions_timing <- tibble(
  schools_start = 37,
  isolation_start = 33,
  events_start = 35,
  lockdown_start = 42,
  distancing_start = 30
)

# we assume that non-adjusted Rt goes down smoothly from 3.5 to 0.7 and stays there after that
R0 <- 3.5
R_final <- 0.7
Rt <- R0 - unlist(
  map(
    1:100,
    function(t) 1 / (1 + exp(-0.33 * (t - 40)))
  )
) * (R0 - R_final)

run_simulation_and_plot_results(model, Rt, interventions_timing)

# interventions_timing2 <- tibble(
#   schools_start = 36,
#   isolation_start = 32,
#   events_start = 34,
#   lockdown_start = 38,
#   distancing_start = 30
# )
# 
# run_simulation_and_plot_results(model, Rt1, interventions_timing2, "2")
# 
# interventions_timing3 <- tibble(
#   schools_start = 36,
#   isolation_start = 32,
#   events_start = 34,
#   lockdown_start = 45,
#   distancing_start = 30
# )
# 
# run_simulation_and_plot_results(model, Rt1, interventions_timing3, "3")
# 
# interventions_timing4 <- tibble(
#   schools_start = 47,
#   isolation_start = 43,
#   events_start = 45,
#   lockdown_start = 51,
#   distancing_start = 40
# )
# 
# # we assume that non-adjusted Rt starts at 3.5 and is reduced by approximately 33% by each intervention except
# # for the lockdown that doesn't affect it, so that it ends up at 0.7 after the last intervention took place
# Rt2 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt2[i] <- Rt2[i] *
#   exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing4$schools_start) *
#   exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing4$isolation_start) *
#   exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing4$events_start) *
#   exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing4$distancing_start) *
#   (1 - 0.2)^(i >= interventions_timing4$lockdown_start)
# }
# 
# run_simulation_and_plot_results(model, Rt2, interventions_timing4, "4")
# 
# # we assume that non-adjusted Rt starts at 3.5 and is reduced by approximately 33% by each intervention except
# # for the lockdown that doesn't affect it, so that it ends up at 0.7 after the last intervention took place
# Rt3 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt3[i] <- Rt3[i] -
#     0.65*(i >= interventions_timing4$schools_start) -
#     0.65*(i >= interventions_timing4$isolation_start) -
#     0.65*(i >= interventions_timing4$events_start) -
#     0.65*(i >= interventions_timing4$distancing_start) -
#     0.2*(i >= interventions_timing4$lockdown_start)
# }
# 
# run_simulation_and_plot_results(model, Rt3, interventions_timing4, "5")
# 
# # we assume that non-adjusted Rt starts at 3.5 and is reduced by approximately 33% by each intervention except
# # for the lockdown that doesn't affect it, so that it ends up at 0.7 after the last intervention took place
# Rt4 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt4[i] <- Rt4[i] -
#     0.65 / 3 * (i >= interventions_timing4$schools_start) -
#     0.65 / 3 * (i >= interventions_timing4$schools_start + 1) -
#     0.65 / 3 * (i >= interventions_timing4$schools_start + 2) -
#     0.65 / 3 *(i >= interventions_timing4$isolation_start) -
#     0.65 / 3 *(i >= interventions_timing4$isolation_start + 1) -
#     0.65 / 3 *(i >= interventions_timing4$isolation_start + 2) -
#     0.65 / 3 *(i >= interventions_timing4$events_start) -
#     0.65 / 3 *(i >= interventions_timing4$events_start + 1) -
#     0.65 / 3 *(i >= interventions_timing4$events_start + 2) -
#     0.65 / 3 * (i >= interventions_timing4$distancing_start) -
#     0.65 / 3 * (i >= interventions_timing4$distancing_start + 1) -
#     0.65 / 3 * (i >= interventions_timing4$distancing_start + 2) -
#     0.2 / 2 *(i >= interventions_timing4$lockdown_start) -
#     0.2 / 2 *(i >= interventions_timing4$lockdown_start + 1)
# }
# 
# run_simulation_and_plot_results(model, Rt4, interventions_timing4, "6")
# 
# interventions_timing5 <- tibble(
#   schools_start = 43,
#   isolation_start = 34,
#   events_start = 39,
#   lockdown_start = 45,
#   distancing_start = 30
# )
# 
# # we assume that non-adjusted Rt starts at 3.5 and is reduced by approximately 33% by each intervention except
# # for the lockdown that doesn't affect it, so that it ends up at 0.7 after the last intervention took place
# Rt5 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt5[i] <- Rt5[i] -
#     0.65 / 5 * (i >= interventions_timing5$schools_start) -
#     0.65 / 5 * (i >= interventions_timing5$schools_start + 1) -
#     0.65 / 5 * (i >= interventions_timing5$schools_start + 2) -
#     0.65 / 5 * (i >= interventions_timing5$schools_start + 3) -
#     0.65 / 5 * (i >= interventions_timing5$schools_start + 4) -
#     0.65 / 5 *(i >= interventions_timing5$isolation_start) -
#     0.65 / 5 *(i >= interventions_timing5$isolation_start + 1) -
#     0.65 / 5 *(i >= interventions_timing5$isolation_start + 2) -
#     0.65 / 5 *(i >= interventions_timing5$isolation_start + 3) -
#     0.65 / 5 *(i >= interventions_timing5$isolation_start + 4) -
#     0.65 / 5 *(i >= interventions_timing5$events_start) -
#     0.65 / 5 *(i >= interventions_timing5$events_start + 1) -
#     0.65 / 5 *(i >= interventions_timing5$events_start + 2) -
#     0.65 / 5 *(i >= interventions_timing5$events_start + 3) -
#     0.65 / 5 *(i >= interventions_timing5$events_start + 4) -
#     0.65 / 5 * (i >= interventions_timing5$distancing_start) -
#     0.65 / 5 * (i >= interventions_timing5$distancing_start + 1) -
#     0.65 / 5 * (i >= interventions_timing5$distancing_start + 2) -
#     0.65 / 5 * (i >= interventions_timing5$distancing_start + 3) -
#     0.65 / 5 * (i >= interventions_timing5$distancing_start + 4) -
#     0.2 / 2 *(i >= interventions_timing5$lockdown_start) -
#     0.2 / 2 *(i >= interventions_timing5$lockdown_start + 1)
# }
#   
# run_simulation_and_plot_results(model, Rt5, interventions_timing5, "7")
# 
# interventions_timing6 <- tibble(
#   schools_start = 41,
#   isolation_start = 42,
#   events_start = 40,
#   lockdown_start = 44,
#   distancing_start = 43
# )
# 
# # we assume that non-adjusted Rt starts at 3.5 and is reduced by approximately 33% by each intervention except
# # for the lockdown that doesn't affect it, so that it ends up at 0.7 after the last intervention took place
# Rt6 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt6[i] <- Rt6[i] *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$schools_start) *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$isolation_start) *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$events_start) *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$distancing_start) *
#     (1 - 0.2)^(i >= interventions_timing6$lockdown_start)
# }
# 
# run_simulation_and_plot_results(model, Rt6, interventions_timing6, "8")
# 
# # we assume that non-adjusted Rt starts at 3.5 and is reduced by approximately 33% by each intervention except
# # for the lockdown that doesn't affect it, so that it ends up at 0.7 after the last intervention took place
# Rt7 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt7[i] <- Rt7[i] *
#     exp(log(0.875 / 3.5))^(i >= interventions_timing6$lockdown_start)
# }
# 
# run_simulation_and_plot_results(model, Rt7, interventions_timing6, "9")
# 
# Rt8 <- rep(3.5, 100)
# for (i in 1:100) {
#   Rt8[i] <- Rt8[i] *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$schools_start) *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$isolation_start) *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$events_start) *
#     exp(log(0.875 / 3.5) / 4)^(i >= interventions_timing6$distancing_start)
# }
# 
# run_simulation_and_plot_results(model, Rt8, interventions_timing6, "10")
# 
# interventions_timing7 <- tibble(
#   schools_start = 51,
#   isolation_start = 47,
#   events_start = 49,
#   lockdown_start = 53,
#   distancing_start = 45
# )
# 
# run_simulation_and_plot_results(model, Rt1, interventions_timing7, "11")
# 
# interventions_timing8 <- tibble(
#   schools_start = 25,
#   isolation_start = 21,
#   events_start = 23,
#   lockdown_start = 30,
#   distancing_start = 18
# )
# 
# run_simulation_and_plot_results(model, Rt1, interventions_timing8, "12")
