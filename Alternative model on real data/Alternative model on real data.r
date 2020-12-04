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

source("utils/process-covariates.r")

# Read which countires to use
countries <- readRDS("data/regions.rds")
# Read deaths data for regions
d <- readRDS("data/COVID-19-up-to-date.rds")
# Read IFR and pop by country
ifr.by.country <- readRDS("data/popt-ifr.rds")

# Read interventions
interventions <- readRDS("data/interventions.rds")

forecast <- 0 # increase to get correct number of days to simulate
# Maximum number of days to simulate
N2 <- (max(d$DateRep) - min(d$DateRep) + 1 + forecast)[[1]]

processed_data <- process_covariates(
  countries = countries,
  interventions = interventions,
  d = d,
  ifr.by.country = ifr.by.country,
  N2 = N2
  )

stan_data <- processed_data$stan_data

dates <- processed_data$dates
deaths_by_country <- processed_data$deaths_by_country
reported_cases <- processed_data$reported_cases
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
m <- stan_model(paste0("stan-models/model.stan"))

fit <- sampling(
  m,
  data = stan_data,
  iter = 2000,
  warmup = 500,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
  )

out = rstan::extract(fit)
prediction = out$prediction
estimated.deaths = out$E_deaths
estimated.deaths.cf = out$E_deaths0

JOBID <- Sys.getenv("PBS_JOBID")
if(JOBID == "")
  JOBID <- as.character(abs(round(rnorm(1) * 1000000)))
print(sprintf("Jobid = %s", JOBID))

countries <- countries$Regions
save(
  fit,
  prediction,
  dates,
  reported_cases,
  deaths_by_country,
  countries,
  estimated.deaths,
  estimated.deaths.cf,
  stan_data,
  file = paste0("results/model-", JOBID, "-stanfit.Rdata")
  )

filename <- paste0("model-", JOBID)

print("Generating mu, rt plots")
mu <- as.matrix(out$mu)
colnames(mu) <- countries
g <- mcmc_intervals(mu, prob = .9)
ggsave(
  sprintf("figures/%s-mu.png", filename),
  g,
  width = 4,
  height = 6
  )
tmp <- lapply(
  1:length(countries),
  function(i) (out$Rt_adj[,stan_data$N[i],i])
  )
Rt_adj <- do.call(cbind, tmp)
colnames(Rt_adj) <- countries
g <- mcmc_intervals(Rt_adj, prob = .9)
ggsave(
  sprintf("figures/%s-final-rt.png", filename),
  g,
  width = 4,
  height = 6
  )

print("Generate 3-panel plots")
source("utils/plot-3-panel.r")
make_three_panel_plot(filename)