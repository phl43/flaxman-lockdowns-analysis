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
library(tidyverse)
library(gt)

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

fit = sampling(
  m,
  data = stan_data,
  iter = 8000,
  warmup = 2000,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 10)
  )

out <- rstan::extract(fit)
prediction <- out$prediction
estimated.deaths <- out$E_deaths
estimated.deaths.cf <- out$E_deaths0

JOBID = Sys.getenv("PBS_JOBID")
if(JOBID == "")
  JOBID = as.character(abs(round(rnorm(1) * 1000000)))
print(sprintf("Jobid = %s",JOBID))

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
  file=paste0("results/model-", JOBID, "-stanfit.Rdata")
  )

filename <- paste0("model-", JOBID)

countries <- gsub("_", " ", countries)

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
Rt_adj <- do.call(cbind,tmp)
colnames(Rt_adj) <- countries
g <- mcmc_intervals(Rt_adj,prob = .9)
ggsave(
  sprintf("figures/%s-final-rt.png", filename),
  g,
  width = 4,
  height = 6
  )

source("utils/plot-3-panel.r")
make_three_panel_plot(out, filename)

source("utils/make-table.r")
make_table(out, filename)

interventions <- c(
  "School Closure",
  "Self Isolation",
  "Public Events",
  "First Intervention",
  "Lockdown",
  "Social distancing"
)

for (i in 1:length(countries)) {
  alpha <- as.matrix(out$alpha[,1:6,i])
  alpha <- data.frame(alpha)
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
  
  if (countries[i] == "Sweden") {
    individual_effects <- individual_effects %>%
      filter(parameter != "Lockdown" & parameter != "School Closure")
  } 
  
  # plot the effects
  
  ggplot(individual_effects, aes(x = m, y = parameter)) +
    geom_point(size = 2, position = position_dodge(-.5), aes(color = type)) +
    geom_errorbarh(aes(xmin = ll, xmax = hh, height = .2, color = type), position = position_dodge(-.5)) +
    theme_bw() +
    ggtitle(
      as.expression(bquote("Individual effects of interventions on " ~ R[t] ~ .(paste(" in ", countries[i], " with 95% credibility interval"))))
    ) +
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
    ggsave(paste0("figures/Individual effects of interventions in ", countries[i], " - ", filename, ".png"), width = 12, height = 6)
}

# create a table showing the number of deaths averted by country depending on the counterfactual used

averted_deaths <- read_csv(paste0("results/deaths-averted-", filename, ".csv")) %>%
  mutate(Country = gsub("_", " ", Country)) %>%
  select(
    Country,
    `Original counterfactual`,
    `Alternative counterfactual`
  )

averted_deaths %>%
  mutate(Country = factor(Country, levels = c(sort(Country[Country != "Total"]), "Total"))) %>%
  arrange(Country) %>%
  gt() %>%
  cols_align(
    align = "center",
    columns = vars(
      `Original counterfactual`,
      `Alternative counterfactual`
    )
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(
      columns = vars(
        Country,
        `Original counterfactual`,
        `Alternative counterfactual`
      )
    )
  ) %>%
  cols_label(
    Country = "Country",
    `Original counterfactual` = "Number of deaths averted relative to Flaxman et al.'s counterfactual",
    `Alternative counterfactual` = "Number of deaths averted relative to a counterfactual in which every country adopts Sweden's policies using Flaxman et al.'s estimate of how much those policies reduced transmission (non-pooled model)"
  ) %>%
  tab_header(
    title = "Number of deaths averted by lockdown with non-pooled model"
  ) %>%
  gtsave("figures/Number of deaths averted by lockdowns (non-pooled model).png")
