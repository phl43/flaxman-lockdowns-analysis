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
library(extraDistr)
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

# increase to get correct number of days to simulate
forecast <- 0
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
  iter = 8000,
  warmup = 2000,
  chains = 4,
  thin = 1,
  control = list(adapt_delta = 0.99, max_treedepth = 10)
  )

out <- rstan::extract(fit)
prediction <- out$prediction
estimated.deaths <- out$E_deaths
estimated.deaths.cf <- out$E_deaths0

JOBID <- Sys.getenv("PBS_JOBID")
if(JOBID == "")
  JOBID <- as.character(abs(round(rnorm(1) * 1000000)))
print(sprintf("Jobid = %s",JOBID))

filename <- paste0("model-", JOBID)

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

print("Generating mu, rt plots")

mu = (as.matrix(out$mu))
colnames(mu) = countries
g = mcmc_intervals(mu, prob = .9)
ggsave(
  sprintf("figures/%s-mu.png", filename),
  g,
  width = 4,
  height = 6
  )

tmp = lapply(1:length(countries), function(i) (out$Rt_adj[, stan_data$N[i],i]))
Rt_adj = do.call(cbind, tmp)
colnames(Rt_adj) = countries
g = mcmc_intervals(Rt_adj, prob = .9)
ggsave(
  sprintf("figures/%s-final-rt.png", filename),
  g,
  width = 4,
  height = 6
  )

print("Generate 3-panel plots")
source("utils/plot-3-panel.r")
make_three_panel_plot(out, filename)

print("Making table")
source("utils/make-table.r")
make_table(out, filename)

# extract the posterior draws of the intervention effects

interventions <- c(
  "School Closure",
  "Self Isolation",
  "Public Events",
  "First Intervention",
  "Lockdown",
  "Social distancing"
)

alpha <- as.matrix(out$alpha)
alpha <- data.frame(alpha)
colnames(alpha) <- interventions

sum_alpha_nonsweden <- rowSums(alpha[,1:6])
sum_alpha_sweden <- rowSums(as.matrix(out$alpha)[,c(1,2,3,4,6)])

# extract the posterior draws for the country-specific effects

beta <- cbind(
  as.matrix(out$lockdown)[,1:9],
  as.matrix(out$last_intervention)[,10]
) %>%
  cbind(as.matrix(out$lockdown)[,11])
colnames(beta) <- countries
beta <- data.frame(beta)

# compute the sum of the effects, including the country-specific effect, for each country

results <- matrix(nrow = nrow(beta), ncol = ncol(beta))
colnames(results) <- countries
results <- data.frame(results)

for (i in 1:nrow(beta)) {
  for (j in 1:ncol(beta)) {
    if (j != 10) {
      results[i,j] <- sum_alpha_nonsweden[i] + beta[i,j]
    }
    else {
      results[i,j] <- sum_alpha_sweden[i] + beta[i,j]
    }
  }
}

# compute the mean of the posterior draws of the overall effect for each country and
# a 95% credible interval and prepare the data for display in a chart

overall_effects <- mcmc_intervals_data(
  results,
  prob_outer = 0.95,
  transformation = function(x) 1 - exp(-x),
  point_est = "mean"
)

levels(overall_effects$parameter) <- gsub("t(", "", levels(overall_effects$parameter), fixed = TRUE)
levels(overall_effects$parameter) <- gsub(")", "", levels(overall_effects$parameter), fixed = TRUE)
overall_effects$parameter <- factor(
  as.character(overall_effects$parameter),
  levels = sort(countries)
)

# compute the mean of the posterior draws of the country-specific effect for each country and
# a 95% credible interval and prepare the data for display in a chart

country_specific_effects <- mcmc_intervals_data(
  beta,
  prob_outer = 0.95,
  transformation = function(x) 1 - exp(-x),
  point_est = "mean"
)

levels(country_specific_effects$parameter) <- gsub("t(", "", levels(country_specific_effects$parameter), fixed = TRUE)
levels(country_specific_effects$parameter) <- gsub(")", "", levels(country_specific_effects$parameter), fixed = TRUE)
country_specific_effects$parameter <- factor(
  as.character(country_specific_effects$parameter),
  levels = sort(countries)
)

# compute the mean of the posterior draws of the effect for each intervention depending on whether
# the intervention was the first intervention or there was already other interventions in place
# when it started and a 95% credible interval and prepare the data for display in a chart

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

ggplot(overall_effects, aes(x = m, y = parameter)) +
  geom_errorbarh(aes(xmin = ll, xmax = hh, height = .2)) +
  geom_point(size = 2, color = "steelblue") +
  theme_bw() +
  ggtitle(expression(paste("Overall effects of interventions on ", R[t], " with 95% credibility interval"))) +
  xlab(expression(paste("Reduction in ", R[t]))) +
  ylab("") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
  ) +
  scale_y_discrete(limits = rev(levels(overall_effects$parameter))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggsave("figures/Overall effects of interventions.png", width = 12, height = 6)

ggplot(country_specific_effects, aes(x = m, y = parameter)) +
  geom_point(size = 2, color = "steelblue") +
  geom_errorbarh(aes(xmin = ll, xmax = hh, height = .2)) +
  theme_bw() +
  ggtitle(expression(paste("Country-specific effects of the last intervention on ", R[t], " with 95% credibility interval"))) +
  xlab(expression(paste("Reduction in ", R[t]))) +
  ylab("") +
  scale_x_continuous(
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_discrete(limits = rev(levels(country_specific_effects$parameter))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggsave("figures/Country-specific effects of the last intervention.png", width = 12, height = 6)

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

# numerically estimate the probability that Sweden's country-specific effect would be at least as large

set.seed(21)

n_simulations <- 1e7
gamma_draws <- rhnorm(n_simulations, 0.2)

beta_draws <- rep(0, n_simulations)

for (i in 1:n_simulations) {
  beta_draws[i] <- rnorm(1, sd = gamma_draws[i])
}

mean_beta_fit <- mean(beta$Sweden)

probability_sweden_effect <- mean(beta_draws >= mean_beta_fit)

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
    `Alternative counterfactual` = "Number of deaths averted relative to a counterfactual in which every country adopts Sweden's policies using Flaxman et al.'s estimate of how much those policies reduced transmission"
  ) %>%
  tab_header(
    title = "Number of deaths averted by lockdown with pooled model"
  ) %>%
  gtsave("figures/Number of deaths averted by lockdowns (pooled model).png")