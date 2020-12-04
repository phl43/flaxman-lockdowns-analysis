library(tidyr)
library(dplyr)
library(rstan)
library(data.table)
library(lubridate)
library(gdata)
library(EnvStats)
library(matrixStats)
library(scales)
library(gridExtra)
library(ggpubr)
library(bayesplot)
library(cowplot)
library(ggplot2)

source("utils/geom-stepribbon.r")

make_model_fit_three_panel_plot <- function(data, covariates, filename_suffix) {
  covariates_long <- gather(covariates, key = "key", value = "value")
  covariates_long$x <- rep(0, length(covariates_long$key))
  un_dates <- unique(covariates_long$value)
  
  for (k in 1:length(un_dates)){
    idxs <- which(covariates_long$value == un_dates[k])
    max_val <- round(max(data$rt_ui)) + 0.3
    for (j in idxs){
      covariates_long$x[j] <- max_val
      max_val <- max_val - 0.3
    }
  }
  
  covariates_long$value <- as_date(covariates_long$value)
  
  data_cases_95 <- data.frame(data$date, data$predicted_cases_li, data$predicted_cases_ui)
  names(data_cases_95) <- c("date", "cases_min", "cases_max")
  data_cases_95$key <- rep("nintyfive", length(data_cases_95$date))
  data_cases_50 <- data.frame(data$date, data$predicted_cases_li2, data$predicted_cases_ui2)
  names(data_cases_50) <- c("date", "cases_min", "cases_max")
  data_cases_50$key <- rep("fifty", length(data_cases_50$date))
  data_cases <- rbind(data_cases_95, data_cases_50)
  levels(data_cases$key) <- c("ninetyfive", "fifty")
  
  p1 <- ggplot(data) +
    geom_bar(data = data, aes(x = date, y = reported_cases), fill = "coral4", stat='identity', alpha = 0.5) + 
    geom_ribbon(data = data_cases, aes(x = date, ymin = cases_min, ymax = cases_max, fill = key)) +
    xlab("") +
    ylab("Daily number of infections\n") +
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    scale_y_continuous(expand = c(0, 0), labels = comma) + 
    scale_fill_manual(
      name = "", labels = c("50%", "95%"),
      values = c(alpha("deepskyblue4", 0.55), alpha("deepskyblue4", 0.45))
      ) + 
    theme_pubr(base_family="sans") + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "None"
      ) +
    ggtitle("Simulation (model)") +
    guides(fill = guide_legend(ncol = 1))
  
  data_deaths_95 <- data.frame(data$date, data$predicted_deaths_li, data$predicted_deaths_ui)
  names(data_deaths_95) <- c("date", "death_min", "death_max")
  data_deaths_95$key <- rep("nintyfive", length(data_deaths_95$date))
  data_deaths_50 <- data.frame(data$date, data$predicted_deaths_li2, data$predicted_deaths_ui2)
  names(data_deaths_50) <- c("date", "death_min", "death_max")
  data_deaths_50$key <- rep("fifty", length(data_deaths_50$date))
  data_deaths <- rbind(data_deaths_95, data_deaths_50)
  levels(data_deaths$key) <- c("ninetyfive", "fifty")
  
  p2 <- ggplot(data, aes(x = date)) +
    geom_bar(data = data, aes(y = reported_deaths, fill = "reported"), fill = "coral4", stat='identity', alpha = 0.5) +
    geom_ribbon(
      data = data_deaths,
      aes(ymin = death_min, ymax = death_max, fill = key)) +
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) +
    scale_y_continuous(expand = c(0, 0), labels = comma) + 
    scale_fill_manual(
      name = "",
      labels = c("50%", "95%"),
      values = c(alpha("deepskyblue4", 0.55), alpha("deepskyblue4", 0.45))
      ) + 
    ylab("Daily number of deaths\n") + 
    xlab("") +
    theme_pubr(base_family="sans") + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "None"
      ) + 
    guides(fill = guide_legend(ncol = 1))
  
  plot_labels <- c("Complete lockdown", 
                   "Public events banned",
                   "School closure",
                   "Self isolation",
                   "Social distancing")
  
  data_rt_95 <- data.frame(data$date, data$rt_li, data$rt_ui)
  names(data_rt_95) <- c("date", "rt_min", "rt_max")
  data_rt_95$key <- rep("nintyfive", length(data_rt_95$date))
  data_rt_50 <- data.frame(data$date, data$rt_li2, data$rt_ui2)
  names(data_rt_50) <- c("date", "rt_min", "rt_max")
  data_rt_50$key <- rep("fifty", length(data_rt_50$date))
  data_rt <- rbind(data_rt_95, data_rt_50)
  levels(data_rt$key) <- c("ninetyfive", "fifth")
  
  p3 <- ggplot(data) +
    geom_stepribbon(data = data_rt, aes(x = date, ymin = rt_min, ymax = rt_max, group = key, fill = key)) +
    geom_hline(yintercept = 1, color = "black", size = 0.1) + 
    geom_segment(
      data = covariates_long,
      aes(x = value, y = 0, xend = value, yend = max(x)),
      linetype = "dashed",
      colour = "grey",
      alpha = 0.75
      ) +
    geom_point(
      data = covariates_long,
      aes(x = value, y = x, group = key, shape = key, col = key),
      size = 2
      ) +
    xlab("") +
    ylab(expression(R[t])) +
    scale_fill_manual(
      name = "",
      labels = c("50%", "95%"),
      values = c(alpha("seagreen", 0.75), alpha("seagreen", 0.5))
      ) + 
    scale_shape_manual(
      name = "Interventions",
      labels = plot_labels,
      values = c(21, 22, 23, 24, 25, 12)
      ) + 
    scale_colour_discrete(name = "Interventions", labels = plot_labels) + 
    scale_x_date(
      date_breaks = "weeks",
      labels = date_format("%e %b"),
      limits = c(data$date[1], data$date[length(data$date)])
      ) +
    scale_y_continuous(expand = c(0,0.1)) + 
    theme_pubr(base_family="sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position="right")
  
  p <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 2))
  save_plot(
    filename = paste0("figures/Simulation - epidemic - model -", filename_suffix, ".png"),
    p,
    base_width = 14
    )
}

make_simulation_data_three_panel_plot <- function(data, covariates, filename_suffix) {
  
  covariates_country <- covariates
  covariates_long <- gather(covariates_country, key = "key", value = "value")
  covariates_long$x <- rep(0, length(covariates_long$key))
  un_dates <- unique(covariates_long$value)
  
  for (k in 1:length(un_dates)){
    idxs <- which(covariates_long$value == un_dates[k])
    max_val <- max(data$Rt_adjusted) + 0.3
    for (j in idxs){
      covariates_long$x[j] <- max_val
      max_val <- max_val - 0.3
    }
  }
  
  covariates_long$value <- as_date(covariates_long$value)
  
  p1 <- ggplot(data, aes(x = date, y = cases)) +
    geom_bar(aes(y = reported_cases), fill = "coral4", stat = "identity", alpha = 0.5) + 
    geom_line(size = 1, color = "steelblue") +
    xlab("") +
    ylab("Daily number of infections\n") +
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) + 
    scale_y_continuous(expand = c(0, 0), labels = comma) + 
    theme_pubr(base_family = "sans") + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "None"
      ) +
    ggtitle("Simulation (data)") +
    guides(fill = guide_legend(ncol = 1))
  
  p2 <- ggplot(data, aes(x = date, y = deaths)) +
    geom_bar(aes(y = reported_deaths), fill = "coral4", stat = "identity", alpha = 0.5) +
    geom_line(size = 1, color = "steelblue") +
    scale_x_date(date_breaks = "weeks", labels = date_format("%e %b")) +
    scale_y_continuous(expand = c(0, 0), labels = comma) + 
    ylab("Daily number of deaths\n") + 
    xlab("") +
    theme_pubr(base_family = "sans") + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "None"
      ) + 
    guides(fill = guide_legend(ncol = 1))
  
  plot_labels <- c("Complete lockdown",
                   "Public events banned",
                   "School closure",
                   "Self isolation",
                   "Social distancing")
  
  p3 <- ggplot(data, aes(x = date, y = Rt_adjusted)) +
    geom_line(size = 1, color = "steelblue") +
    geom_hline(yintercept = 1, color = "black", size = 0.1) + 
    geom_segment(
      data = covariates_long,
      aes(x = value, y = 0, xend = value, yend = max(x)),
      linetype = "dashed", colour = "grey", alpha = 0.75
      ) +
    geom_point(
      data = covariates_long,
      aes(x = value, y = x, group = key, shape = key, col = key),
      size = 2
      ) +
    xlab("") +
    ylab(expression(R[t])) +
    scale_shape_manual(name = "Interventions", labels = plot_labels, values = c(21, 22, 23, 24, 25, 12)) + 
    scale_colour_discrete(name = "Interventions", labels = plot_labels) + 
    scale_x_date(
      date_breaks = "weeks",
      labels = date_format("%e %b"),
      limits = c(data$date[1], data$date[length(data$date)])
      ) +
    scale_y_continuous(expand = c(0, 0.1)) + 
    theme_pubr(base_family = "sans") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "right")
  
  p <- plot_grid(p1, p2, p3, ncol = 3, rel_widths = c(1, 1, 2))
  save_plot(
    filename = paste0("figures/Simulation - epidemic - actual -", filename_suffix, ".png"),
    p,
    base_width = 14
    )
}