# This is copyrighted by Michael Betancourt and licensed under the new BSD (3-clause) license:
#  
#  https://opensource.org/licenses/BSD-3-Clause

############################################################
# Configure Graphics
############################################################

c_light <- c("#DCBCBC")
c_light_highlight <- c("#C79999")
c_mid <- c("#B97C7C")
c_mid_highlight <- c("#A25050")
c_dark <- c("#8F2727")
c_dark_highlight <- c("#7C0000")

c_light_teal <- c("#6B8E8E")
c_mid_teal <- c("#487575")
c_dark_teal <- c("#1D4F4F")

par(family="serif", las=1, bty="l", 
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

library(colormap)
library(scales)
nom_colors <- c("#DCBCBC", "#C79999", "#B97C7C", "#A25050", "#8F2727", "#7C0000")

############################################################
# Set up Stan
############################################################

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

util <- new.env()
wd <- '~/Documents/git/phenology/'

source(file.path(wd, 'r/functions_phenology.R'))
source(file.path('~/Documents/git/mcmc_visualization_tools/r/mcmc_analysis_tools_rstan.R'), local=util)
source(file.path('~/Documents/git/mcmc_visualization_tools/r/mcmc_visualization_tools.R'), local=util)

##################################################
# Prep data
##################################################

mod_data <- prep_data(wd)


##################################################
# original priors 
##################################################

# getting a sense for the range of bloom days offset by the pre-period
util$day_num_of_year('2024-03-20', base_date="2023-10-01") # 172
util$day_num_of_year('2024-06-20', base_date="2023-10-01") # 264

# higher confidence 
mean_Z_chill <- 3
sd_Z_chill <- 4 / 2.32 # 0 to 6
mean_Z_forcing <- 6 
sd_Z_forcing <- 4 / 2.32 # 3 to 9

mean_alpha <- 35
sd_alpha <- 7.45 / 2.32
mean_beta1 <- -0.5
sd_beta1 <- 0.21 / 2.32
mean_beta2 <- 0
sd_beta2 <- 0.006 / 2.32
sd_sigma <- 3

# center of covariate neighborhood 
x0_chill_days <- 38

mod_data$mean_alpha <- mean_alpha
mod_data$sd_alpha<- sd_alpha
mod_data$mean_beta1 <- mean_beta1
mod_data$sd_beta1 <- sd_beta1
mod_data$mean_beta2 <- mean_beta2
mod_data$sd_beta2 <- sd_beta2
mod_data$sd_sigma<- sd_sigma

mod_data$mean_Z_chill <- mean_Z_chill
mod_data$sd_Z_chill <- sd_Z_chill
mod_data$mean_Z_forcing <- mean_Z_forcing
mod_data$sd_Z_forcing <-sd_Z_forcing
mod_data$x0_chill_days <- x0_chill_days


simu_inits <- function(chain_id) {
  out <- list("alpha" = rnorm(1, mean_alpha, sd_alpha), 
              "beta1" = rnorm(1, mean_beta1, sd_beta1), 
              "beta2" = rnorm(1, mean_beta2, sd_beta2), 
              "Z_chill" = rnorm(1, mean_Z_chill, sd_Z_chill),
              "Z_forcing" = rnorm(1, mean_Z_forcing, sd_Z_forcing),
              "sigma" = abs(rnorm(1, 0, sd_sigma)))
  print(paste0('Chain ', chain_id, '...'))
  print(out)
  return(out)
}

mod_surv_fit_prior <- stan(file=file.path("~/Documents/git/peak-bloom-prediction/stan_programs/pheno_chill_nonlinear_prior.stan"), 
                           data=mod_data, seed=8438339, algorithm="Fixed_param",
                           iter=1000, chains=1, warmup=0,init=simu_inits
)

extract_mod_prior <- extract(mod_surv_fit_prior)

# The diagnostics for all parameter expectands are clear.
mod_surv_samples_prior <- util$extract_expectand_vals(mod_surv_fit_prior)
names <- c('alpha', 'beta1', 'beta2', 'sigma',  "Z_forcing", "Z_chill") #c(daily_forcings_plot, daily_chill_plot)
base_samples <- util$filter_expectands(mod_surv_samples_prior, names)
util$check_all_expectand_diagnostics(base_samples)



# Prior predictive
par(mfrow=c(1, 1))
pred_names <- grep('pred_events', names(mod_surv_samples_prior), value=TRUE)
util$plot_hist_quantiles(samples = mod_surv_samples_prior, val_name_prefix = "pred_events",bin_delta = 5, baseline_values = mod_data$pheno_data_prior_df$doy, xlab = "Phenology Event")

util$my_quantile(extract_mod_prior$pred_events)
util$my_quantile(mod_data$pheno_data_prior_df$doy)

#####################################################################
# add forecasts
#####################################################################
mean_delta <- 0
sd_delta <- 2/2.32 # 10 / 2.32
sd_tau <- 2/2.32#7 / 2.32

mod_data$mean_delta <- mean_delta
mod_data$sd_delta <- sd_delta
mod_data$sd_tau <- sd_tau

simu_inits <- function(chain_id) {
  out <- list("alpha" = rnorm(1, mean_alpha, sd_alpha), 
              "beta1" = rnorm(1, mean_beta1, sd_beta1), 
              "beta2" = rnorm(1, mean_beta2, sd_beta2), 
              "Z_chill" = rnorm(1, mean_Z_chill, sd_Z_chill),
              "Z_forcing" = rnorm(1, mean_Z_forcing, sd_Z_forcing),
              "sigma" = abs(rnorm(1, 0, sd_sigma)),
              "delta" = rnorm(1, mean_delta, sd_delta),
              "tau" = abs(rnorm(1, 0, sd_tau))
              )
  print(paste0('Chain ', chain_id, '...'))
  print(out)
  return(out)
}

mod_surv_fit_w_forecasts_prior <- stan(file=file.path("~/Documents/git/peak-bloom-prediction/stan_programs/pheno_chill_nonlinear_with_forecasts_prior.stan"), 
                           data=mod_data, seed=8438339, algorithm="Fixed_param",
                           iter=1000, chains=1, warmup=0,init=simu_inits
)


# The diagnostics for all parameter expectands are clear.
mod_surv_samples_prior_forecast <- util$extract_expectand_vals(mod_surv_fit_w_forecasts_prior)
names <- c('alpha', 'beta1', 'beta2', 'sigma',  "Z_forcing", "Z_chill", "delta", "tau") #c(daily_forcings_plot, daily_chill_plot)
base_samples <- util$filter_expectands(mod_surv_samples_prior_forecast, names)
util$check_all_expectand_diagnostics(base_samples)

# Prior predictive
par(mfrow=c(1, 1))
pred_names <- grep('pred_events', names(mod_surv_samples_prior_forecast), value=TRUE)
pred_forecast_names <- grep('pred_events_forecast', names(mod_surv_samples_prior_forecast), value=TRUE)
util$plot_hist_quantiles(samples = mod_surv_samples_prior_forecast, val_name_prefix = "pred_events",bin_delta = 5, baseline_values = mod_data$pheno_data_prior_df$doy, xlab = "Phenology Event")
util$plot_hist_quantiles(samples = mod_surv_samples_prior_forecast, val_name_prefix = "pred_events_forecast",bin_delta = 5, baseline_values = NULL, xlab = "Phenology Event")


util$my_quantile(extract_mod_prior$pred_events)
util$my_quantile(mod_data$pheno_data_prior_df$doy)

pred_events_forecast <- as.data.frame(extract(mod_surv_fit_w_forecasts_prior, pred_forecast_names))
hist(pred_events_forecast[,5])



#####################################################################
# test model
#####################################################################

mod_surv_fit <- stan(file=file.path("~/Documents/git/peak-bloom-prediction/stan_programs/pheno_chill_nonlinear_with_forecasts.stan"), 
                     data=mod_data, init=simu_inits, seed=8438339,
                     warmup=1500, iter=3000, control = list(adapt_delta = 0.8))#, refresh=0)

extract_mod <- extract(mod_surv_fit)
# All Hamiltonian Monte Carlo diagnostics are clear.
diagnostics <- util$extract_hmc_diagnostics(mod_surv_fit)
util$check_all_hmc_diagnostics(diagnostics)

# The diagnostics for all parameter expectands are clear.
mod_surv_samples <- util$extract_expectand_vals(mod_surv_fit)
names <- c('alpha', 'beta1', 'beta2', 'sigma','Z_forcing', 'Z_chill') #c(daily_forcings_plot, daily_chill_plot)
base_samples <- util$filter_expectands(mod_surv_samples, names)
util$check_all_expectand_diagnostics(base_samples)


par(mfrow=c(1,1))
#x_minus_x0 <- extract_mod$sum_chill_days_vector - mod_data$x0_chill_days
x_minus_x0 <- seq(0,80) - mod_data$x0_chill_days
psi0 <- matrix(rep(extract_mod$alpha, length(x_minus_x0)), nrow =length(extract_mod$alpha), ncol = length(x_minus_x0)) +
  matrix(extract_mod$beta1)%*%x_minus_x0 +
  matrix(extract_mod$beta2)%*%(x_minus_x0^2)
for (k in 1:length(x_minus_x0)) {
  mod_surv_samples[[paste0('psi0[',k,']')]] <- matrix(psi0[,k], nrow=4, ncol=1500)
}
psi0_names <- grep('psi0', names(mod_surv_samples), value=TRUE)[1:length(x_minus_x0)]
util$plot_hist_quantiles(mod_surv_samples, val_name_prefix = "psi0")  
util$plot_conditional_mean_quantiles(mod_surv_samples, psi0_names,x_minus_x0)

# Posterior predictive
par(mfrow=c(1, 1))
pred_names <- grep('pred_events', names(mod_surv_samples), value=TRUE)
util$plot_hist_quantiles(samples = mod_surv_samples, val_name_prefix = "pred_events",bin_delta = 5, baseline_values = mod_data$events, xlab = "Phenology Event")

util$my_quantile(extract_mod$pred_events)
util$my_quantile(mod_data$events)
mean(extract_mod_prior$pred_events == min(extract_mod_prior$pred_events))


mean_temps <- rep(NA, mod_data$N_site_years)
mean_chilling <- rep(NA, mod_data$N_site_years)
mean_forcing <- rep(NA, mod_data$N_site_years)

for(i in 1:(mod_data$N_site_years)) {
  temps <- mod_data$obs_temps[mod_data$temp_start_idxs[i]:(mod_data$temp_start_idxs[i]+mod_data$N_days[i]-1)]
  mean_temps[i] <- mean(temps)
  mean_chilling[i] <- mean(temps[70:mod_data$precursor_events[i]])
  mean_forcing[i] <- mean(temps[(mod_data$precursor_events[i]+1):200])#mod_data$N_days[i]])
}

par(mfrow=c(1,2))
# events vs. events
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mod_data$events, baseline_values = mod_data$events, bin_delta = 3)
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mod_data$events, residual = T, baseline_values = mod_data$events, bin_delta=3)
# full-year
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_temps, baseline_values = mod_data$events, bin_delta = 0.25)
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_temps, residual = T, baseline_values = mod_data$events, bin_delta = 0.25)
# chilling
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_chilling, baseline_values = mod_data$events, bin_delta = 0.5)
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_chilling, residual = T, baseline_values = mod_data$events, bin_delta = 0.5)
# forcing
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_forcing, baseline_values = mod_data$events, bin_delta = 0.5)
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_forcing, residual = T, baseline_values = mod_data$events, bin_delta = 0.5)
# years
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mod_data$sites_and_years_df[,"year"], baseline_values = mod_data$events, bin_delta = 4)
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mod_data$sites_and_years_df[,"year"], residual = T, baseline_values = mod_data$events, bin_delta = 4)
# chilling/forcing combined
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_forcing - mean_chilling, baseline_values = mod_data$events)
util$plot_conditional_median_quantiles(mod_surv_samples, names=pred_names, obs_xs = mean_forcing - mean_chilling, residual = T, baseline_values = mod_data$events)

# psi0
util$plot_conditional_median_quantiles(mod_surv_samples, psi0_names,obs_xs = x_minus_x0)


# by location
locations <- unique(mod_data$sites_and_years_df$location)
for (l in locations) {
  bin_min <- NULL
  bin_max <- NULL
  idx <- mod_data$sites_and_years_df$location == l
  if (sum(idx) == 1) {
    bin_min <- mod_data$events[idx]-1
    bin_max <- mod_data$events[idx]+1
  }
  util$plot_conditional_median_quantiles(mod_surv_samples, pred_names[idx],obs_xs = mod_data$events[idx], baseline_values = mod_data$events[idx], main = l, bin_min=bin_min, bin_max=bin_max)
  util$plot_conditional_median_quantiles(mod_surv_samples, pred_names[idx],obs_xs = mod_data$events[idx], baseline_values = mod_data$events[idx], residual = T, main=l, bin_min=bin_min, bin_max=bin_max)
}

locations <- unique(mod_data$sites_and_years_df$location)
for (l in locations) {
  bin_min <- NULL
  bin_max <- NULL
  idx <- mod_data$sites_and_years_df$location == l
  if (sum(idx) == 1) {
    bin_min <- colMeans(x_minus_x0)[idx]-1
    bin_max <- colMeans(x_minus_x0)[idx]+1
  }
  util$plot_conditional_median_quantiles(mod_surv_samples, pred_names[idx],obs_xs = colMeans(x_minus_x0)[idx], baseline_values = mod_data$events[idx], main = l, bin_min=bin_min, bin_max=bin_max)
  util$plot_conditional_median_quantiles(mod_surv_samples, pred_names[idx],obs_xs = colMeans(x_minus_x0)[idx], baseline_values = mod_data$events[idx], residual = T, main=l, bin_min=bin_min, bin_max=bin_max)
}

locations <- unique(mod_data$sites_and_years_df$location)
for (l in locations) {
  bin_min <- NULL
  bin_max <- NULL
  idx <- mod_data$sites_and_years_df$location == l
  if (sum(idx) == 1) {
    bin_min <- mod_data$sites_and_years_df[idx, "year"]-1
    bin_max <- mod_data$sites_and_years_df[idx, "year"]+1
  }
  util$plot_conditional_median_quantiles(mod_surv_samples, pred_names[idx],obs_xs = mod_data$sites_and_years_df[idx, "year"], baseline_values = mod_data$events[idx], main = l, bin_min=bin_min, bin_max=bin_max)
  util$plot_conditional_median_quantiles(mod_surv_samples, pred_names[idx],obs_xs = mod_data$sites_and_years_df[idx, "year"], baseline_values = mod_data$events[idx], residual = T, main=l, bin_min=bin_min, bin_max=bin_max)
}

# by location
locations <- unique(mod_data$sites_and_years_df$location)
for (l in locations) {
  bin_min <- NULL
  bin_max <- NULL
  idx <- mod_data$sites_and_years_df$location == l
  if (sum(idx) == 1) {
    bin_min <- mod_data$events[idx]-1
    bin_max <- mod_data$events[idx]+1
  }
  mod_surv_samples_sub <- lapply(pred_names[which(idx)], function(name) mod_surv_samples[[name]])
  names(mod_surv_samples_sub) <- pred_names[which(idx)]
  util$plot_hist_quantiles(samples =  mod_surv_samples_sub, val_name_prefix = "pred_events",bin_delta = 5, baseline_values = mod_data$events[idx], xlab = "Phenology Event", main=l)
}



# In contrast to the standard survival model, which lead to inferences of
# narrow forcing functions, the modified survival model suggests that wider
# forcing functions are in fact a better fit to the observed data.
par(mfrow=c(2, 2))

util$plot_expectand_pushforward(mod_surv_samples[["Z_chill"]], 25,
                                display_name="Z_chill", flim=c(-3, 10))
xs <- seq(-3, 10, 0.01)
ys <- dnorm(xs,mod_data$mean_Z_chill , mod_data$sd_Z_chill)
lines(xs, ys, lwd=2, col=c_light_teal)


util$plot_expectand_pushforward(mod_surv_samples[["Z_forcing"]], 25,
                                display_name="Z_forcing", flim=c(-3, 10))
xs <- seq(-3, 10, 0.01)
ys <- dnorm(xs, mod_data$mean_Z_forcing , mod_data$sd_Z_forcing)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["alpha"]], 25,
                                display_name="alpha", flim=c(10, 50))
xs <- seq(10, 50, 0.01)
ys <- dnorm(xs, mean_alpha, sd_alpha)
lines(xs, ys, lwd=2, col=c_light_teal)
abline(v=alpha_true)

util$plot_expectand_pushforward(mod_surv_samples[["beta1"]], 25,
                                display_name="beta1", flim=c(-2, 1))
xs <- seq(-2, 1, 0.01)
ys <- dnorm(xs, mean_beta1, sd_beta1)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["beta2"]], 25,
                                display_name="beta2", flim=c(-0.01, 0.02))
xs <- seq(-0.03, 0.03, 0.0001)
ys <- dnorm(xs, mean_beta2, sd_beta2)
lines(xs, ys, lwd=2, col=c_light_teal)


util$plot_expectand_pushforward(mod_surv_samples[["sigma"]], 25,
                                display_name="sigma", flim=c(0, 15))
xs <- seq(0, 15, 0.1)
ys <- dnorm(xs, 0, sd_sigma)
lines(xs, ys, lwd=2, col=c_light_teal)
abline(v=sigma_true)


#####################################################################
# model the temp forecasts along with it
#####################################################################









