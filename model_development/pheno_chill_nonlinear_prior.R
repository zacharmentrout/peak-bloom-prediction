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
# Apply domain expertise to contain alpha
##################################################

# getting a sense for the range of bloom days offset by the pre-period
util$day_num_of_year('2024-03-20', base_date="2023-10-01") # 172
util$day_num_of_year('2024-06-20', base_date="2023-10-01") # 264

# higher confidence 
mean_Z_chill <- 3
sd_Z_chill <- 4 / 2.32 # 0 to 6
mean_Z_forcing <- 6 
sd_Z_forcing <- 4 / 2.32 # 3 to 9

# placeholders for now
mean_alpha <- 35
sd_alpha <- 5 / 2.32
mean_beta1 <- -0.5
sd_beta1 <- 1 / 2.32
mean_beta2 <- 0
sd_beta2 <- 0.01
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

#####################################################################
# plot temps and daily chill days / forcings
#####################################################################
sites_and_years <- mod_data$sites_and_years_df
accum <- T

local_temps <- matrix(NA, nrow(sites_and_years), 366)
for (i in 1:nrow(local_temps)) {
  temps <- mod_data$obs_temps[mod_data$temp_start_idxs[i]:
                                mod_data$temp_end_idxs[i]]
  
  local_temps[i,1:length(temps)] <- temps 
}

xlim <- c(1, 250)  # X-axis range


daily_chill_matrix <- extract_mod_prior$daily_chill_matrix[,,1:124]
daily_forcings_matrix <- extract_mod_prior$daily_forcings_matrix[,,125:xlim[2]]

if(accum) {
  for (j in 1:dim(daily_chill_matrix)[1]) {
    for (k in 1:dim(daily_chill_matrix)[2]) {
      daily_chill_matrix[j,k,] <- cumsum(daily_chill_matrix[j, k, ])
    }
  }
  
  for (j in 1:dim(daily_forcings_matrix)[1]) {
    for (k in 1:dim(daily_forcings_matrix)[2]) {
      daily_forcings_matrix[j,k,] <- cumsum(daily_forcings_matrix[j, k, ])
    }
  }
}

daily_chill_quantiles <- array(NA, dim=c(dim(daily_chill_matrix)[3], 9))
daily_forcings_quantiles <- array(NA, dim=c(dim(daily_forcings_matrix)[3], 9))

for (i in 1:dim(daily_chill_quantiles)[1]) {
  daily_chill_quantiles[i,] <- util$my_quantile(daily_chill_matrix[,,i],na.rm=T)
}

x0_chill <- daily_chill_quantiles[124,]

for (i in 1:dim(daily_forcings_quantiles)[1]) {
  daily_forcings_quantiles[i,] <- util$my_quantile(daily_forcings_matrix[,,i],na.rm=T)
}


daily_chill_plot <- daily_chill_quantiles
daily_forcings_plot <- daily_forcings_quantiles

ylim_right <- c(min(c(daily_forcings_plot, daily_chill_plot)), max(c(daily_forcings_plot, daily_chill_plot))) # c(0, 1)  # Range for right axis
tlim <- c(min(c(daily_forcings_plot, daily_chill_plot)), max(c(daily_forcings_plot, daily_chill_plot))) # c(0, 1)  # Range for right axis
tlim <- c(0, 150)

# Plot the main temperature data on the left axis
par(mar = c(5, 4, 4, 4) + 0.1, mfrow=c(1,1))  # Expand right margin for the second y-axis
plot("", main="Daily Temperature and Chill/Forcing Data",
     type="p", col="black", pch=16, cex=0.5, 
     xlim=xlim, xlab="Day", ylim=tlim, ylab="Cumulative Chill and Forcing Units",yaxt="n")

# Y-axis
axis(2, at = seq(0, 150, by=5))
axis(4, at = seq(0, 150, by=5))

# Overlay polygons for daily forcing quantiles (x = 125:366)
x_chill <- 1:124
x_forcing <- 125:xlim[2]

all_q <- c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
quantiles <- c(0.1, 0.25, 0.5, 0.75, 0.9)  # Selected quantiles


for (i in 1:floor(length(quantiles)/2)) {
  q1 <- quantiles[i]
  q2 <- quantiles[length(quantiles)-i+1]
  polygon(c(x_chill, rev(x_chill)), 
          #tlim[1] + c(daily_chill_plot[, i], rev(daily_chill_plot[, i + 1])) / diff(ylim_right) * diff(tlim),
          c(daily_chill_plot[, which(all_q == q1)], rev(daily_chill_plot[, which(all_q==q2)])),
          col=nom_colors[i], border=NA)
}

lines(x_chill,
      # tlim[1] + (daily_chill_plot[, which(quantiles == 0.5)] / diff(ylim_right)) * diff(tlim), 
      daily_chill_plot[, which(all_q == 0.5)],
      col=c_dark_highlight, lwd=2)  # Highlight the media

for (i in 1:floor(length(quantiles)/2)) {
  q1 <- quantiles[i]
  q2 <- quantiles[length(quantiles)-i+1]
  polygon(c(x_forcing, rev(x_forcing)), 
          # tlim[1] + c(daily_forcings_plot[, i], rev(daily_forcings_plot[, i + 1])) / diff(ylim_right) * diff(tlim),
          c(daily_forcings_plot[, which(all_q == q1)], rev(daily_forcings_plot[,which(all_q == q2)])),
          col=nom_colors[i], border=NA)
}

lines(x_forcing,
      # tlim[1] + (daily_forcings_plot[, which(quantiles == 0.5)] / diff(ylim_right)) * diff(tlim), 
      daily_forcings_plot[, which(all_q == 0.5)] ,
      col=c_dark_highlight, lwd=2)  # Highlight the median

util$plot_line_hist(mod_data$pheno_data_prior_df$doy, bin_delta = 5, bin_min=0, bin_max=250, add = T)
event_prior_quantiles <- quantile(mod_data$pheno_data_prior_df$doy,c(0.1,0.5, 0.9))
abline(v=event_prior_quantiles,lty='longdash')

median_forcings <- daily_forcings_plot[, which(all_q == 0.5)]

y_ref1 <- median_forcings[which(abs(event_prior_quantiles[1] - x_forcing) == min(abs(event_prior_quantiles[1]- x_forcing)))]
y_ref2 <- median_forcings[which(abs(event_prior_quantiles[2] - x_forcing) == min(abs(event_prior_quantiles[2]- x_forcing)))]
y_ref3 <- median_forcings[which(abs(event_prior_quantiles[3] - x_forcing) == min(abs(event_prior_quantiles[3]- x_forcing)))]

abline(h=c(y_ref1, y_ref2, y_ref3))

print(round(c(y_ref1, y_ref2, y_ref3)))

mean_alpha <- y_ref2
sd_alpha <- (y_ref3 - mean_alpha)/2.32

print(c(mean_alpha, y_ref3 - mean_alpha))


# conclusion: alpha: mean 35, sd 9.4/2.32


##################################################
# Contain beta1
##################################################

daily_chill_quantiles[124,]
# choosing an x neighboor based on chill day quantiles

x0_chill <- 38
dx_chill <- 35 # roughly covers most of the distribution

# range of forcings to motivate containment
q_low_forcings <- daily_forcings_plot[, which(all_q == 0.25)]
q_high_forcings <- daily_forcings_plot[, which(all_q == 0.75)]

# check where a broader range of forcings align with the end quantiles of unmodeled bloom days
y_ref1 <- median_forcings[which(abs(event_prior_quantiles[1] - x_forcing) == min(abs(event_prior_quantiles[1]- x_forcing)))]
y_ref2 <- median_forcings[which(abs(event_prior_quantiles[2] - x_forcing) == min(abs(event_prior_quantiles[2]- x_forcing)))]
y_ref3 <- median_forcings[which(abs(event_prior_quantiles[3] - x_forcing) == min(abs(event_prior_quantiles[3]- x_forcing)))]

dy_psi0 <- y_ref3 - y_ref2

# constrain slope parameter beta1 based on
dy_psi0 / dx_chill

# center beta1 near zero but bias negative based on prior findings of a negative relationship
# between chill days and the forcing threshold ("heating requirement")

mean_beta1 <- -0.5
sd_beta1 <- 0.21 / 2.32

# contain quadratic slope parameter beta2 based on
dy_psi0 / (dx_chill^2)

mean_beta2 <- 0
sd_beta2 <- 0.006 / 2.32





#####################################################################
# try simulating from prior predictive with new values
#####################################################################


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

mean(extract_mod_prior$pred_events == min(extract_mod_prior$pred_events))


#####################################################################
# test model
#####################################################################

mod_surv_fit <- stan(file=file.path("~/Documents/git/peak-bloom-prediction/stan_programs/pheno_chill_nonlinear.stan"), 
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
x_minus_x0 <- extract_mod$sum_chill_days_vector - mod_data$x0_chill_days
psi0 <- matrix(extract_mod$alpha, ncol = mod_data$N_site_years, nrow = length(extract_mod$alpha), byrow = F) +
  matrix(extract_mod$beta1, ncol = mod_data$N_site_years, nrow = length(extract_mod$beta1), byrow = F)*x_minus_x0 +
  matrix(extract_mod$beta2, ncol = mod_data$N_site_years, nrow = length(extract_mod$beta2), byrow = F)*(x_minus_x0^2)
for (k in 1:mod_data$N_site_years) {
  mod_surv_samples[[paste0('psi0[',k,']')]] <- matrix(psi0[,k], nrow=4, ncol=1500)
}
psi0_names <- grep('psi0', names(mod_surv_samples), value=TRUE)
util$plot_hist_quantiles(mod_surv_samples, val_name_prefix = "psi0")  
util$plot_conditional_mean_quantiles(mod_surv_samples, psi0_names,1:255)

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
util$plot_conditional_median_quantiles(mod_surv_samples, psi0_names,obs_xs = colMeans(x_minus_x0))

# chill day means
util$plot_conditional_median_quantiles(mod_surv_samples, pred_names,obs_xs = colMeans(x_minus_x0), baseline_values = mod_data$events)
util$plot_conditional_median_quantiles(mod_surv_samples, pred_names,obs_xs = colMeans(x_minus_x0), baseline_values = mod_data$events, residual = T)

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



