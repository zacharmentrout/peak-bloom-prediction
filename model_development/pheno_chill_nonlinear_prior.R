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
# Set chill and forcing thresholds
##################################################

day_num_of_year('2024-03-20', base_date="2023-10-01") # 172
day_num_of_year('2024-06-20', base_date="2023-10-01") # 264

# higher confidence 
mean_Z_chill <- 3
sd_Z_chill <- 4 / 2.32 # -1 to 7
mean_Z_forcing <- 6 
sd_Z_forcing <- 4 / 2.32 # 2 to 10

# placeholders for now
mean_alpha <- 35
sd_alpha <- 5 / 2.32
mean_beta <- -0.5
sd_beta <- 1 / 2.32
sd_sigma <- 3


mod_data$mean_alpha <- mean_alpha
mod_data$sd_alpha<- sd_alpha
mod_data$mean_beta <- mean_beta
mod_data$sd_beta <- sd_beta
mod_data$sd_sigma<- sd_sigma

mod_data$mean_Z_chill <- mean_Z_chill
mod_data$sd_Z_chill <- sd_Z_chill
mod_data$mean_Z_forcing <- mean_Z_forcing
mod_data$sd_Z_forcing <-sd_Z_forcing


simu_inits <- function(chain_id) {
  out <- list("alpha" = rnorm(1, mean_alpha, sd_alpha), 
              "beta" = rnorm(1, mean_beta, sd_beta), 
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
names <- c('alpha', 'beta', 'sigma',  "Z_forcing", "Z_chill") #c(daily_forcings_plot, daily_chill_plot)
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
event_prior_quantiles <- quantile(mod_data$pheno_data_prior_df$doy,c(0.05,0.5, 0.95))
abline(v=event_prior_quantiles,lty='longdash')

median_forcings <- daily_forcings_plot[, which(all_q == 0.5)]

y_ref1 <- median_forcings[which(abs(event_prior_quantiles[1] - x_forcing) == min(abs(event_prior_quantiles[1]- x_forcing)))]
y_ref2 <- median_forcings[which(abs(event_prior_quantiles[2] - x_forcing) == min(abs(event_prior_quantiles[2]- x_forcing)))]
y_ref3 <- median_forcings[which(abs(event_prior_quantiles[3] - x_forcing) == min(abs(event_prior_quantiles[3]- x_forcing)))]

abline(h=c(y_ref1, y_ref2, y_ref3))

print(round(c(y_ref1, y_ref2, y_ref3)))

test_mean_alpha <- (y_ref1 + y_ref3)/2
test_sd_alpha <- (y_ref3 - test_mean_alpha)/2.32

print(c(test_mean_alpha, test_sd_alpha))

##################################################
# Modified survival model
##################################################

set.seed(49384838)

day_num_of_year('2024-03-20', base_date="2023-10-01") # 172
day_num_of_year('2024-06-20', base_date="2023-10-01") # 264

mean_alpha <- 35
sd_alpha <- 5 / 2.32
mean_beta <- -0.5
sd_beta <- 1 / 2.32
mean_Z <- 5
sd_Z <- 5 / 2.32
sd_sigma <- 3

mean_Z_chill <- 3
sd_Z_chill <- 4 / 2.32 # -1 to 7
mean_Z_forcing <- 6 
sd_Z_forcing <- 4 / 2.32 # 2 to 10


mod_data$mean_alpha <- mean_alpha
mod_data$sd_alpha<- sd_alpha
mod_data$mean_beta <- mean_beta
mod_data$sd_beta <- sd_beta
mod_data$mean_Z <- mean_Z
mod_data$sd_Z <- sd_Z
mod_data$sd_sigma<- sd_sigma

mod_data$mean_Z_chill <- mean_Z_chill
mod_data$sd_Z_chill <- sd_Z_chill
mod_data$mean_Z_forcing <- mean_Z_forcing
mod_data$sd_Z_forcing <-sd_Z_forcing
mod_data$Z_chill_fixed <- 5


simu_inits <- function(chain_id) {
  out <- list("alpha" = rnorm(1, mean_alpha, sd_alpha), 
              "beta" = rnorm(1, mean_beta, sd_beta), 
              #"Z" = rnorm(1, mean_Z, sd_Z), 
              "Z_chill" = rnorm(1, mean_Z_chill, sd_Z_chill),
              "Z_forcing" = rnorm(1, mean_Z_forcing, sd_Z_forcing),
              "sigma" = abs(rnorm(1, 0, sd_sigma)))
  print(paste0('Chain ', chain_id, '...'))
  print(out)
  return(out)
}

mod_surv_fit_prior <- stan(file=file.path(wd, "stan_files/pheno_modified_survival_chill_nonlinear_prior.stan"), 
                       data=mod_data, seed=8438339, algorithm="Fixed_param",
                       iter=1000, chains=1, warmup=0,init=simu_inits
)
mod_surv_fit <- mod_surv_fit_prior
extract_mod <- extract(mod_surv_fit)

(extract_mod$sum_chill_days_vector)
util$plot_hist_quantiles(mod_surv_samples,"sum_chill_days_vector")
my_quantile(extract_mod$sum_chill_days_vector)
mean(extract_mod$sum_chill_days_vector)


#####
## fit to prior predictive to see if there's a model issue
#####
prior_samples <- util$extract_expectands(mod_surv_fit_prior)
draw_num <- 28
event_draw <- sapply(1:mod_data$N_site_years, function(k) { prior_samples[[paste0('pred_events[', k,']')]][draw_num]})
# Z_chill_true <- prior_samples[['Z_chill']][draw_num]
alpha_true <- prior_samples[['alpha']][draw_num]
beta_true <- prior_samples[['beta']][draw_num]
Z_forcing_true <- prior_samples[['Z_forcing']][draw_num]
sigma_true <- prior_samples[['sigma']][draw_num]

mod_data$events <- event_draw

##########

# 1 to 200ish range?
mod_surv_fit <- stan(file=file.path(wd,
                                    #"stan_files/pheno_modified_survival_chill_one_thresh.stan"), 
                                    
                                    "stan_files/pheno_modified_survival_simplified_with_chill.stan"), 
                     data=mod_data, init=simu_inits, seed=8438339,
                     warmup=1000, iter=2000)#, refresh=0)

extract_mod <- extract(mod_surv_fit)


# All Hamiltonian Monte Carlo diagnostics are clear.
diagnostics <- util$extract_hmc_diagnostics(mod_surv_fit)
util$check_all_hmc_diagnostics(diagnostics)


# The diagnostics for all parameter expectands are clear.
mod_surv_samples <- util$extract_expectand_vals(mod_surv_fit)
names <- c('alpha', 'beta', 'sigma',  "Z_forcing", "Z_chill") #c(daily_forcings_plot, daily_chill_plot)
base_samples <- util$filter_expectands(mod_surv_samples, names)
util$check_all_expectand_diagnostics(base_samples)


daily_chill_matrix <- extract_mod$daily_chill_matrix[,,1:124]
daily_forcings_matrix <- extract_mod$daily_forcings_matrix[,,125:366]

accum <- T
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

daily_chill_quantiles <- array(NA, dim=c(dim(daily_chill_matrix)[2:3], 9))
daily_forcings_quantiles <- array(NA, dim=c(dim(daily_forcings_matrix)[2:3], 9))


for (i in 1:dim(daily_forcings_quantiles)[1]) {
  for (j in 1:dim(daily_forcings_quantiles)[2]) {
    daily_forcings_quantiles[i,j,] <- my_quantile(daily_forcings_matrix[,i,j],na.rm=T)
  }
}

for (i in 1:dim(daily_chill_quantiles)[1]) {
  for (j in 1:dim(daily_chill_quantiles)[2]) {
    daily_chill_quantiles[i,j,] <- my_quantile(daily_chill_matrix[,i,j],na.rm=T)
  }
}



#####################################################################
# plot temps and daily chill days / forcings
#####################################################################
sites_and_years <- mod_data$sites_and_years_df
year <- 2022
location <- 'kyoto'
index <- which(sites_and_years$location == location & sites_and_years$year == year)


local_temps <- mod_data$obs_temps[mod_data$temp_start_idxs[index]:
                     mod_data$temp_end_idxs[index]]
xlim <- c(1, length(local_temps))  # X-axis range

# Define plot limits
tlim <- range(local_temps, na.rm = TRUE)  # Temperature limits for left axis

daily_chill_plot <- daily_chill_quantiles[index,1:124,]
daily_forcings_plot <- daily_forcings_quantiles[index, 1:(xlim[2]-124),]

ylim_right <- c(min(c(daily_forcings_plot, daily_chill_plot)), max(c(daily_forcings_plot, daily_chill_plot))) # c(0, 1)  # Range for right axis


# Plot the main temperature data on the left axis
par(mar = c(5, 4, 4, 4) + 0.1, mfrow=c(1,1))  # Expand right margin for the second y-axis
plot("", main="Daily Temperature and Chill/Forcing Data",
     type="p", col="black", pch=16, cex=0.5, 
     xlim=xlim, xlab="Day", ylim=tlim, ylab="Recorded Temperature (C)")


# Add the right axis directly (no need to rescale)
axis(4, at=seq(tlim[1], tlim[2], length.out=6), 
     labels=round(seq(ylim_right[1], ylim_right[2], length.out=6)/10)*10)



# Overlay polygons for daily forcing quantiles (x = 125:366)
x_chill <- 1:124
x_forcing <- 125:xlim[2]

quantiles <- c(0.1, 0.25, 0.5, 0.75, 0.9)  # Selected quantiles

for (i in 1:(length(quantiles) - 1)) {
  polygon(c(x_chill, rev(x_chill)), 
          tlim[1] + c(daily_chill_plot[, i], rev(daily_chill_plot[, i + 1])) / diff(ylim_right) * diff(tlim),
          col=nom_colors[i], border=NA)
}

lines(x_chill, tlim[1] + (daily_chill_plot[, which(quantiles == 0.5)] / diff(ylim_right)) * diff(tlim), 
      col=c_dark_highlight, lwd=2)  # Highlight the media

for (i in 1:(length(quantiles) - 1)) {
  polygon(c(x_forcing, rev(x_forcing)), 
          tlim[1] + c(daily_forcings_plot[, i], rev(daily_forcings_plot[, i + 1])) / diff(ylim_right) * diff(tlim),
          col=nom_colors[i], border=NA)
}

lines(x_forcing, tlim[1] + (daily_forcings_plot[, which(quantiles == 0.5)] / diff(ylim_right)) * diff(tlim), 
      col=c_dark_highlight, lwd=2)  # Highlight the median

points(xlim[1]:xlim[2], local_temps,
     type="p", col="black", pch=16, cex=0.5, 
     xlim=xlim, xlab="Day", ylim=tlim)

x1 <- mod_data$precursor_events[index]
x2 <- mod_data$events[index]
y <- tlim[1]+1

lines(c(x1, x1), c(y - 0.5, y + 0.5), lwd=2, col=c_mid_teal)
lines(c(x1, x2), c(y, y), lwd=2, col=c_mid_teal)
lines(c(x2, x2), c(y - 0.5, y + 0.5), lwd=2, col=c_mid_teal)

#####################################################################
# plot grid cumulative forcings
#####################################################################
# Plot setup
plot(NA,
     xlim = c(0, 216),
     ylim = c(min(forcings_cumsum_quantiles), max(forcings_cumsum_quantiles)),
     xlab = "",
     ylab = "",
     axes = FALSE)

# X-axis with sensible spacing
axis(1, at = seq(0, 200, by = 20), labels = seq(0, 200, by = 20))
axis(2)

# Plot filled areas for each quantile
for(i in 1:nrow(forcings_cumsum_quantiles)) {
  # Reverse order for polygon to create filled area
  polygon(
    c(1:ncol(forcings_cumsum_quantiles), rev(1:ncol(forcings_cumsum_quantiles))),
    c(as.numeric(forcings_cumsum_quantiles[i,]), rep(min(forcings_cumsum_quantiles), ncol(forcings_cumsum_quantiles))),
    col = color_palette[i],
    border = NA
  )
}

#####################################################################
# plot daily forcings for event prediction
#####################################################################
# Plot setup
plot(NA, 
     xlim = c(0, 366),
     ylim = c(min(daily_forcings_quantiles), 1),
     xlab = "",
     ylab = "",
     axes = FALSE)

# X-axis with sensible spacing
axis(1, at = seq(0, 350, by = 50), labels = seq(0, 350, by = 50))
axis(2)

# Plot filled areas for each quantile
for(i in 1:nrow(daily_forcings_quantiles)) {
  # Reverse order for polygon to create filled area
  polygon(
    c(1:ncol(daily_forcings_quantiles), rev(1:ncol(daily_forcings_quantiles))),
    c(as.numeric(daily_forcings_quantiles[i,]), rep(min(daily_forcings_quantiles), ncol(daily_forcings_quantiles))),
    col = color_palette[i],
    border = NA
  )
}

#####################################################################
# plot daily cumulative forcings for event prediction
#####################################################################
# Plot setup
plot(NA,
     xlim = c(0, 366),
     ylim = c(min(accum_forcings_quantiles), max(accum_forcings_quantiles)),
     xlab = "",
     ylab = "",
     axes = FALSE)

# X-axis with sensible spacing
axis(1, at = seq(0, 350, by = 50), labels = seq(0, 350, by = 50))
axis(2)

# Plot filled areas for each quantile
for(i in 1:nrow(accum_forcings_quantiles)) {
  # Reverse order for polygon to create filled area
  polygon(
    c(1:ncol(accum_forcings_quantiles), rev(1:ncol(accum_forcings_quantiles))),
    c(as.numeric(accum_forcings_quantiles[i,]), rep(min(accum_forcings_quantiles), ncol(accum_forcings_quantiles))),
    col = color_palette[i],
    border = NA
  )
}





# With the ability to suppress the accumulated forcings until the threshold
# the modified survival model exhibits a far better posterior retrodictive
# check.  There appers to be some fine-scale structure in the observed data
# that the modified survival model isn't able to accommodate, but this isn't
# surprising given the crude temperature measurements and no accounting for
# other local ecological conditions across sites and across years.
par(mfrow=c(1, 1))
pred_names <- grep('pred_events', names(mod_surv_samples), value=TRUE)
util$plot_hist_quantiles(samples = mod_surv_samples, val_name_prefix = "pred_events",bin_delta = 5, baseline_values = mod_data$events, xlab = "Phenology Event (new)")

my_quantile(extract_mod$pred_events)
my_quantile(mod_data$events)
mean(extract_mod$pred_events == min(extract_mod$pred_events))

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

# TODO plot actual chilling/forcing vs. pheno
par(mfrow=c(1,1))
plot(mean_chilling, mean_forcing, col=c_dark, pch=16, cex=1.5,
     xlim=c(0.9*min(mean_chilling), 1.1*max(mean_chilling)),
     ylim=c(0.9*min(mean_forcing), c(1.1*max(mean_forcing))))

par(mfrow=c(1,1))
plot(mean_forcing - mean_chilling, mod_data$events, col=c_dark, pch=16, cex=1.5)
x <- mean_forcing - mean_chilling
linear <- lm(mod_data$events ~ x)
xs <- seq(0,6,0.1)
lines(xs, predict(linear, newdata=data.frame(x=xs)))
plot(x, predict(linear, newdata=data.frame(x=x)) - mod_data$events, ylim=c(-6,7), xlim=c(-1,7))
abline(h=0)



# In contrast to the standard survival model, which lead to inferences of
# narrow forcing functions, the modified survival model suggests that wider
# forcing functions are in fact a better fit to the observed data.
par(mfrow=c(2, 2))

util$plot_expectand_pushforward(mod_surv_samples[["Z"]], 25,
                                display_name="Z", flim=c(-5, 15))
xs <- seq(-5, 15, 0.01)
ys <- dnorm(xs,mean_Z , sd_Z)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["Z_chill"]], 25,
                                display_name="Z_chill", flim=c(-5, 15))
xs <- seq(-5, 15, 0.01)
ys <- dnorm(xs,mean_Z , sd_Z)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["Z_forcing"]], 25,
                                display_name="Z_forcing", flim=c(-5, 15))
xs <- seq(-5, 15, 0.01)
ys <- dnorm(xs,mean_Z , sd_Z)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["alpha"]], 25,
                                display_name="alpha", flim=c(0, 100))
xs <- seq(0, 70, 0.01)
ys <- dnorm(xs, mean_alpha, sd_alpha)
lines(xs, ys, lwd=2, col=c_light_teal)
abline(v=alpha_true)

util$plot_expectand_pushforward(mod_surv_samples[["beta"]], 25,
                                display_name="beta", flim=c(-5, 2))
xs <- seq(-5, 2, 0.01)
ys <- dnorm(xs, mean_beta, sd_beta)
lines(xs, ys, lwd=2, col=c_light_teal)
abline(v=beta_true)



util$plot_expectand_pushforward(mod_surv_samples[["sigma"]], 25,
                                display_name="sigma", flim=c(0, 15))
xs <- seq(0, 15, 0.1)
ys <- dnorm(xs, 0, 8 / 2.57)
lines(xs, ys, lwd=2, col=c_light_teal)
abline(v=sigma_true)

# The higher-temperatures in late spring and summer now actively contribute to
# the accumulated forcing of versaison, but the thresholding prevents veraison 
# events from starting until a minimal accumulation is reached.
par(mfrow=c(1, 1))

hist(mod_data$obs_temp, prob=T, breaks=50, col=c_dark_teal, border=c_light_teal,
     main="", xlim=c(-20, 35), xlab="", xaxt='n', 
     ylim=c(0, 0.1), ylab="", yaxt='n')

par(new = T)

f_names <- grep('forcings', names(mod_surv_samples), value=TRUE)
fs <- sapply(f_names, function(f_name) c(t(mod_surv_samples[[f_name]]), recursive=TRUE))

plot_realizations(mod_data$temp_grid, fs, 
                  name="Posterior Forcing Function Realizations",
                  xname="Temperature (C)", 
                  yname="Forcing (FU)", display_ylim=c(0, 1.1))

############################################################
# Retrodictive Comparisons
############################################################

par(mfrow=c(1, 2))
pred_names <- grep('pred_events', names(surv_samples), value=TRUE)
hist_retro(mod_data$events, surv_samples, pred_names, 125, 366, 15, "Phenology Event")

pred_names <- grep('pred_events', names(mod_surv_samples), value=TRUE)
hist_retro(mod_data$events, mod_surv_samples, pred_names, 50, 150, 5, "Phenology Event")

# Save output for paper graphics
d <- 25
B <- 366 / d + 1

idx <- rep(1:B, each=2)
x <- sapply(1:length(idx), function(b) if(b %% 2 == 0) d * idx[b] - 0.5 
            else d * (idx[b] - 1) - 0.5)
breaks <- d * (0:B) - 0.5
bin_min <- min(breaks)
bin_max <- max(breaks)

obs_counts <- hist(mod_data$events, breaks=breaks, plot=FALSE)$counts
pad_obs_counts <- do.call(cbind, lapply(idx, function(n) obs_counts[n]))

cat(sprintf("%.3f/%.3f,", x, pad_obs_counts), "\n")

x1 <- x[2 * (1:25) - 1]
x2 <- x[2 * (1:25)]

pred_names <- grep('pred_events', names(surv_samples), value=TRUE)
pred <- sapply(pred_names,
               function(name) c(t(surv_samples[[name]]), recursive=TRUE))
N <- dim(pred)[1]
pred_counts <- sapply(1:N,
                      function(n) hist(pred[n,][bin_min < pred[n,] & pred[n,] < bin_max],
                                       breaks=breaks,
                                       plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:B,
               function(b) quantile(pred_counts[b,], probs=probs))

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[1,], cred[9,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[2,], cred[8,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[3,], cred[7,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[4,], cred[6,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f,",      x1, x2, cred[5,]), "\n")

pred_names <- grep('pred_events', names(mod_surv_samples), value=TRUE)
pred <- sapply(pred_names,
               function(name) c(t(mod_surv_samples[[name]]), recursive=TRUE))
N <- dim(pred)[1]
pred_counts <- sapply(1:N,
                      function(n) hist(pred[n,][bin_min < pred[n,] & pred[n,] < bin_max],
                                       breaks=breaks,
                                       plot=FALSE)$counts)
probs = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
cred <- sapply(1:B,
               function(b) quantile(pred_counts[b,], probs=probs))

cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[1,], cred[9,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[2,], cred[8,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[3,], cred[7,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f/%.3f,", x1, x2, cred[4,], cred[6,]), "\n")
cat(sprintf("%.3f/%.3f/%.3f,",      x1, x2, cred[5,]), "\n")

############################################################
# Survival Function Comparisons
############################################################

forcing <- function(T_in, T_min, T_opt, T_max, delta) {
  if (T_in < T_min)
    return(0)
  if (T_in > T_max)
    return(0)
  
  if (T_opt > 0.5 * (T_min + T_max)) {
    phi <- (T_max - T_opt) / (T_opt - T_min)
    gamma <- (delta * T_max + T_opt - (1 + delta) * T_min) / (T_max - T_opt)
    a <- (T_in - T_min) / (T_opt - T_min)
    
    if (T_in > T_opt) {
      b <- (T_max - T_in) / (T_max - T_opt);
      c <- b**phi
      return((a * c)**gamma)
    } else {
      b <- (T_max - T_opt) / (T_max - T_in)
      c <- b**(-phi)
      return((a * c)**gamma)
    }
  } else {
    phi <- (T_opt - T_min) / (T_max - T_opt)
    gamma <-  ((1 + delta) * T_max - T_opt - delta * T_min) / (T_opt - T_min)
    b <- (T_max - T_in) / (T_max - T_opt)
    
    if (T_in < T_opt) {
      a <- (T_in - T_min) / (T_opt - T_min)
      c <- a**phi
      return((c * b)**gamma)
    } else {
      a <- (T_opt - T_min) / (T_in - T_min)
      c <- a**(-phi)
      return((c * b)**gamma);
    }
  }
}

logistic <- function(x) {
  if (x > 0) {
    a <- exp(-x)
    1 / (1 + a)
  } else {
    a <- exp(+x)
    a / (1 + a)
  }
}

n <- 82
sample_idxs <- 1:2000
e <- 1

# Accumulated Forcings
nom_samples <- sapply(c("T_min", "T_opt", "T_max", "delta", "gamma"),
                  function(name) c(t(surv_samples[[name]]), recursive=TRUE))

mod_samples <- sapply(c("T_min", "T_opt", "T_max", "delta", "Psi0", "sigma"),
                      function(name) c(t(mod_surv_samples[[name]]), recursive=TRUE))

#nom_post_daily_forcings <- matrix(NA, nrow=4000, ncol=mod_data$N_days[n]) 
mod_post_daily_forcings <- matrix(NA, nrow=4000, ncol=mod_data$N_days[n]) 

for (d in 1:mod_data$N_days[n]) {
  for (r in sample_idxs) {
    # nom_post_daily_forcings[r, d] <- forcing(mod_data$obs_temp[mod_data$temp_start_idxs[n] + d - 1], 
    #                                          nom_samples[r, 1], nom_samples[r, 2], 
    #                                          nom_samples[r, 3], nom_samples[r, 4])
    mod_post_daily_forcings[r, d] <- forcing(mod_data$obs_temp[mod_data$temp_start_idxs[n] + d - 1], 
                                             mod_samples[r, 1], mod_samples[r, 2], 
                                             mod_samples[r, 3], mod_samples[r, 4])
  }
}

nom_post_accum_forcings <- t(sapply(sample_idxs, function(r) cumsum(nom_post_daily_forcings[r,])))
mod_post_accum_forcings <- t(sapply(sample_idxs, function(r) cumsum(mod_post_daily_forcings[r,])))

local_precursor_events <- mod_data$precursor_events[mod_data$event_start_idxs[n]:mod_data$event_end_idxs[n]]
local_events <- mod_data$events[mod_data$event_start_idxs[n]:mod_data$event_end_idxs[n]]

precursor <- local_precursor_events[e]
event <- local_events[e]

nom_post_init_forcing <- nom_post_accum_forcings[,mod_data$precursor_events[e]]
mod_post_init_forcing <- mod_post_accum_forcings[,mod_data$precursor_events[e]]

nom_post_surv <- matrix(NA, nrow=4000, ncol=mod_data$N_days[n]) 
mod_post_surv <- matrix(NA, nrow=2000, ncol=mod_data$N_days[n]) 

for (r in sample_idxs) {
  for (d in 1:mod_data$N_days[n]) {
    if (d < precursor) {
     # nom_post_surv[r, d] = 1
      mod_post_surv[r, d] = 1
    } else {
     # accum <- nom_post_accum_forcings[r, d] - nom_post_init_forcing[r]
     # nom_post_surv[r, d] = exp(- nom_samples[r, 5] * accum)
      
      accum <- mod_post_accum_forcings[r, d] - mod_post_init_forcing[r]
      mod_post_surv[r, d] = 1 - logistic( (accum - mod_samples[r, 5]) / mod_samples[r, 6])
    }
  }
}

par(mfrow=c(1, 1), oma = c(0, 0, 2, 0))

f_names <- grep('forcings', names(surv_samples), value=TRUE)
fs <- sapply(f_names, function(f_name) c(t(surv_samples[[f_name]]), recursive=TRUE))

plot_realizations(1:mod_data$N_days[n], nom_post_surv, 
                  name="Nominal Survival Model",
                  xname="Day of Year", 
                  yname="Survival Function", display_ylim=c(0, 1.1))

abline(v=precursor, lwd=3, col="white")
abline(v=precursor, lwd=2, col=c_dark_teal)

abline(v=event, lwd=3, col="white")
abline(v=event, lwd=2, col=c_dark_teal)

plot_realizations(1:mod_data$N_days[n], mod_post_surv, 
                  name="Modified Survival Model",
                  xname="Day of Year", 
                  yname="Survival Function", display_ylim=c(0, 1.1))

abline(v=precursor, lwd=3, col="white")
abline(v=precursor, lwd=2, col=c_dark_teal)

abline(v=event, lwd=3, col="white")
abline(v=event, lwd=2, col=c_dark_teal)

mtext(paste0("Year ", years[n], ", Event ", e), side=3, line=0.5, outer=TRUE)

# Save output for paper graphics
idx <- 400 * (1:10)
cat(sprintf("%.3f,", nom_post_surv[idx[1],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[2],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[3],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[4],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[5],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[6],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[7],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[8],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[9],]), "\n")
cat(sprintf("%.3f,", nom_post_surv[idx[10],]), "\n")

cat(sprintf("%.3f,", mod_post_surv[idx[1],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[2],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[3],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[4],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[5],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[6],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[7],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[8],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[9],]), "\n")
cat(sprintf("%.3f,", mod_post_surv[idx[10],]), "\n")

############################################################
# Inference Comparisons
############################################################

par(mfrow=c(4, 2))

# T_min
util$plot_expectand_pushforward(surv_samples[["T_min"]], 25,
                                display_name="T_min", flim=c(-5, 15))
xs <- seq(-5, 15, 0.01)
ys <- dnorm(xs, 5, 5 / 2.32)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["T_min"]], 25,
                                display_name="T_min", flim=c(-5, 15))
xs <- seq(-5, 15, 0.01)
ys <- dnorm(xs, 5, 5 / 2.32)
lines(xs, ys, lwd=2, col=c_light_teal)

# T_opt
util$plot_expectand_pushforward(surv_samples[["T_opt"]], 25,
                                display_name="T_opt", flim=c(15, 35))
xs <- seq(15, 35, 0.01)
ys <- dnorm(xs, 25, 5 / 2.32)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["T_opt"]], 25,
                                display_name="T_opt", flim=c(15, 35))
xs <- seq(15, 35, 0.01)
ys <- dnorm(xs, 25, 5 / 2.32)
lines(xs, ys, lwd=2, col=c_light_teal)

# T_max
util$plot_expectand_pushforward(surv_samples[["T_max"]], 25,
                                display_name="T_max", flim=c(20, 45))
xs <- seq(20, 45, 0.01)
ys <- dnorm(xs, 35, 5 / 2.32)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["T_max"]], 25,
                                display_name="T_max", flim=c(20, 45))
xs <- seq(20, 45, 0.01)
ys <- dnorm(xs, 35, 5 / 2.32)
lines(xs, ys, lwd=2, col=c_light_teal)

# delta
util$plot_expectand_pushforward(surv_samples[["delta"]], 25,
                                display_name="delta", flim=c(0, 15))
xs <- seq(0, 15, 0.01)
ys <- dlnorm(xs, 0.5, 0.5)
lines(xs, ys, lwd=2, col=c_light_teal)

util$plot_expectand_pushforward(mod_surv_samples[["delta"]], 25,
                                display_name="delta", flim=c(0, 15))
xs <- seq(0, 15, 0.01)
ys <- dlnorm(xs, 0.5, 0.5)
lines(xs, ys, lwd=2, col=c_light_teal)

############################################################
# Inferred Forcing Comparisons
############################################################



mod_fs <- sapply(f_names, function(f_name) c(t(mod_surv_samples[[f_name]]), recursive=TRUE))
plot_realizations(mod_data$temp_grid, mod_fs, 
                  name="Posterior\nForcing Function Realizations",
                  xname="Temperature (C)", 
                  yname="Forcing (FU)", display_ylim=c(0, 1.1))



