

source(file.path('~/Documents/git/mcmc_visualization_tools/r/mcmc_analysis_tools_rstan.R'), local=util)
source(file.path('~/Documents/git/mcmc_visualization_tools/r/mcmc_visualization_tools.R'), local=util)

library(rstan)
rstan_options(auto_write = TRUE)            # Cache compiled Stan programs
options(mc.cores = parallel::detectCores()) # Parallelize chains
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")

library(openmeteo)
library(stringr)
library(colormap)
library(scales)


par(family="serif", las=1, bty="l", 
    cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

prep_base_data <- function() {
  d_wash <- read.csv("data/washingtondc.csv")
  d_wash$dataset <- 'washingtondc'
  d_wash$bloom_percentage <- 0.70
  
  d_liestal <- read.csv("data/liestal.csv")
  d_liestal$dataset <- 'liestal'
  d_liestal$bloom_percentage <- 0.25
  
  d_kyoto <- read.csv("data/kyoto.csv")
  d_kyoto$dataset <- 'kyoto'
  d_kyoto$bloom_percentage <- 0.80
  
  d_vancouver <- read.csv("data/vancouver.csv")
  d_vancouver$dataset <- 'vancouver'
  d_vancouver$bloom_percentage <- 0.70
  
  d_nyc <- read.csv("data/nyc.csv")
  d_nyc$dataset <- 'nyc'
  d_nyc$bloom_percentage <- 0.70
  
  d_bloom <- rbind(
    d_wash,
    d_liestal,
    d_kyoto,
    d_vancouver,
    d_nyc
  )
  
  # get min and max years for each location
  locations_all <- data.table::as.data.table(d_bloom)[,j=list(
    n=.N,
    min_year=min(year),
    max_year=max(year)), 
    by=c('location', 'lat', 'long')]

  locations_all <- as.data.frame(locations_all)
  
  list(locations_init=locations_all, d_init=d_bloom)
  
}


pull_daily_temperature_data <- function(locations_df, out_file, max_date, min_date=NA, write_to_file=F) {
  if (write_to_file == F) {
    return('already pulled daily temperatures...')
  }
  # this loop takes several hours over multiple days to complete given the API
  # limits
  temp_daily_list <- list()
  k <- 1
  total_iter <- 0
  max_iter <- 1000
  while(k <= nrow(locations_df)) {
    row <- locations_df[k,]
    loc <- as.numeric(row[,c('lat','long')])
    if (is.na(min_date)) {
      date_min <- max(paste0(str_pad(row[,'min_year']-1, width=4, 'left', pad='0'), '-01-01'), '1940-01-01')
    } else {
      date_min <- min_date
    }
    if (is.na(max_date)) {
      date_max <- max(c(paste0(row[,'max_year'], '-12-31'), '1940-01-01'))
    } else {
      date_max <- max_date
    }
    tryCatch( {
      temp <- get_temperature(loc, date_min, date_max)
      temp$location <- row[,'location']
      temp$timestamp_run <- Sys.time()

      
      if(!file.exists(out_file)) {
        print('creating temperature file...')
        write.csv(temp, file = out_file, row.names = F)
      } else {
        write.table(temp, sep=',', file = out_file, append=T, row.names = F, col.names = F)
      }
      k <- k+1
    },
    error=function(cond) {
      m <- cond$message
      print(m)
      if(grepl('minutely', m, ignore.case = T)) {
        print(paste0('waiting 60 seconds on location ', k, '...'))
        Sys.sleep(60)
      }
      if(grepl('hourly', m, ignore.case = T)) {
        print(paste0('waiting an hour on location ', k, '...'))
        Sys.sleep(60*60)
      }
      if (grepl('daily', m, ignore.case = T)) {
        print(paste0('waiting a day on location ', k, '...'))
        Sys.sleep(60*60*24)          
      }
    }
    )
    total_iter <- total_iter+1
    temp <- NULL
    print(paste0('total iterations = ', total_iter, '...'))
    if (total_iter > max_iter) {
      stop('max iterations reached')
    }
  }
  
}



pull_forecasted_temperature_data <- function(locations_df, out_file, min_date='2025-02-04', max_date='2025-06-01', write_to_file=F) {
  if (write_to_file == F) {
    return('already pulled daily forecasted temperatures...')
  }
  # this loop takes several hours over multiple days to complete given the API
  # limits
  temp_daily_list <- list()
  k <- 1
  total_iter <- 0
  max_iter <- 1000
  while(k <= nrow(locations_df)) {
    row <- locations_df[k,]
    loc <- as.numeric(row[,c('lat','long')])
    
    if (is.na(min_date)) {
      date_min <- max(paste0(str_pad(row[,'min_year']-1, width=4, 'left', pad='0'), '-01-01'), '1950-01-01')
    } else {
      date_min <- min_date
    }
    if (is.na(max_date)) {
      date_max <- max(c(paste0(row[,'max_year'], '-12-31'), '1950-01-01'))
    } else {
      date_max <- max_date
    }

    tryCatch( {
      temp <- get_forecasted_temperature(loc, date_min, date_max)
      temp$location <- row[,'location']
      temp$timestamp_run <- Sys.time()
      
      if(!file.exists(out_file)) {
        print('creating temperature file...')
        write.csv(temp, file = out_file, row.names = F)
      } else {
        write.table(temp, sep=',', file = out_file, append=T, row.names = F, col.names = F)
      }
      k <- k+1
    },
    error=function(cond) {
      m <- cond$message
      print(m)
      if(grepl('minutely', m, ignore.case = T)) {
        print(paste0('waiting 60 seconds on location ', k, '...'))
        Sys.sleep(60)
      }
      if(grepl('hourly', m, ignore.case = T)) {
        print(paste0('waiting an hour on location ', k, '...'))
        Sys.sleep(60*60)
      }
      if (grepl('daily', m, ignore.case = T)) {
        print(paste0('waiting a day on location ', k, '...'))
        Sys.sleep(60*60*24)          
      }
    }
    )
    total_iter <- total_iter+1
    temp <- NULL
    print(paste0('total iterations = ', total_iter, '...'))
    if (total_iter > max_iter) {
      stop('max iterations reached')
    }
  }
  
}


get_temperature <- function (loc, date_min, date_max) {
  temp <- as.data.frame(openmeteo::weather_history(location=loc, start=date_min, end=date_max, daily=c(
    "temperature_2m_mean"
  )))
  
  temp
}

get_forecasted_temperature <- function(loc, date_min, date_max, models=NA) {
  if (is.na(models)) {
    models <- c(
      'CMCC_CM2_VHR4',
      'FGOALS_f3_H',
      'HiRAM_SIT_HR',
      'MRI_AGCM3_2_S',
      'EC_Earth3P_HR',
      'MPI_ESM1_2_XR',
      'NICAM16_8S'
    )
  }
  temp_out <- data.frame(matrix(ncol=3, nrow=0))
  names(temp_out) <- c("date", "daily_temperature_2m_mean", "model")
  for (j in models) {
    temp <- as.data.frame(openmeteo::climate_forecast(location=loc, start=date_min, end=date_max, daily=c(
      "temperature_2m_mean"
    ), model = j))
    temp$model <- j
    
    temp_out <- rbind(temp_out, temp)
  }

  
  temp_out
}






my_quantile <- function(x, na.rm=F) {
  if(any(is.na(x)))
    print('removing NAs')
  quantile(x, c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99), na.rm=na.rm)
}

# Convert celsius to fahrenheit
c_to_f <- function(celsius) {
  # Input validation
  if (!is.numeric(celsius)) {
    stop("Input must be numeric")
  }
  
  # Conversion formula: °F = (°C × 9/5) + 32
  fahrenheit <- (celsius * 9/5) + 32
  
  return(fahrenheit)
}

# March 19 - June 20
day_num_of_year <- function(date, base_date = NA) {
  if(is.na(base_date)) {
    base_date <- as.Date(paste0(substr(date, 1, 4), '-01-01'))
  }
  as.integer(as.Date(date) - as.Date(base_date)) + 1
}

# logistic function
logit <- function(x, x0=0, k=1) {
  1 / (1 + exp(-k*(x-x0)))
}


# Concatenate ragged site-year data into a single array 
format_site_data <- function(site_temp_data, site_pheno_data, site_years, site_years_predict, site_temp_data_predict, site_forecast_temp_data) {
  N_site_years <- nrow(site_years)
  
  max_chill_days <- max(site_pheno_data$startdate)
  
  # Concatenate ragged temperature data into a monolithic array.
  N_days <- sapply(1:nrow(site_years), 
                   function(d) length(site_temp_data$day[site_temp_data$year == site_years[d,2] & site_temp_data$location == site_years[d,1]]))
  
  temp_start_idxs <- rep(1, N_site_years)
  temp_end_idxs <- rep(1, N_site_years)
  
  N_obs_temp <- sum(N_days)
  obs_temp <- rep(NA, N_obs_temp)
  
  start_idx <- 1
  end_idx <- start_idx + N_days[1] - 1
  
  temp_start_idxs[1] <- start_idx
  temp_end_idxs[1] <- end_idx
  obs_temp[start_idx:end_idx] <- 
    site_temp_data$value[site_temp_data$year == site_years[1,2] & site_temp_data$location == site_years[1,1]]
  
  for (n in 2:N_site_years) {
    start_idx <- temp_end_idxs[n - 1] + 1
    end_idx <- start_idx + N_days[n] - 1
    
    temp_start_idxs[n] <- start_idx
    temp_end_idxs[n] <- end_idx
    
    obs_temp[start_idx:end_idx] <- 
      site_temp_data$value[site_temp_data$year == site_years[n,2] & site_temp_data$location == site_years[n,1]]
  }
  
  # observed temperatures for prediction (no historical event)
  N_site_years_predict <- nrow(site_years_predict)
  
  N_days_predict <- sapply(1:nrow(site_years_predict), 
                   function(d) length(site_temp_data_predict$day[site_temp_data_predict$year == site_years_predict[d,2] & site_temp_data_predict$location == site_years_predict[d,1]]))
  
  temp_start_idxs_predict <- rep(1, N_site_years_predict)
  temp_end_idxs_predict <- rep(1, N_site_years_predict)
  
  N_obs_temp_predict <- sum(N_days_predict)
  obs_temp_predict <- rep(NA, N_obs_temp_predict)
  
  start_idx_predict <- 1
  end_idx_predict <- start_idx_predict + N_days_predict[1] - 1
  
  temp_start_idxs_predict[1] <- start_idx_predict
  temp_end_idxs_predict[1] <- end_idx_predict
  obs_temp_predict[start_idx_predict:end_idx_predict] <- 
    site_temp_data_predict$value[site_temp_data_predict$year == site_years_predict[1,2] & site_temp_data_predict$location == site_years_predict[1,1]]
  
  for (n in 2:N_site_years_predict) {
    start_idx_predict <- temp_end_idxs_predict[n - 1] + 1
    end_idx_predict <- start_idx_predict + N_days_predict[n] - 1
    
    temp_start_idxs_predict[n] <- start_idx_predict
    temp_end_idxs_predict[n] <- end_idx_predict
    
    obs_temp_predict[start_idx_predict:end_idx_predict] <- 
      site_temp_data_predict$value[site_temp_data_predict$year == site_years_predict[n,2] & site_temp_data_predict$location == site_years_predict[n,1]]
  }
  
  # observed temperatures for prediction (no historical event)
  N_site_years_predict <- nrow(site_years_predict)
  
  N_days_predict <- sapply(1:nrow(site_years_predict), 
                   function(d) length(site_temp_data_predict$day[site_temp_data_predict$year == site_years_predict[d,2] & site_temp_data_predict$location == site_years_predict[d,1]]))
  
  temp_start_idxs_predict <- rep(1, N_site_years_predict)
  temp_end_idxs_predict <- rep(1, N_site_years_predict)
  
  N_obs_temp_predict <- sum(N_days_predict)
  obs_temp_predict <- rep(NA, N_obs_temp_predict)
  
  start_idx_predict <- 1
  end_idx_predict <- start_idx_predict + N_days_predict[1] - 1
  
  temp_start_idxs_predict[1] <- start_idx_predict
  temp_end_idxs_predict[1] <- end_idx_predict
  obs_temp_predict[start_idx_predict:end_idx_predict] <- 
    site_temp_data_predict$value[site_temp_data_predict$year == site_years_predict[1,2] & site_temp_data_predict$location == site_years_predict[1,1]]
  
  for (n in 2:N_site_years_predict) {
    start_idx_predict <- temp_end_idxs_predict[n - 1] + 1
    end_idx_predict <- start_idx_predict + N_days_predict[n] - 1
    
    temp_start_idxs_predict[n] <- start_idx_predict
    temp_end_idxs_predict[n] <- end_idx_predict
    
    obs_temp_predict[start_idx_predict:end_idx_predict] <- 
      site_temp_data_predict$value[site_temp_data_predict$year == site_years_predict[n,2] & site_temp_data_predict$location == site_years_predict[n,1]]
  }
  
  # Concatenate ragged event data into a monolithic array.
  N_events <- rep(NA, N_site_years)
  event_start_idxs <- rep(NA, N_site_years)
  event_end_idxs <- rep(NA, N_site_years)
  
  N <- nrow(site_pheno_data)
  precursor_events <- rep(NA, N)
  events <- rep(NA, N)
  
  N_events[1] <- nrow(site_pheno_data[site_pheno_data$year == site_years[1,2] & site_pheno_data$location == site_years[1,1],])
  start_idx <- 1
  end_idx <- start_idx + N_events[1] - 1
  
  event_start_idxs[1] <- start_idx
  event_end_idxs[1] <- end_idx
  precursor_events[start_idx:end_idx] <- 
    site_pheno_data$startdate[site_pheno_data$year == site_years[1,2] & site_pheno_data$location == site_years[1,1]]
  events[start_idx:end_idx] <- 
    site_pheno_data$doy[site_pheno_data$year == site_years[1,2] & site_pheno_data$location == site_years[1,1]]
  
  for (n in 2:N_site_years) {
    N_events[n] <- nrow(site_pheno_data[site_pheno_data$year == site_years[n,2] & site_pheno_data$location == site_years[n,1],])
    start_idx <- event_end_idxs[n - 1] + 1
    end_idx <- start_idx + N_events[n] - 1
    
    event_start_idxs[n] <- start_idx
    event_end_idxs[n] <- end_idx
    
    precursor_events[start_idx:end_idx] <- 
      site_pheno_data$startdate[site_pheno_data$year == site_years[n,2] & site_pheno_data$location == site_years[n,1]]
    events[start_idx:end_idx] <- 
      site_pheno_data$doy[site_pheno_data$year == site_years[n,2] & site_pheno_data$location == site_years[n,1]]
  }
  
  N_events_predict <- rep(1, length.out=length(unique(site_temp_data$location)))
  
  
  # find common obs and forecasted temps
  forecast_prediction_start_date <- as.Date(max(site_temp_data_predict$date))+1
  obs_temp_w_forecast_df <- merge(site_temp_data[,c('location', 'date', 'value')], site_forecast_temp_data, by=c('location', 'date'),all=F, suffixes = c("_obs", "_forecast"))
  pred_temp_w_forecast_df <- merge(site_temp_data_predict[,c('location', 'date', 'value')], site_forecast_temp_data, by=c('location', 'date'),all=F, suffixes = c("_obs", "_forecast"))
  forecast_w_no_obs_temp_df <- merge(site_temp_data_predict[,c('location', 'date', 'value')], site_forecast_temp_data, by=c('location', 'date'),all.x=F, all.y=T, suffixes = c("_obs", "_forecast"))
  forecast_w_no_obs_temp_df <- forecast_w_no_obs_temp_df[is.na(forecast_w_no_obs_temp_df$value_obs) & forecast_w_no_obs_temp_df$date >= forecast_prediction_start_date,]
  
  
  
  # historical obs temps and forecasts with events
  forecast_mean_w_obs_temp <- obs_temp_w_forecast_df$value_forecast[,'mean']
  forecast_sd_w_obs_temp <- obs_temp_w_forecast_df$value_forecast[,'sd']
  obs_temp_w_forecast <- obs_temp_w_forecast_df$value_obs
  
  # historical obs temps and forecasts NO events
  forecast_mean_w_pred_temp <- pred_temp_w_forecast_df$value_forecast[,'mean']
  forecast_sd_w_pred_temp <- pred_temp_w_forecast_df$value_forecast[,'sd']
  pred_temp_w_forecast <- pred_temp_w_forecast_df$value_obs
  
  # forecasts NO events or observed temps
  forecast_mean_w_no_obs_temp <- forecast_w_no_obs_temp_df$value_forecast[,'mean']
  forecast_sd_w_no_obs_temp <- forecast_w_no_obs_temp_df$value_forecast[,'sd']

  
  
  # get forecast indices
  # forecasts for prediction (NO historical event + NO obs temp)
  N_site_years_predict_forecast <- nrow(site_years_predict)
  
  N_days_predict_forecast <- sapply(1:nrow(site_years_predict), 
                           function(d) length(forecast_w_no_obs_temp_df$day[forecast_w_no_obs_temp_df$year == site_years_predict[d,2] & forecast_w_no_obs_temp_df$location == site_years_predict[d,1]]))
  
  temp_start_idxs_predict_forecast <- rep(1, N_site_years_predict_forecast)
  temp_end_idxs_predict_forecast <- rep(1, N_site_years_predict_forecast)
  
  N_forecast_temp_predict <- sum(N_days_predict_forecast)
  forecast_temp_predict <- rep(NA, N_forecast_temp_predict)
  
  start_idx_predict_forecast <- 1
  end_idx_predict_forecast <- start_idx_predict_forecast + N_days_predict_forecast[1] - 1
  
  temp_start_idxs_predict_forecast[1] <- start_idx_predict_forecast
  temp_end_idxs_predict_forecast[1] <- end_idx_predict_forecast

  for (n in 2:N_site_years_predict_forecast) {
    start_idx_predict_forecast <- temp_end_idxs_predict_forecast[n - 1] + 1
    end_idx_predict_forecast <- start_idx_predict_forecast + N_days_predict_forecast[n] - 1
    
    temp_start_idxs_predict_forecast[n] <- start_idx_predict_forecast
    temp_end_idxs_predict_forecast[n] <- end_idx_predict_forecast
   }

   list("N_days" = N_days, 
       "temp_start_idxs" = temp_start_idxs, "temp_end_idxs" = temp_end_idxs,
       "obs_temp" = obs_temp,
       "N_events" = N_events,
       "event_start_idxs" = event_start_idxs, "event_end_idxs" = event_end_idxs,
       "precursor_events" = precursor_events, "events" = events,
       "max_chill_days" = max_chill_days,
       "N_days_predict" = N_days_predict,
       "temp_start_idxs_predict" = temp_start_idxs_predict, "temp_end_idxs_predict" = temp_end_idxs_predict,
       "obs_temp_predict" = obs_temp_predict,
       "N_events_predict" = N_events_predict,
       "N_site_years_predict" = N_site_years_predict,
       "N_site_years_predict_forecast" = N_site_years_predict_forecast,

       # historical obs temps and forecasts with events
       "obs_temp_w_forecast" = obs_temp_w_forecast, # historical observed temps with forecast
       "pred_temp_w_forecast" = pred_temp_w_forecast, # historical observed temp for prediction (no event) with forecast
       "forecast_mean_w_obs_temp" = forecast_mean_w_obs_temp, 
       "forecast_sd_w_obs_temp" = forecast_sd_w_obs_temp,
       "obs_temp_w_forecast" = obs_temp_w_forecast,
       
       # historical obs temps and forecasts NO events
       "forecast_mean_w_pred_temp" = forecast_mean_w_pred_temp,
       "forecast_sd_w_pred_temp" = forecast_sd_w_pred_temp,
       "pred_temp_w_forecast" = pred_temp_w_forecast ,
       
       # forecasts NO events or observed temps
       "forecast_mean_w_no_obs_temp" = forecast_mean_w_no_obs_temp,
       "forecast_sd_w_no_obs_temp" = forecast_sd_w_no_obs_temp,
       "N_forecast_temp_predict" =   N_forecast_temp_predict,
       "temp_start_idxs_predict_forecast" = temp_start_idxs_predict_forecast,
       "temp_end_idxs_predict_forecast" = temp_end_idxs_predict_forecast,
       "N_days_predict_forecast" = N_days_predict_forecast,
       
       "obs_temp_w_forecast_df" = obs_temp_w_forecast_df,
       "pred_temp_w_forecast_df" = pred_temp_w_forecast_df,
       "forecast_w_no_obs_temp_df" = forecast_w_no_obs_temp_df
  )
}


# Plot posterior samples of functional behavior.
plot_realizations <- function(xs, fs, name="", 
                              xname=NULL, display_xlim=NULL, 
                              yname=NULL, display_ylim=NULL) {
  if (is.null(display_xlim)) {
    display_xlim <- c(xs[1], xs[length(xs)])
  }
  
  if (is.null(display_ylim)) {
    display_ylim <- range(c(fs, recursive=TRUE))
  }
  
  I <- length(fs[,1])
  J <- min(50, I)
  plot_idx <- sapply(1:J, function(j) (I %/% J) * (j - 1) + 1)
  
  line_colors <- colormap(colormap=nom_colors, nshades=J)
  
  plot(1, type="n", xlab=xname, ylab=yname, main=name,
       xlim=display_xlim, ylim=display_ylim)
  for (j in 1:J)
    lines(xs, fs[plot_idx[j],], col=line_colors[j], lwd=2)
}

prep_mod_data <- function(pheno_data_init, temp_data_init, forecast_temp_data_init, plot_temp=F, plot_forecast_temp=F) {
  
  years <- 1940:2024
  years <- years[(years - 1) %in% years]
  
  # Read in phenology data
  bloom_cycle_start_month_day <- "10-01"
  pheno_data <- pheno_data_init
  pheno_data$year = as.numeric(substr(pheno_data$bloom_date, 1, 4))
  pheno_data$is_leap <- (pheno_data$year %% 4 == 0 & (pheno_data$year %% 100 != 0 | pheno_data$year %% 400 == 0))
  
  pheno_data$bloom_cycle_start_date <- as.Date(paste0(pheno_data$year - 1, '-', bloom_cycle_start_month_day))
  pheno_data$bloom_cycle_end_date <- as.Date(pheno_data$bloom_cycle_start_date + ifelse(pheno_data$is_leap, 366, 365) - 1)
  pheno_data$doy_init <- pheno_data$bloom_doy
  pheno_data$doy <- as.integer(as.Date(pheno_data$bloom_date) - pheno_data$bloom_cycle_start_date + 1)
  
  # for priors
  pheno_data_no_temp <- pheno_data[pheno_data$year < min(years),]
  
  pheno_data <- pheno_data[pheno_data$year %in% years,]
  pheno_data$startdate <- 124 # as.integer(as.Date(paste0(pheno_data$year, '-02-01')) -  pheno_data$bloom_cycle_start_date + 1)

  # temperature data.
  temp_data <- temp_data_init[temp_data_init$location %in% unique(pheno_data$location),]
  temp_data$year_init <- as.numeric(substr(temp_data$date, 1, 4))
  temp_data$month <- as.numeric(substr(temp_data$date, 6, 7))
  temp_data$dom <- as.numeric(substr(temp_data$date, 9, 10))
  temp_data$value <- temp_data$daily_temperature_2m_mean
  temp_data$day <- as.numeric(format(as.Date(temp_data$date), "%j"))
  
  unique_temp_dates <- sort(unique(temp_data$date))
  map_temp_dates_to_year <- data.frame(date=as.Date(unique_temp_dates), year=sapply(unique_temp_dates, function(d) {unique(
    pheno_data$year[pheno_data$bloom_cycle_start_date <= d & d <= pheno_data$bloom_cycle_end_date])[1]}))
  
  for (k in which(is.na(map_temp_dates_to_year$year))) {
    year_k <- as.integer(format(as.Date(map_temp_dates_to_year$date[k]), "%Y"))
    if (year_k >= 2024)
      { map_temp_dates_to_year$year[k] <- 2025 }
  }
  map_temp_dates_to_year <- na.omit(map_temp_dates_to_year)
  
  # prediction year temps
  temp_data_predict <- merge(temp_data, map_temp_dates_to_year[map_temp_dates_to_year$year == 2025,], by='date', all.x=F, all.y=F)
  temp_data_predict <- temp_data_predict[temp_data_predict$year == 2025,] 
  temp_data_predict <- temp_data_predict[,c('year_init', 'year', 'month', 'dom', 'value', 'day', 'location', 'date')]
  
  # historical bloom year temps
  temp_data <- merge(temp_data, map_temp_dates_to_year[map_temp_dates_to_year$year < 2025,], by='date', all.x=T, all.y=F)
  temp_data <- temp_data[temp_data$year %in% years,]
  temp_data <- temp_data[,c('year_init', 'year', 'month', 'dom', 'value', 'day', 'location', 'date')]
  
  # forecasted temp data
  forecast_temp_data <- forecast_temp_data_init[forecast_temp_data_init$location %in% unique(pheno_data$location),]
  forecast_temp_data$year_init <- as.numeric(substr(forecast_temp_data$date, 1, 4))
  forecast_temp_data$month <- as.numeric(substr(forecast_temp_data$date, 6, 7))
  forecast_temp_data$dom <- as.numeric(substr(forecast_temp_data$date, 9, 10))
  forecast_temp_data$value <- forecast_temp_data$daily_temperature_2m_mean
  forecast_temp_data$day <- as.numeric(format(as.Date(forecast_temp_data$date), "%j"))
  
  unique_forecast_temp_dates <- sort(unique(forecast_temp_data$date))
  map_forecast_temp_dates_to_year <- data.frame(date=as.Date(unique_forecast_temp_dates), year=sapply(unique_forecast_temp_dates, function(d) {unique(
    pheno_data$year[pheno_data$bloom_cycle_start_date <= d & d <= pheno_data$bloom_cycle_end_date])[1]}))
  
  for (k in which(is.na(map_forecast_temp_dates_to_year$year))) {
    year_k <- as.integer(format(as.Date(map_forecast_temp_dates_to_year$date[k]), "%Y"))
    if (year_k >= 2024)
    { map_forecast_temp_dates_to_year$year[k] <- 2025 }
  }
  
  
  forecast_temp_data <- merge(forecast_temp_data, map_forecast_temp_dates_to_year, by='date', all.x=T, all.y=F)
  forecast_temp_data <- forecast_temp_data[,c('year_init', 'year', 'month', 'dom', 'value', 'day', 'location', 'date')]  
  
  # Order temp and pheno data
  temp_data <- temp_data[order(temp_data$location, temp_data$date),]
  temp_data_predict <- temp_data_predict[order(temp_data_predict$location, temp_data_predict$date),]
  forecast_temp_data <- forecast_temp_data[order(forecast_temp_data$location, forecast_temp_data$date),]
  pheno_data <- pheno_data[order(pheno_data$location, pheno_data$year),]

  # Set valid years.
  sites_and_years <- unique(pheno_data[,c('location', 'year')])
  sites_and_years <- sites_and_years[order(sites_and_years$location, sites_and_years$year),]
  
  sites_and_years_predict <- unique(temp_data_predict[,c('location', 'year')])
  
  # Concatenate yearly data.
  site_data <- format_site_data(site_temp_data = temp_data,
                                site_pheno_data = pheno_data,
                                site_years =  sites_and_years,
                                site_years_predict = sites_and_years_predict,
                                site_temp_data_predict = temp_data_predict,
                                site_forecast_temp_data = forecast_temp_data)
  site_data$N_site_years <- nrow(sites_and_years)
  
  if (plot_temp) {
    # Visualize data.
    tlim <- range(site_data$obs_temp)
    # xlim <- c(min(site_data$precursor_events) - 5, max(site_data$events) + 5)
    xlim <- c(1, max(site_data$events)+1)
    
    par(mfrow=c(2, 2))
    
    for (n in 1:site_data$N_site_years) {
      local_temps <- site_data$obs_temp[site_data$temp_start_idxs[n]:
                                          site_data$temp_end_idxs[n]]
      local_precursor_events <- site_data$precursor_events[site_data$event_start_idxs[n]:
                                                             site_data$event_end_idxs[n]]
      local_events <- site_data$events[site_data$event_start_idxs[n]:
                                         site_data$event_end_idxs[n]]
      
      title <- paste0("Site: ", sites_and_years$location[n], ", Year: ", sites_and_years$year[n])
      
      plot(1:site_data$N_days[n], local_temps, main=title,
           type="p", col="black", pch=16, cex=0.5, 
           xlim=xlim, xlab="Day", ylim=tlim, ylab="Recorded Temperature (C)")
      
      for (e in 1:site_data$N_events[n]) {
        x1 <- local_precursor_events[e]
        x2 <- local_events[e]
        
        y <- tlim[1] + (2 * e + -1)
        
        lines(c(x1, x1), c(y - 0.5, y + 0.5), lwd=2, col=c_mid_teal)
        lines(c(x1, x2), c(y, y), lwd=2, col=c_mid_teal)
        lines(c(x2, x2), c(y - 0.5, y + 0.5), lwd=2, col=c_mid_teal)
      }
    }
  }
  
  if (plot_forecast_temp) {
    # Visualize data.
    tlim <- range(
                  site_data$forecast_mean_w_obs_temp + 2*site_data$forecast_sd_w_obs_temp,
                  site_data$forecast_mean_w_obs_temp - 2*site_data$forecast_sd_w_obs_temp,
                  site_data$forecast_mean_w_no_obs_temp + 2*site_data$forecast_sd_w_no_obs_temp,
                  site_data$forecast_mean_w_no_obs_temp - 2*site_data$forecast_sd_w_no_obs_temp
                  )
    # xlim <- c(min(site_data$precursor_events) - 5, max(site_data$events) + 5)
    xlim <- c(1, 366)
    
    par(mfrow=c(2, 2))
    
    forecast_comparison_df <- rbind(site_data$obs_temp_w_forecast_df, site_data$pred_temp_w_forecast_df)
    
    sites_years_forecast_comparison <- unique(forecast_comparison_df[,c('location', 'year')])
    for (n in 1:nrow(sites_years_forecast_comparison)) {
      local_temps <- forecast_comparison_df[forecast_comparison_df$location == sites_years_forecast_comparison$location[n] & 
                                              forecast_comparison_df$year == sites_years_forecast_comparison$year[n],]
      
      title <- paste0("Site: ",  sites_years_forecast_comparison$location[n], ", Year: ", sites_years_forecast_comparison$year[n])
      
      plot(1:nrow(local_temps), main=title,
           col="black", pch=16, cex=0.75, 
           xlim=xlim, xlab="Day", ylim=tlim, ylab="Recorded Temperature (C)", type='n')
      
      # Draw the polygon (shaded confidence region)
      polygon(c(1:nrow(local_temps), rev(1:nrow(local_temps))),
              c(local_temps$value_forecast[,'mean'] + 2 * local_temps$value_forecast[,'sd'],
                rev(local_temps$value_forecast[,'mean'] - 2 * local_temps$value_forecast[,'sd'])),
              col = c_light, border = NA)
      
      lines(1:nrow(local_temps), local_temps$value_forecast[,'mean'], main=title,
            col=c_mid, pch=16, lwd=3, 
            xlim=xlim, xlab="Day", ylim=tlim)
      
      
      points(1:nrow(local_temps), local_temps$value_obs, main=title,
           type="p", col="black", pch=16, cex=0.75, 
           xlim=xlim, xlab="Day", ylim=tlim, ylab="Recorded Temperature (C)")
 
      
    }
  }
 
  
  data <- list(
               "N_site_years" = site_data$N_site_years, 
               "N_obs_temps" = length(site_data$obs_temp),
               "obs_temps" = site_data$obs_temp,
               "temp_start_idxs" = site_data$temp_start_idxs,
               "temp_end_idxs" = site_data$temp_end_idxs,
               "N_days" = site_data$N_days,
               "N" = length(site_data$events),
               "N_events" = site_data$N_events,
               "precursor_events" = site_data$precursor_events, 
               "events" = site_data$events,
               "event_start_idxs" = site_data$event_start_idxs,
               "event_end_idxs" = site_data$event_end_idxs,
               "N_site_years_predict" = site_data$N_site_years_predict,
               "N_obs_temps_predict" = length(site_data$obs_temp_predict),
               "obs_temps_predict" = site_data$obs_temp_predict,
               "temp_start_idxs_predict" = site_data$temp_start_idxs_predict,
               "temp_end_idxs_predict" = site_data$temp_end_idxs_predict,
               "N_days_predict" = site_data$N_days_predict,
               "N_predict" = length(site_data$N_days_predict),
               "D" = site_data$max_chill_days,
               "temp_data_df" = temp_data,
               "pheno_data_df" = pheno_data,
               "sites_and_years_df" = sites_and_years,
               "pheno_data_prior_df" = pheno_data_no_temp,
               "temp_data_predict_df" = temp_data_predict,
               "forecast_temp_data_df" = forecast_temp_data
  )
  
  return(data)
  
  
}









