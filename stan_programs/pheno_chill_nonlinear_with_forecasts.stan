functions {
  // Forcing function - logistic
  real forcing(real T, real Z, real k) {
    if (T - Z >= 0) {
      return 1 / (1 + exp(k*(Z - T)));
    } else {
      return exp(k*(T - Z)) / (1 + exp(k*(T - Z)));
    }
  }
  
  // Log of forcing function
  real log_forcing(real T, real Z, real k) {
    return -log1p_exp(k*(Z - T));
  }
  
  
  // Chilling function
  real chill_days(real T, real Z) {
    if (Z - T >= 0) {
      return 1 / (1 + exp(T - Z));
    } else {
      return exp(Z - T) / (1 + exp(Z - T));
    }
  }
  
}

data {

  int<lower=1> N_site_years; // Number of year-site instances
  
  int<lower=1> N_obs_temps;   // Number of year-site days
  real obs_temps[N_obs_temps]; // Recorded temperature (C) for each day
  int<lower=1, upper=N_obs_temps> temp_start_idxs[N_site_years];
  int<lower=1> N_days[N_site_years];
  
  // Number of phenological events across all sites and years
  int<lower=1> N;
  
  // Day when phenological transition starts
  // Max equals 366 to account for leap years
  int<lower=1, upper=366> precursor_events[N];
  
  // Day when phenological transition ends
  // Max equals 366 to account for leap years
  int<lower=1, upper=366> events[N];
  
  int<lower=1, upper=N> event_start_idxs[N_site_years];
  int<lower=1> N_events[N_site_years];
  
  // prior hyperparameters
  real mean_Z_chill;
  real<lower=0> sd_Z_chill;

  real mean_Z_forcing;
  real<lower=0> sd_Z_forcing;
  
  real<lower=0> mean_alpha;
  real<lower=0> sd_alpha;
  
  real mean_beta1;
  real<lower=0> sd_beta1;
  
  real mean_beta2;
  real<lower=0> sd_beta2;
  
  real<lower=0> sd_sigma;
  
  real<lower=0> x0_chill_days;


  // temp forecasting
  real mean_delta;
  real sd_delta;
  real<lower=0> sd_tau;
  
  // event prediction
  int<lower=1> N_site_years_predict; // Number of year-site instances to predict
  
  int precursor_events_predict[N_site_years_predict];
  
  // predict temps (observed)
  int<lower=1> N_obs_temps_predict; // observed temps
  real obs_temps_predict[N_obs_temps_predict]; // Recorded temperature (C) for each day
  int<lower=1, upper=N_obs_temps_predict> temp_start_idxs_predict[N_site_years_predict];
  int<lower=1> N_days_predict[N_site_years_predict];
  
  // predict temps (forecast)
  int<lower=1> N_forecast_temp_predict; // unobserved temps
  real forecast_mean_w_no_obs_temp[N_forecast_temp_predict]; 
  real forecast_sd_w_no_obs_temp[N_forecast_temp_predict]; 
  int<lower=1, upper=N_forecast_temp_predict> temp_start_idxs_predict_forecast[N_site_years_predict];
  int<lower=1> N_days_predict_forecast[N_site_years_predict];
  
  // all obs temps for historical events with forecast
  int<lower=1> N_obs_temp_w_forecast; // unobserved temps
  real obs_temp_w_forecast[N_obs_temp_w_forecast];
  real forecast_mean_w_obs_temp[N_obs_temp_w_forecast]; 
  real forecast_sd_w_obs_temp[N_obs_temp_w_forecast]; 
  
  // all obs temps for predict events with forecast
  int<lower=1> N_pred_temp_w_forecast; // unobserved temps
  real pred_temp_w_forecast[N_pred_temp_w_forecast];
  real forecast_mean_w_pred_temp[N_pred_temp_w_forecast]; 
  real forecast_sd_w_pred_temp[N_pred_temp_w_forecast]; 

}

parameters {
  // Forcing/Chilling Parameters
  real Z_chill; // Maximum temperature for non-zero chill days
  real Z_forcing; // Minimum temperature for non-zero forcing (C)
  real<lower=0> alpha; // intercept for phenology HR threshold
  real beta1; // HR linear coefficient
  real beta2; // HR quadratic coefficient
  real<lower=0> sigma; // Overall scale gamma^{-1} * alpha (Forcing Units)
  

  // Temperature model priors
  real delta; // mean forecast bias
  real<lower=0> tau; // baseline epistemic sd (independent of forecast ensemble uncertainty)
}

model {
  // Prior model - Phenology
  Z_chill ~ normal(mean_Z_chill, sd_Z_chill);
  Z_forcing ~ normal(mean_Z_forcing, sd_Z_forcing);
  alpha ~ normal(mean_alpha, sd_alpha);
  beta1 ~ normal(mean_beta1, sd_beta1);
  beta2 ~ normal(mean_beta2, sd_beta2);
  sigma ~ normal(0, sd_sigma); 
  
  // Prior model - Temperature forecasts
  delta ~ normal(mean_delta, sd_delta);
  tau ~ normal(0, sd_tau);
  
  vector[N_obs_temp_w_forecast] obs_temp_means = delta + to_vector(forecast_mean_w_obs_temp);
  vector[N_obs_temp_w_forecast] obs_temp_sds = sqrt(tau^2 + square(to_vector(forecast_sd_w_obs_temp)));
  target += normal_lpdf(to_vector(obs_temp_w_forecast) | obs_temp_means, obs_temp_sds);

  vector[N_pred_temp_w_forecast] pred_temp_means = delta + to_vector(forecast_mean_w_pred_temp);
  vector[N_pred_temp_w_forecast] pred_temp_sds = sqrt(tau^2 + square(to_vector(forecast_sd_w_pred_temp)));
  target += normal_lpdf(to_vector(pred_temp_w_forecast) | pred_temp_means, pred_temp_sds);

  
  
  // Phenology event modeling
  real k; 
  k = 1;
  
  for (i in 1:N_site_years) {
    // Extract temperature data
    int temp_start_idx = temp_start_idxs[i];
    int temp_end_idx = temp_start_idxs[i] + N_days[i] - 1;
    
    real local_temps[N_days[i]] = obs_temps[temp_start_idx:temp_end_idx];
      
    // Extract event data
    int event_start_idx = event_start_idxs[i];
    int event_end_idx = event_start_idxs[i] + N_events[i] - 1;

    int local_precursor_event = precursor_events[i]; 
    int local_event = events[i];
    
    // Compute needed daily forcings
    real sum_chill_days;
    int first_precursor_event = local_precursor_event;
    int last_event = local_event;
    vector[last_event] daily_forcings;
    real Psi0;
    real x_minus_x0;
    
    
    sum_chill_days = 0;
    
    for (n in 1:first_precursor_event) {
      sum_chill_days += chill_days(local_temps[n], Z_chill);
    }

    for (n in first_precursor_event:last_event) {
      daily_forcings[n] = forcing(local_temps[n], Z_forcing, k);
    }
    
    // Centered chill accumulation (CA)
    x_minus_x0 = (sum_chill_days - x0_chill_days); 
    
    // Heat Requirement
    Psi0 = alpha + beta1*x_minus_x0 + beta2*x_minus_x0^2;
  
    // Increment density for each event
    real Psi = sum(daily_forcings[local_precursor_event:local_event]);
    // Log of derivative of aggregated forcings
    real log_dPsidt = log_forcing(obs_temps[local_event], 
                                  Z_forcing, k);
    
    // Add event density to model
    target += logistic_lpdf(Psi | Psi0, sigma * sqrt(3) / pi()) + log_dPsidt;
  }
}

generated quantities {
  vector[N_site_years] sum_chill_days_vector;
  int<lower=1, upper=366> pred_events[N];
  real k;
  k = 1;
  
  // Temperature forecasts
  real mod_pred_temp[N_obs_temp_w_forecast + N_pred_temp_w_forecast + N_forecast_temp_predict];
  real mod_pred_temp_sd[N_obs_temp_w_forecast + N_pred_temp_w_forecast + N_forecast_temp_predict];
  real mod_pred_temp_forecast[N_forecast_temp_predict];

  // obs temp w/ forecast
  for (i in 1:N_obs_temp_w_forecast) {
    mod_pred_temp_sd[i] = sqrt(tau^2 + forecast_sd_w_obs_temp[i]^2);
    mod_pred_temp[i] = normal_rng(delta + forecast_mean_w_obs_temp[i], mod_pred_temp_sd[i]);
  }
  
  // pred (obs.) temp w/ forecast
  for (i in 1:N_pred_temp_w_forecast) {
    int j = i + N_obs_temp_w_forecast;
    mod_pred_temp_sd[j] = sqrt(tau^2 + forecast_sd_w_pred_temp[i]^2);
    mod_pred_temp[j] = normal_rng(delta + forecast_mean_w_pred_temp[i], mod_pred_temp_sd[j]);
  }
  
  // forecast temp (no obs temp)
  for (i in 1:N_forecast_temp_predict) {
    int j = i + N_obs_temp_w_forecast + N_pred_temp_w_forecast;
    mod_pred_temp_sd[j] = sqrt(tau^2 + forecast_sd_w_no_obs_temp[i]^2);
    mod_pred_temp[j] = normal_rng(delta + forecast_mean_w_no_obs_temp[i], mod_pred_temp_sd[j]);
    mod_pred_temp_forecast[i] = mod_pred_temp[j];
  }
  
  
  // Observed phenology events
  for (i in 1:N_site_years) {
    // Extract temperature data
    int temp_start_idx = temp_start_idxs[i];
    int temp_end_idx = temp_start_idxs[i] + N_days[i] - 1;
    
    real local_temps[N_days[i]] = obs_temps[temp_start_idx:temp_end_idx];
      
    // Extract event data
    int event_start_idx = event_start_idxs[i];
    int event_end_idx = event_start_idxs[i] + N_events[i] - 1;
    
    // int local_precursor_events[N_events[i]]
    //   = precursor_events[event_start_idx:event_end_idx];
    int local_precursor_event = precursor_events[i];
    
    // Compute needed daily forcings
    int first_precursor_event = local_precursor_event; // min(local_precursor_events);
    vector[N_days[i]] daily_forcings = rep_vector(0, N_days[i]);
    vector[N_days[i]] accum_forcings;
    real sum_chill_days;
    real Psi0;  // Phenology threshold (Forcing Units)
    real x_minus_x0;
    
    sum_chill_days = 0;
    
    for (n in 1:first_precursor_event) {
      sum_chill_days += chill_days(local_temps[n], Z_chill);
    }
    
    sum_chill_days_vector[i] = sum_chill_days;

    for (n in first_precursor_event:N_days[i]) {
      daily_forcings[n] = forcing(local_temps[n], Z_forcing, k);
    }
    accum_forcings = cumulative_sum(daily_forcings);
    
    x_minus_x0 = (sum_chill_days - x0_chill_days);
    // threshold 
    Psi0 = alpha + beta1*x_minus_x0 + beta2*x_minus_x0^2;
  
    // Simulate retrodictions for observed events
    for(n in 1:N_events[i]) {
      int event_idx = event_start_idxs[i] + n - 1;
      int precursor = local_precursor_event; 
      real init_forcing = accum_forcings[precursor];
      
      real l = logistic_rng(Psi0, sigma * sqrt(3) / pi());
      
      pred_events[event_idx] = N_days[i];
      for (d in precursor:N_days[i]) {
        real Psi = accum_forcings[d] - init_forcing;
        if (Psi >= l) {
          pred_events[event_idx] = d;
          break;
        }
      }
    }
  }
  
  
  
  
  
  // Forecasts
  int<lower=1, upper=366> pred_events_forecast[N_site_years_predict];
  
  for (i in 1:N_site_years_predict) {
    //print("i = ", i);
    // Extract temperature data
    int temp_start_idx_predict = temp_start_idxs_predict[i];
    int temp_end_idx_predict = temp_start_idxs_predict[i] + N_days_predict[i] - 1;
    real local_temps_predict[N_days_predict[i]] = pred_temp_w_forecast[temp_start_idx_predict:temp_end_idx_predict];

    int temp_start_idx_predict_forecast = temp_start_idxs_predict_forecast[i];
    int temp_end_idx_predict_forecast = temp_start_idxs_predict_forecast[i] + N_days_predict_forecast[i] - 1;
    real local_temps_predict_forecast[N_days_predict_forecast[i]] = mod_pred_temp_forecast[temp_start_idx_predict_forecast:temp_end_idx_predict_forecast];

    
    // Compute needed daily forcings
    int first_precursor_event = precursor_events_predict[i];
    vector[N_days_predict[i] + N_days_predict_forecast[i]] daily_forcings = rep_vector(0, N_days_predict[i] + N_days_predict_forecast[i]); // N_days[i]);
    vector[N_days_predict[i] + N_days_predict_forecast[i]] accum_forcings = rep_vector(0, N_days_predict[i] + N_days_predict_forecast[i]);
    real sum_chill_days;
    
    real Psi0;  // Phenology threshold (Forcing Units)
    real x_minus_x0;

    for (n in first_precursor_event:N_days_predict[i]) {
      daily_forcings[n] = forcing(local_temps_predict[n], Z_forcing, k);
    }
    
    for (n in 1:N_days_predict_forecast[i]) {
      int m = n + N_days_predict[i] ;
      daily_forcings[m] = forcing(local_temps_predict_forecast[n], Z_forcing, k);
    }
    
    sum_chill_days = 0;

    for (n in 1:first_precursor_event) {
      real chill_days_n = chill_days(local_temps_predict[n], Z_chill);
      sum_chill_days += chill_days_n;
    }
    
    accum_forcings = cumulative_sum(daily_forcings);
    
    x_minus_x0 = (sum_chill_days - x0_chill_days);
    
    // threshold 
    Psi0 = alpha + beta1*x_minus_x0 + beta2*x_minus_x0^2;

    // Simulate retrodictions for observed events
    int precursor = first_precursor_event;
    real init_forcing = accum_forcings[precursor];
    
    real l = logistic_rng(Psi0, sigma * sqrt(3) / pi());
    //real l = normal_rng(Psi0, sigma);
    
    pred_events_forecast[i] = N_days_predict[i] + N_days_predict_forecast[i];
    for (d in precursor:(N_days_predict[i] + N_days_predict_forecast[i])) {
      //print("d = ", d);
      real Psi = accum_forcings[d] - init_forcing;
      if (Psi >= l) {
        pred_events_forecast[i] = d;
        // print("Psi = ", l, " exceeded!");
        break;
      }
    }
  }
  
  
  
  
  
  
  
  
}
