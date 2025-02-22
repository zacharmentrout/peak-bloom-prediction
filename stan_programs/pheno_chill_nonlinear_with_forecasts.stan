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
  
  // real Z_chill_fixed; 
}

parameters {
  // Forcing/Chilling Parameters
  real Z_chill; // Maximum temperature for non-zero chill days
  real Z_forcing; // Minimum temperature for non-zero forcing (C)
  real<lower=0> alpha; // intercept for phenology threshold
  real beta1; // coefficient of chill day fraction covariate for phenology threshold +/- 1ish but bias negative
  real beta2; // quadratic term 
  real<lower=0> sigma; // Overall scale gamma^{-1} * alpha (Forcing Units)
  

  // Temperature model priors
  delta ~ normal(0, 2); // mean forecast bias
  tau ~ normal(0, 2); // baseline forecast sd (independent of forecast ensemble uncertainty)
}

model {
  // Prior model
  Z_chill ~ normal(mean_Z_chill, sd_Z_chill);
  Z_forcing ~ normal(mean_Z_forcing, sd_Z_forcing);
  alpha ~ normal(mean_alpha, sd_alpha);
  beta1 ~ normal(mean_beta1, sd_beta1);
  beta2 ~ normal(mean_beta2, sd_beta2);
  sigma ~ normal(0, sd_sigma);     // 0   <~ sigma (FU) <~ 10
  
  
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
    
    int local_precursor_events[N_events[i]]
      = precursor_events[event_start_idx:event_end_idx];
    int local_events[N_events[i]]
      = events[event_start_idx:event_end_idx];
    
    // Compute needed daily forcings
    real sum_chill_days;
    int first_precursor_event = min(local_precursor_events);
    int last_event = max(local_events);
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
    
    x_minus_x0 = (sum_chill_days - x0_chill_days);
    
    Psi0 = alpha + beta1*x_minus_x0 + beta2*x_minus_x0^2;
    
  
    // Increment density for each event
    for(n in 1:N_events[i]) {
      // Aggregated forcings during phenology interval
      real Psi = sum(daily_forcings[local_precursor_events[n]:local_events[n]]);
      
      // Log of derivative of aggregated forcings
      real log_dPsidt = log_forcing(obs_temps[local_events[n]], 
                                    Z_forcing, k);
      
      // Add event density to model
      target += logistic_lpdf(Psi | Psi0, sigma * sqrt(3) / pi()) + log_dPsidt;
      //target += normal_lpdf(Psi | Psi0, sigma) + log_dPsidt;
    }
  }
}

generated quantities {
  vector[N_site_years] sum_chill_days_vector;
  int<lower=1, upper=366> pred_events[N];
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
    
    int local_precursor_events[N_events[i]]
      = precursor_events[event_start_idx:event_end_idx];
    
    // Compute needed daily forcings
    int first_precursor_event = min(local_precursor_events);
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
      int precursor = local_precursor_events[n];
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
}
