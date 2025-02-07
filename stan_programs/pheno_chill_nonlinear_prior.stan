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
  
  // // Forcing Function
  // real forcing(real T, real Z) {
  //   //return log(1 + exp(T - T_c));
  //   return log1p(exp(T - Z));
  // }
  // 
  // real log_forcing(real T, real Z) {
  //   return log(log1p(exp(T - Z)));
  // }
  
  // Chilling function
  real chill_days(real T, real Z) {
    if (Z - T >= 0) {
      return 1 / (1 + exp(T - Z));
    } else {
      return exp(Z - T) / (1 + exp(Z - T));
    }
  }
  
  // Log of chilling function
  real log_chill_days(real T, real Z) {
    return -log1p_exp(T - Z);
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
  
  real mean_Z_chill;
  real<lower=0> sd_Z_chill;

  real mean_Z_forcing;
  real<lower=0> sd_Z_forcing;
  
  real<lower=0> mean_alpha;
  real<lower=0> sd_alpha;
  
  real mean_beta;
  real<lower=0> sd_beta;
  
  real<lower=0> sd_sigma;
  
  // real Z_chill_fixed;
}


generated quantities {
  // Forcing Function Parameters
  real Z_chill; // Maximum temperature for non-zero chill days

  real Z_forcing; // Minimum temperature for non-zero forcing (C)
  real<lower=0> alpha; // intercept for phenology threshold
  real alpha_init;
  real beta; // coefficient of chill day fraction covariate for phenology threshold

  // Phenological Transition Parameters
  real Psi0;  // Phenology threshold (Forcing Units)
  real<lower=0> sigma; // Overall scale gamma^{-1} * alpha (Forcing Units)


  alpha_init = normal_rng(mean_alpha, sd_alpha);
  while (alpha_init < 0) {
    alpha_init = normal_rng(mean_alpha, sd_alpha);
  }
  
  alpha = alpha_init;
  beta =  normal_rng(mean_beta, sd_beta);
  Z_forcing = normal_rng(mean_Z_forcing, sd_Z_forcing);
  Z_chill = normal_rng(mean_Z_chill, sd_Z_chill);
  sigma = abs(normal_rng(0, sd_sigma));     // 0   <~ sigma (FU) <~ 10

  int<lower=1, upper=366> pred_events[N];
  
  array[N_site_years, 366] real daily_chill_matrix;
  daily_chill_matrix = rep_array(0, N_site_years, 366);
  vector[N_site_years] sum_chill_days_vector;

  array[N_site_years, 366] real accum_forcings_matrix; 
  array[N_site_years, 366] real daily_forcings_matrix;

  for (i in 1:N_site_years) {
    //print("i = ", i);
    // Extract temperature data
    int temp_start_idx = temp_start_idxs[i];
    int temp_end_idx = temp_start_idxs[i] + N_days[i] - 1;
    real k;
    
    real local_temps[N_days[i]] = obs_temps[temp_start_idx:temp_end_idx];
      
    // Extract event data
    int event_start_idx = event_start_idxs[i];
    int event_end_idx = event_start_idxs[i] + N_events[i] - 1;
    
    int local_precursor_events[N_events[i]]
      = precursor_events[event_start_idx:event_end_idx];
    
    k = 1;
    // Compute needed daily forcings
    int first_precursor_event = min(local_precursor_events);
    vector[366] daily_forcings = rep_vector(0, 366); // N_days[i]);
    vector[366] accum_forcings = rep_vector(0, 366);
    real sum_chill_days;

    for (n in first_precursor_event:N_days[i]) {
      daily_forcings[n] = forcing(local_temps[n], Z_forcing, k);
    }
    
    sum_chill_days = 0;

    for (n in 1:first_precursor_event) {
      real chill_days_n = chill_days(local_temps[n], Z_chill);
      daily_chill_matrix[i, n] = chill_days_n;
      sum_chill_days += chill_days_n;
    }
    
    sum_chill_days_vector[i] = sum_chill_days;

    accum_forcings = cumulative_sum(daily_forcings);
    
    for (n in 1:366) {
      daily_forcings_matrix[i,n] = daily_forcings[n];
      accum_forcings_matrix[i,n] = accum_forcings[n];
    }
    
    // Psi0 = fma(beta, sum_chill_days - 30, alpha);  // threshold(alpha, psi1, psi2, chill); 
    Psi0 = fma(beta, 0, alpha); // use for setting priors on alpha and beta


    // Simulate retrodictions for observed events
    for(n in 1:N_events[i]) {
      int event_idx = event_start_idxs[i] + n - 1;
      int precursor = local_precursor_events[n];
      real init_forcing = accum_forcings[precursor];
      
      real l = logistic_rng(Psi0, sigma * sqrt(3) / pi());
      
      pred_events[event_idx] = N_days[i];
      for (d in precursor:N_days[i]) {
        //print("d = ", d);
        real Psi = accum_forcings[d] - init_forcing;
        if (Psi >= l) {
          pred_events[event_idx] = d;
          // print("Psi = ", l, " exceeded!");
          break;
        }
      }
    }
  }
}
