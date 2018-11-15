// Stan model for estimation of temporal changes

data {
  int < lower = 1 > N; // Number of Bones and teeth
  int < lower = 1 > T; // Number of time intervals
  vector[N] y_mean; // Means of Normal distributions of isotopic values
  vector[N] y_sigma; // SD of Normal distributions of isotopic values
  real t[T]; // index indicating the time
  matrix[N, T] x; // predictor matrix containing the renewal percentage for each interval and bone
  
  // Hyperparameters
  int < lower = 0 > mu_df;
  real mu_mean;
  real < lower = 0 > mu_sd;
  real rho_mean;
  real < lower = 0 > rho_sd;
  real alpha_mean;
  real < lower = 0 > alpha_sd;
}

parameters {
  real mu;
  real<lower=0> alpha;
  real<lower=0> rho;
  vector[T] interval; // regression coefficients (estimators of isotopic values for each interval)
}

model {
  // priors
  mu ~ student_t(1,0,10); // df, mean, sd
  rho ~ normal(1,0.25);
  alpha ~ normal(2, 0.5);
  interval ~ multi_normal_cholesky(rep_vector(mu, T), 
                                   cholesky_decompose(cov_exp_quad(t, alpha, rho)));
  // likelihood
  y_mean ~ normal(x * interval , y_sigma); 
}
