# Load Rstan
library("rstan")

# Set a timer
tmp <- Sys.time()

# Load data
dat <- read.csv("./final/data1.txt", header = T, sep = '\t')

# Parallelize
options(mc.cores = parallel::detectCores())

# Create a binary indicator for control samples (1 if enhancerID or guideID == -1, 0 otherwise)
dat$is_control <- ifelse(dat$enhancerID == -1 | dat$guideID == -1, 1, 0)

## Indicate control data
y_control <- dat[dat$is_control == 1,]$expression
y_experiment <- dat[dat$is_control == 0,]$expression

# Make sure we have the variables factored, then converted to integer representations
for (nm in colnames(dat)[colnames(dat) != "expression"]) dat[[nm]] <- as.integer(as.factor(dat[[nm]]))

# Take a look
print("nrow, ncol, head, str:")
print(nrow(dat))
print(ncol(dat))
print(head(dat))
print(str(dat))


# Prepare data for Stan
stan_data <- list(
  N = nrow(dat),
  X = dat$expression,
  gene = dat$geneID,
  enhancer = dat$enhancerID,
  guideID = dat$guideID,
  G = length(unique(dat$guideID)),  # Number of unique guides
  E = length(unique(dat$enhancerID)),
  NGENE = length(unique(dat$geneID)),
  is_control = dat$is_control,
  y_control = y_control,
  NCONTROL = length(y_control),
  NEXPERIMENT = length(y_experiment),
  y_experiment = y_experiment
)

stan_code <- "
data {
  int<lower=0> NCONTROL;                 // Number of control samples
  int<lower=0> NEXPERIMENT;              // Number of experiment samples
  int<lower=0> N;                     // Total number of observations
  int<lower=1> G;                     // Number of unique guides
  int<lower=1> E;                     // Number of unique enhancers
  real<lower=0> X[N];                 // Expression counts
  int<lower=1,upper=G> guideID[N];    // GuideID for each observation
  int<lower=1,upper=E> enhancer[N];   // EnhancerID for each observation
  int<lower=0> is_control[N];         // Control indicator (1=control, 2=experiment)
  int<lower=0> y_control[NCONTROL];             // Control samples
  int<lower=0> y_experiment[NEXPERIMENT];          // Experiment samples
}

transformed data {
  real<lower=0> mu_control_observed = mean(y_control);            // Empirical mean of y_control
  real<lower=0> sigma_control_observed = sd(y_control);           // Empirical standard deviation of y_control
  real<lower=0> mu_experiment_observed = mean(y_experiment);      // Empirical mean of y_experiment
  real<lower=0> sigma_experiment_observed = sd(y_experiment);     // Empirical standard deviation of y_experiment
  real beta_observed = mu_experiment_observed / mu_control_observed;  // Empirical effect size
}

parameters {
  real<lower=0,upper=1> r;            // Mixture probability for experimental component
  vector<lower=0,upper=1>[E] beta;    // Beta coefficients for enhancers
  real mu_experimental;            // Means for the experimental component by guideID
  real<lower=0> sigma_experimental;   // Shared experimental component standard deviation
  real mu_control[G];                    // Mean for the control component
  real<lower=0> sigma_control;        // Control standard deviation
}

model {
  // Priors
  r ~ beta(2, 2);                     // Beta Prior for r
  beta ~ beta(2, 2);                  // Priors for beta coefficients
  mu_experimental ~ normal(mu_experiment_observed, sigma_experiment_observed);    // Priors for experimental means
  sigma_experimental ~ inv_gamma(1, 1); // Prior for experimental standard deviation
  mu_control ~ normal(mu_control_observed, sigma_control_observed);         // Prior for control mean
  sigma_control ~ inv_gamma(1, 1);    // Prior for control standard deviation
  
  // Likelihood
  for (i in 1:N) {
    if (is_control[i] == 2) {
      // Control samples follow the control distribution
      X[i] ~ normal(mu_control, sigma_control);
    } else {
      // Experimental samples are a mixture
      target += log_sum_exp(
        log(r) + normal_lpdf(X[i] | mu_control[guideID[i]] * beta[enhancer[i]], sigma_experimental),
        log(1 - r) + normal_lpdf(X[i] | mu_control, sigma_control)
      );
    }
  }
}

generated quantities {
  real X_pred[N];                     // Predicted expression for each observation
  for (i in 1:N) {
    if (is_control[i] == 2) {
      // Predict from control distribution for control samples
      X_pred[i] = normal_rng(mu_control[guideID[i]], sigma_control);
    } else {
      // Predict from mixture for experimental samples
      if (bernoulli_rng(r) == 1) {
        X_pred[i] = normal_rng(mu_control[guideID[i]] * beta[enhancer[i]], sigma_experimental);
      } else {
        X_pred[i] = normal_rng(mu_control[guideID[i]], sigma_control);
      }
    }
  }
}
"

# Compile and fit the model
fit <- stan(
  model_code = stan_code,
  data = stan_data,
  iter = 1000,
  warmup = 250,
  chains = 4,
  seed = 1337,
  diagnostic_file = "FINAL_diagnostics.csv"
)
print("Time:")
print(Sys.time() - tmp)

saveRDS(fit, file = "FINAL_stan_fit.rds")
