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
    int<lower=0> N;                        // Total number of observations
    int<lower=0> NCONTROL;                 // Number of control samples
    int<lower=0> NEXPERIMENT;              // Number of experiment samples
    int<lower=1> NGENE;                    // Number of unique genes
    real X[N];                             // expression data
    int<lower=0> enhancer[N];              // Index for enhancer
    int<lower=0> is_control[N];            // Indicator for control samples
    int<lower=0> y_control[NCONTROL];             // Control samples
    int<lower=0> y_experiment[NEXPERIMENT];          // Experiment samples
    int<lower=0> E; //Enhancer ID
}

transformed data {
  real<lower=0> mu_control_observed = mean(y_control);            // Empirical mean of y_control
  real<lower=0> sigma_control_observed = sd(y_control);           // Empirical standard deviation of y_control
  real<lower=0> mu_experiment_observed = mean(y_experiment);      // Empirical mean of y_experiment
  real<lower=0> sigma_experiment_observed = sd(y_experiment);     // Empirical standard deviation of y_experiment
  real beta_observed = mu_experiment_observed / mu_control_observed;  // Empirical effect size
}

parameters {
    vector<lower=0, upper=1>[E] beta; //Beta for each enhancer
}

model {
    // Priors
    beta ~ beta(2, 2);  // Beta distribution prior

    // Likelihood
    for (i in 1:N) {
        if (is_control[i] == 2) {
            X[i] ~ normal(mu_control_observed, sigma_control_observed);
        } else {
            X[i] ~ normal(mu_control_observed*beta[enhancer[i]], sigma_experiment_observed);
        }
    }
}

generated quantities {
    real X_pred[N];  // Predicted expression values for each observation
    for (i in 1:N) {
        if (is_control[i] == 2) {
            X_pred[i] = normal_rng(mu_control_observed, sigma_control_observed);
        } else {
            X_pred[i] = normal_rng(mu_control_observed * beta[enhancer[i]], sigma_experiment_observed);
        }
    }
}

"

# Compile and fit the model
fit <- stan(
  model_code = stan_code,
  data = stan_data,
  iter = 1000,
  warmup = 200,
  chains = 4,
  diagnostic_file = "baseline_diagnostics.csv",
  cores = 4
)

print("Time:")
print(Sys.time() - tmp)

#saveRDS(fit, file = "BASELINE_stan_fit.rds")

# Extract predicted samples (posterior predictive)
posterior_predictive <- extract(fit, pars = "X_pred")$X_pred

# Compute the mean of posterior predictive samples for each observation
predicted_mean <- apply(posterior_predictive, 2, mean)

plot_data <- data.frame(
  observed = dat$expression,
  predicted = predicted_mean
)

write.csv(plot_data, "BASELINE_expression.csv", row.names = FALSE)

library(ggplot2)

plot_to_save <- ggplot(plot_data, aes(x = observed, y = predicted)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Actual vs Predicted Expression",
    x = "Observed Expression",
    y = "Predicted Expression"
  ) + theme_minimal()

ggsave(plot_to_save, file = "BASELINE_expression.png")

# Extract the posterior samples for beta
beta_samples <- extract(fit, pars = "beta")$beta

beta_mean <- apply(beta_samples, 2, mean)

beta_df <- data.frame(
  enhancerID = 1:length(beta_mean),  # Enhancer IDs (assuming 1-based indexing in Stan)
  beta_mean = beta_mean
)

# Merge beta coefficients with the original dataset
dat_beta <- merge(dat, beta_df, by.x = "enhancerID", by.y = "enhancerID", all.x = TRUE)

true_betas <- read.csv("./final/true-betas1.txt", header = T, sep = '\t')

results <- merge(dat_beta, true_betas, by = c("geneID", "enhancerID"), all.x = TRUE)

# Subset what we need
results <- results[,c("geneID", "enhancerID", "beta_mean", "beta")]
# Drop duplicate observations
results <- results[!duplicated(results),]

# Remove NAs
results <- na.omit(results)

write.csv(results, file = "BASELINE_beta.csv", row.names = F)

# Mean Squared Error
mse <- mean((results$beta - results$beta_mean)^2)

# Mean Absolute Error
mae <- mean(abs(results$beta - results$beta_mean))

cat("Mean Squared Error (MSE):", mse, "\n")
cat("Mean Absolute Error (MAE):", mae, "\n")

                   plot_to_save <- ggplot(results, aes(x = beta, y = beta_mean)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(
    title = "Actual vs Predicted Beta Values",
    x = "True Beta Values",
    y = "Predicted Beta Values"
  ) +
  theme_minimal() + ylim(c(0, 1)) + xlim(c(0, 1))

ggsave(plot_to_save, file = "BASELINE_beta.png")
