// Same as models/model_free_rev.stan, without the posterior predictions (y_rep)

functions {
  // Signum function
  real signum(real x) {return x > 0 ? 1 : x < 0 ? -1 : 0; }

  // Log-likelihood of subject j (without the shifter)
  real subject_bernoulli_logit_glm_lpmf(
    array[] int omission, array[] int singleton, array[] real reward,  
    array[] int valid, array[] int change, array[] int phase,  
    int idx_start, int I, vector alpha_svalues, vector beta_j) {
    
    // Slice all variables for subject j
    array[I] int ommision_j  = segment(omission,  idx_start, I);
    array[I] int singleton_j = segment(singleton, idx_start, I);
    array[I] real reward_j   = segment(reward,    idx_start, I);
    array[I] int valid_j     = segment(valid,     idx_start, I);
    array[I] int change_j    = segment(change,    idx_start, I);
    array[I] int phase_j     = segment(phase,     idx_start, I);
    
    // Valid trials for subject j
    int N_valid = 0;
    for (i in 1:I) if (valid_j[i] == 1) N_valid += 1;
    
    // Empty predictor matrix X
    matrix[I, 1] X;
    // Initial alpha values for subject j (same value for both singletons)
    vector[2] alpha_j = rep_vector(fmax(0.02, alpha_svalues[1]), 2);;
    // Empty vectors of valid omissions and predictors
    array[N_valid] int ommision_v_j;
    matrix[N_valid, 1] X_v;
    int row_idx = 1;

    // Compute alpha value in each trial
    for (i in 1:I) {
      // Store alpha values
      X[i, 1] = alpha_j[singleton_j[i]];

      // No shifter: nothing happens at the phase change
      if (change_j[i] == 1) {
        // no-op
      }
      
      // Store values only if trial is valid
      if (valid_j[i] == 1) {
        ommision_v_j[row_idx] = ommision_j[i]; // omission
        X_v[row_idx, 1] = (X[i, 1] - 0.02) / 0.98; // alpha
        row_idx += 1;
      }
      
      // Rescorla–Wagner-style update with phase-dependent learning
      real alpha_is = alpha_j[singleton_j[i]];
      real lr = beta_j[phase_j[i] + 2]; // 1->3; 2->4
            
      // Compute and store updated alpha values for the next trial
      real alpha_star = alpha_is + lr * abs(alpha_is) * (reward_j[i] - alpha_is);
      alpha_j[singleton_j[i]] = fmax(0.02, alpha_star);
    }
    
    // Compute and return log-likelihood
    return bernoulli_logit_glm_lpmf(ommision_v_j | X_v, beta_j[1], rep_vector(beta_j[2], 1));
  }

  // Probit exact variance in the constraint parametric space
  real probit_var_exact(real mu, real sigma) {
    real s2 = square(sigma);
    real xi  = mu / sqrt(1 + s2);
    real r = inv_sqrt(1 + 2 * s2);
    real E2 = Phi(xi) - 2 * owens_t(xi, r);
    return fmax(1e-12, E2 - square(Phi(xi)));
  }
  // Bivariate normal CDF via Owen's T with quadrant-correct Owen identity.
  real phi2_from_owen(real xi1, real xi2, real rho_eff) {
    // Edge and symmetric cases in rho
    if (rho_eff == 0) return Phi(xi1) * Phi(xi2);
    if (rho_eff ==  1) return fmin(Phi(xi1), Phi(xi2));
    if (rho_eff == -1) return fmax(Phi(xi1) + Phi(xi2) - 1, 0);
    if (xi1 == 0 && xi2 == 0) return 0.25 + asin(rho_eff) / (2 * pi());
    
    // protect from division by zero when xi==0 or xi==0
    real gamma1 = (xi1 == 0) ? 0 : (rho_eff * xi1 - xi2) / (abs(xi1) * sqrt(1 - rho_eff^2));
    real gamma2 = (xi2 == 0) ? 0 : (rho_eff * xi2 - xi1) / (abs(xi2) * sqrt(1 - rho_eff^2));
    
    // Multivariate CDF value
    real z = 0.5 * (Phi(xi1) + Phi(xi2))
             + signum(xi1) * owens_t(abs(xi1), gamma1)
             + signum(xi2) * owens_t(abs(xi2), gamma2);

    // Quadrant correction: subtract 1/2 if (x,y) in opposite signs
    // (or on an axis with negative sum)
    if (xi1 * xi2 < 0 || (xi1 * xi2 == 0 && xi1 + xi2 < 0)) z -= 0.5;
    return z;
}

// Exact covariance on probit-constrained scale for two jointly normal latents:
// U = Phi(X), V = Phi(Y), where X ~ N(mu_x, sig_x^2), Y ~ N(mu_y, sig_y^2),
// Corr(X, Y) = rho_latent. Returns Cov(U, V).
real probit_cov_exact(real mu_x, real sig_x,
                      real mu_y, real sig_y,
                      real rho_latent) {
  real xi1 = mu_x / sqrt(1 + square(sig_x));
  real xi2 = mu_y / sqrt(1 + square(sig_y));

  // "Effective" correlation after adding an independent N(0,1) for the probit trick
  real rho_eff = (sig_x * sig_y * rho_latent) /
                 sqrt((1 + square(sig_x)) * (1 + square(sig_y)));

  return phi2_from_owen(xi1, xi2, rho_eff) - Phi(xi1) * Phi(xi2);
}

}

data {
  int <lower=0> I;
  int <lower=0> J;
  int <lower=0> N;
  array[J] int <lower=1, upper=N> idx_start;
  array[N] int <lower=1, upper=2> singleton;
  array[N] real<lower=0, upper=1> reward;
  array[N] int <lower=0, upper=1> valid;
  array[N] int <lower=0, upper=1> change;
  array[N] int <lower=1, upper=2> phase;
  array[N] int <lower=0, upper=1> Omission;
}

parameters {
  // 1: intercept, 2: slope+, 3–4: LR1/LR2 (correlated), 5: alpha0
  vector[5] mu;
  vector<lower=0>[5] sigma;
  cholesky_factor_corr[2] L_phi;   // LR1–LR2 correlation
  matrix[5, J] Z_beta;
}

transformed parameters {
  matrix[J, 5] beta_j;

  // 1) intercept
  beta_j[,1] = mu[1] + (Z_beta[1,] .* sigma[1])';

  // 2) positive slope, bounded (0,10)
  beta_j[,2] = 10 * Phi_approx(mu[2] + (Z_beta[2,] .* sigma[2])');

  // 3–4) correlated LRs (bivariate latent → Phi_approx)
  beta_j[,3:4] = Phi_approx(
    rep_matrix(mu[3:4]', J) +
    (diag_pre_multiply(sigma[3:4], L_phi) * Z_beta[3:4,])'
  );

  // 5) hierarchical initial alpha (0,1)
  beta_j[,5] = Phi_approx(mu[5] + (Z_beta[5,] .* sigma[5])');
}

model {
  mu    ~ std_normal();
  sigma ~ std_normal();
  L_phi ~ lkj_corr_cholesky(1);
  to_vector(Z_beta) ~ std_normal();
  
  for (j in 1:J) {
    target += subject_bernoulli_logit_glm_lpmf(
      Omission | singleton, reward, valid, change, phase,
      idx_start[j], I,
      rep_vector(beta_j[j,5], 1),   // alpha0
      beta_j[j,1:4]'                // intercept, slope+, LR1, LR2
    );
  }
}

generated quantities {
  // Latent LR1–LR2 correlation
  corr_matrix[2] rho_lr_raw = multiply_lower_tri_self_transpose(L_phi);

  vector[5] mu_beta;
  mu_beta[1] = mu[1];
  mu_beta[2] = 10 * Phi(mu[2] / sqrt(1 + sigma[2]^2));
  mu_beta[3] = Phi(mu[3] / sqrt(1 + sigma[3]^2));
  mu_beta[4] = Phi(mu[4] / sqrt(1 + sigma[4]^2));
  mu_beta[5] = Phi(mu[5] / sqrt(1 + sigma[5]^2)); // alpha0
  
  vector[5] sigma_beta;
  sigma_beta[1] = sigma[1];
  sigma_beta[2] = 10 * sqrt(probit_var_exact(mu[2], sigma[2]));
  sigma_beta[3] = sqrt(probit_var_exact(mu[3], sigma[3]));
  sigma_beta[4] = sqrt(probit_var_exact(mu[4], sigma[4]));
  sigma_beta[5] = sqrt(probit_var_exact(mu[5], sigma[5])); // alpha0

  // LR1–LR2 correlation on the transformed (probit) scale
  real rho_lr = probit_cov_exact(
    mu[3], sigma[3], mu[4], sigma[4], rho_lr_raw[1,2]
  ) / (sigma_beta[3] * sigma_beta[4]);

  // Log-likelihood of each valid trial
  array[to_int(sum(valid))] real log_lik;

  {
    int k = 1;
    for (j in 1:J) {
      array[I] int  omission_j  = segment(Omission,  idx_start[j], I);
      array[I] int  singleton_j = segment(singleton, idx_start[j], I);
      array[I] real reward_j    = segment(reward,    idx_start[j], I);
      array[I] int  valid_j     = segment(valid,     idx_start[j], I);
      array[I] int  change_j    = segment(change,    idx_start[j], I);
      array[I] int  phase_j     = segment(phase,     idx_start[j], I);

      vector[2] alpha_j = rep_vector(fmax(0.02, beta_j[j,5]), 2);

      for (i in 1:I) {
        real X1 = alpha_j[singleton_j[i]];

        // No shifter: nothing happens at the phase change
        if (change_j[i] == 1) {
          // no-op
        }

        if (valid_j[i] == 1) {
          real X1_std = (X1 - 0.02) / 0.98; 
          real eta    = beta_j[j,1] + X1_std * beta_j[j,2];
          log_lik[k]  = bernoulli_logit_lpmf(omission_j[i] | eta);
          k += 1;
        }

        // Rescorla–Wagner-style update with phase-dependent learning
        {
          real alpha_is  = alpha_j[singleton_j[i]];
          real lr        = beta_j[j, phase_j[i] + 2]; // 1->3, 2->4
          real alpha_star = alpha_is + lr * abs(alpha_is) * (reward_j[i] - alpha_is);
          alpha_j[singleton_j[i]] = fmax(0.02, alpha_star);
        }
      }
    }
  }
}

