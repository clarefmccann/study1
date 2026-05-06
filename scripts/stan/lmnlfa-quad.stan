data {
  // dimensions
  int<lower=1> nobs;              // total person-occasion-item observations
  int<lower=1> p;                 // number of items
  int<lower=1> ni;                // number of individuals
  int<lower=1> d;                 // number of time points (waves)

  // indices
  int<lower=1, upper=ni> person[nobs];
  int<lower=1, upper=p> itm[nobs];
  int<lower=1, upper=d> time[nobs];

  // time metric (recommended: centered age, plus centered age^2 computed in R)
  vector[nobs] age_c;
  vector[nobs] age2_c;

  // responses
  int y[nobs];                    // holds both 0/1 and 1..k responses
  int<lower=0, upper=1> is_binary[p]; // 1 if item is 0/1, else 0

  // for ordinal items:
  int<lower=2> k_item[p];         // number of categories for each item
                                  // set k_item=2 for binary items too, harmless
  int<lower=1> k_max;             // max categories across items (e.g., 4)
                                  // ordinal thresholds will be size k_max-1

  // predictors
  int<lower=0> nfpreds;           // # time-invariant predictors
  int<lower=0> ntvpreds;          // # time-varying predictors
  matrix[ni, nfpreds] xf_person;  // person-level covariates
  matrix[nobs, nfpreds] xf;       // time-invariant covariates (race, bmi, baseline age, etc.)
  matrix[nobs, ntvpreds] xtv;     // time-varying covariates (if you use them for dif/impact)

  // dif pattern like the supplement (col1 time-varying, col2 invariant)
  matrix[p, 2] ldf;
  int<lower=0> mtv;               // count of items flagged for time-varying dif
  int<lower=0> mf;                // count of items flagged for invariant dif

  // prior scales (mirror their inputs)
  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;
}

parameters {
  // baseline loadings/intercepts (like Lp, Np)
  vector<lower=0>[p] lp;
  vector[p] np;

  // dif params (vectors that get mapped into item-length vectors via ldf)
  vector[mf] l_diff;         // loading dif invariant
  vector[mtv] l_diftv;       // loading dif time-varying (multiplies age or other tv predictor)
  vector[mf] n_diff;         // intercept dif invariant
  vector[mtv] n_diftv;       // intercept dif time-varying

  // growth: add quadratic term (nonlinear component)
  real mu_slp;               // slope mean (linear)
  real mu_quad;              // quadratic mean (nonlinear)
  real<lower=0> phi_int;
  real<lower=0> phi_slp;
  real<lower=0> phi_quad;
  real<lower=0> eti_sd;

  // impact of invariants on growth factor means/sds
  matrix[3, nfpreds] b_mu;   // mean impact for (int, slp, quad)
  matrix[3, nfpreds] b_phi;  // sd impact  for (int, slp, quad)

  // Cholesky factor for a proper correlation matrix among growth factors
  cholesky_factor_corr[3] L_Omega;

  // non-centered random effects
  matrix[3, ni] fac_dist;    // std normal draws for growth random effects
  matrix[d, ni] fac_eti_raw; // time-specific errors

  // ordinal thresholds for each item (k_max - 1 cutpoints per item)
  // for binary items we will ignore these and use bernoulli_logit
  ordered[k_max - 1] tau[p];
}

transformed parameters {
  // map dif vectors into item-length vectors (like the supplement does) :contentReference[oaicite:2]{index=2}
  vector[p] ldiff = rep_vector(0, p);
  vector[p] ldiftv = rep_vector(0, p);
  vector[p] ndiff = rep_vector(0, p);
  vector[p] ndiftv = rep_vector(0, p);

  // growth factor sd vector
  vector<lower=0>[3] phi_eta;
  phi_eta[1] = phi_int;
  phi_eta[2] = phi_slp;
  phi_eta[3] = phi_quad;

  matrix[d, ni] fac_eti = fac_eti_raw * eti_sd;

  // assign dif values into item vectors based on ldf
  {
    int tmp;

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 2] == 1) {
        tmp += 1;
        ldiff[i] = l_diff[tmp];
      }
    }

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 1] == 1) {
        tmp += 1;
        ldiftv[i] = l_diftv[tmp];
      }
    }

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 2] == 1) {
        tmp += 1;
        ndiff[i] = n_diff[tmp];
      }
    }

    tmp = 0;
    for (i in 1:p) {
      if (ldf[i, 1] == 1) {
        tmp += 1;
        ndiftv[i] = n_diftv[tmp];
      }
    }
  }
}

model {
  // priors (roughly analogous to their choices) :contentReference[oaicite:4]{index=4}
  lp ~ normal(0, sigma_l);
  np ~ normal(0, sigma_nu);
  eti_sd ~ normal(0, sigma_f);

  l_diff ~ normal(0, sigma_di);
  l_diftv ~ normal(0, sigma_di);
  n_diff ~ normal(0, sigma_di);
  n_diftv ~ normal(0, sigma_di);

  to_vector(b_mu) ~ normal(0, sigma_f);
  to_vector(b_phi) ~ normal(0, sigma_f);

  mu_slp ~ normal(0, sigma_f);
  mu_quad ~ normal(0, sigma_f);

  phi_int ~ normal(0, sigma_f);
  phi_slp ~ normal(0, sigma_f);
  phi_quad ~ normal(0, sigma_f);

  L_Omega ~ lkj_corr_cholesky(sigma_cor);
  to_vector(fac_dist) ~ normal(0, 1);
  to_vector(fac_eti_raw) ~ normal(0, 1);
  for (it in 1:p) {
    tau[it] ~ normal(0, 1.5);
  }

  // latent growth factors per person
  matrix[3, ni] fac_gr;

  // time-specific factor scores
  matrix[d, ni] fac_scor;

  // compute fac_gr using non-centered random effects + covariate impact
  // matches their logic, just expanded to 3 growth factors :contentReference[oaicite:5]{index=5}
  for (k in 1:ni) {
    vector[3] mu_eta;
    mu_eta[1] = 0;
    mu_eta[2] = mu_slp;
    mu_eta[3] = mu_quad;

    vector[3] sd_eta = phi_eta .* exp(b_phi * (xf_person[k,]'));
    fac_gr[,k] = mu_eta
               + b_mu * (xf_person[k,]')
               + diag_pre_multiply(sd_eta, L_Omega) * fac_dist[,k];
  }

  // likelihood loop
  // note: their example does fac_scor = int + slp*age + error :contentReference[oaicite:7]{index=7}
  // we add + quad*age^2
  for (j in 1:nobs) {
    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    real eta_j = fac_gr[1, pe]
               + fac_gr[2, pe] * age_c[j]
               + fac_gr[3, pe] * age2_c[j]
               + fac_eti[ti, pe];

    // dif-adjusted intercept and loading
    real nu = np[it] + ndiff[it] * (nfpreds > 0 ? xf[j,1] : 0) + ndiftv[it] * age_c[j];
    real lam = lp[it] * exp(ldiff[it] * xf[j,1] + ldiftv[it] * age_c[j]);
    real linpred = nu + lam * eta_j;

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(linpred);  // same family as their example :contentReference[oaicite:8]{index=8}
    } else {
      // For ordinal items, location is identified through thresholds, not a free intercept plus thresholds.
      y[j] ~ ordered_logistic(nu + lam * eta_j, head(tau[it], k_item[it] - 1));
    }
  }
}

generated quantities {
  // you can add derived objects here (variance matrices, etc.)
  // like they do in generated quantities :contentReference[oaicite:9]{index=9}
}

