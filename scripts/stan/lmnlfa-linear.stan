data {
  int<lower=1> nobs;
  int<lower=1> p;
  int<lower=1> ni;
  int<lower=1> d;

  int<lower=1, upper=ni> person[nobs];
  int<lower=1, upper=p> itm[nobs];
  int<lower=1, upper=d> time[nobs];

  vector[nobs] age_c;        // centered age only (NO age^2)

  int y[nobs];
  int<lower=0, upper=1> is_binary[p];

  int<lower=2> k_item[p];
  int<lower=1> k_max;

  int<lower=0> nfpreds;
  int<lower=0> ntvpreds;
  matrix[nobs, nfpreds] xf;
  matrix[nobs, ntvpreds] xtv;

  matrix[p, 2] ldf;
  int<lower=0> mtv;
  int<lower=0> mf;

  real<lower=0> sigma_l;
  real<lower=0> sigma_nu;
  real<lower=0> sigma_cor;
  real<lower=0> sigma_f;
  real<lower=0> sigma_di;

  real<lower=0> total_var;
}

parameters {
  // measurement
  vector<lower=0>[p] lp;
  vector[p] np;

  vector[mf] l_diff;
  vector[mtv] l_diftv;
  vector[mf] n_diff;
  vector[mtv] n_diftv;

  // growth (LINEAR ONLY)
  real mu_slp;
  real<lower=0> phi_int;
  real<lower=0> phi_slp;

  matrix[2, nfpreds] b_mu;
  matrix[2, nfpreds] b_phi;

  vector[1] z_cor;              // only 1 correlation (int, slope)

  matrix[2, ni] fac_dist;       // 2 growth factors
  matrix[d, ni] fac_eti_raw;

  ordered[k_max - 1] tau[p];
}

transformed parameters {
  vector[p] ldiff = rep_vector(0, p);
  vector[p] ldiftv = rep_vector(0, p);
  vector[p] ndiff = rep_vector(0, p);
  vector[p] ndiftv = rep_vector(0, p);

  vector<lower=0>[2] phi_eta;
  phi_eta[1] = phi_int;
  phi_eta[2] = phi_slp;

  real<lower=0> eti_sd;
  eti_sd = sqrt(fmax(1e-6, total_var - square(phi_int)));

  matrix[d, ni] fac_eti = fac_eti_raw * eti_sd;

  // map DIF parameters
  {
    int tmp;

    tmp = 0;
    for (i in 1:p)
      if (ldf[i,2] == 1) { tmp += 1; ldiff[i] = l_diff[tmp]; }

    tmp = 0;
    for (i in 1:p)
      if (ldf[i,1] == 1) { tmp += 1; ldiftv[i] = l_diftv[tmp]; }

    tmp = 0;
    for (i in 1:p)
      if (ldf[i,2] == 1) { tmp += 1; ndiff[i] = n_diff[tmp]; }

    tmp = 0;
    for (i in 1:p)
      if (ldf[i,1] == 1) { tmp += 1; ndiftv[i] = n_diftv[tmp]; }
  }
}

model {
  // priors
  lp ~ normal(0, sigma_l);
  np ~ normal(0, sigma_nu);

  l_diff ~ normal(0, sigma_di);
  l_diftv ~ normal(0, sigma_di);
  n_diff ~ normal(0, sigma_di);
  n_diftv ~ normal(0, sigma_di);

  to_vector(b_mu) ~ normal(0, sigma_f);
  to_vector(b_phi) ~ normal(0, sigma_f);

  mu_slp ~ normal(0, sigma_f);

  phi_int ~ normal(0, sigma_f);
  phi_slp ~ normal(0, sigma_f);

  z_cor ~ normal(0, 1);

  to_vector(fac_dist) ~ normal(0, 1);
  to_vector(fac_eti_raw) ~ normal(0, 1);

  // correlation matrix (2x2)
  real rho = tanh(z_cor[1]);

  matrix[2,2] r;
  r[1,1] = 1;
  r[2,2] = 1;
  r[1,2] = rho;
  r[2,1] = rho;

  matrix[2,2] l_r = cholesky_decompose(r);

  matrix[2, ni] fac_gr;
  matrix[d, ni] fac_scor;

  // build growth factors
  {
    int k = 1;
    for (i in 1:nobs) {
      if (person[i] == k) {

        vector[2] mu_eta;
        mu_eta[1] = 0;
        mu_eta[2] = mu_slp;

        vector[2] sd_eta = phi_eta .* exp(b_phi * (xf[i,]'));

        fac_gr[,k] =
          mu_eta +
          b_mu * (xf[i,]') +
          diag_pre_multiply(sd_eta, l_r) * fac_dist[,k];

        k += 1;
      }
    }
  }

  // likelihood
  for (j in 1:nobs) {

    int ti = time[j];
    int it = itm[j];
    int pe = person[j];

    real eta_j =
      fac_gr[1, pe] +
      fac_gr[2, pe] * age_c[j] +
      fac_eti[ti, pe];

    real nu = np[it] + ndiff[it] * (nfpreds > 0 ? xf[j,1] : 0)
                        + ndiftv[it] * age_c[j];

    real lam = lp[it] + ldiff[it] * (nfpreds > 0 ? xf[j,1] : 0)
                        + ldiftv[it] * age_c[j];

    real linpred = nu + lam * eta_j;

    if (is_binary[it] == 1) {
      y[j] ~ bernoulli_logit(linpred);
    } else {
      y[j] ~ ordered_logistic(linpred,
                              head(tau[it], k_item[it] - 1));
    }
  }
}