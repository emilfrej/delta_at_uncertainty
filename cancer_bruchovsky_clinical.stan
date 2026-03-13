// LV parameter recovery for sparse, irregular clinical PSA data.
// PSA is normalized by the first measurement so the observation model
// is directly psa_norm ~ normal(TumourSize, sigma) with no scale parameter.
// t=0 corresponds to the first PSA measurement.
// K>1 (typically 2) allows tumour to grow above the normalized initial size.
// n0 = initial TumourSize / K (fraction of carrying capacity), bounded [0,1].

functions {
  vector dpop_dt(real t,
                 vector y,
                 array[] real theta,
                 array[] int D) {
    real S  = y[1];
    real R  = y[2];
    real rS = theta[1];
    real dS = theta[2];
    real dD = theta[3];
    real rR = theta[4];
    real dR = theta[5];
    real K  = theta[6];

    // Treatment from segments: D[1]=N_segs, D[2..]=trt values; theta[7..]=start times.
    int N_segs = D[1];
    int D_t    = D[2];
    for (s in 1:N_segs) {
      if (theta[6 + s] <= t) D_t = D[1 + s];
    }

    vector[2] dy_dt;
    dy_dt[1] = rS * S * (1 - (S + R) / K) * (1 - dD * D_t) - dS * S;
    dy_dt[2] = rR * R * (1 - (S + R) / K) - dR * R;
    return dy_dt;
  }
}

data {
  int<lower=1>    N_obs;
  array[N_obs] real          t_obs;    // days from first PSA (can start at 0)
  array[N_obs] real<lower=0> psa_norm; // PSA / PSA[1] (first value = 1.0)

  int<lower=1>    N_segs;
  array[N_segs] real           t_segs;    // segment start times (shifted; can be < 0)
  array[N_segs] int<lower=0,upper=1> trt_segs;

  real<lower=0>   rS;
  real<lower=0>   K;    // carrying capacity — set > 1 to allow growth above initial size
  real<lower=0>   dD;
}

transformed data {
  array[1 + N_segs] int D_packed;
  D_packed[1] = N_segs;
  for (s in 1:N_segs) D_packed[1 + s] = trt_segs[s];
}

parameters {
  real<lower=0, upper=1> cost;
  real<lower=0, upper=1> turnover;
  real<lower=0, upper=1> n0;     // initial TumourSize / K at t=0
  real<lower=0, upper=1> rFrac;
  real<lower=0>          sigma;
}

transformed parameters {
  real rR = (1 - cost)  * rS;
  real dS = turnover    * rS;
  real dR = turnover    * rS;

  // Initial conditions: TumourSize(0) = K * n0
  vector[2] y0;
  y0[1] = K * n0 * (1 - rFrac);
  y0[2] = K * n0 * rFrac;

  array[6 + N_segs] real theta;
  theta[1] = rS;  theta[2] = dS;  theta[3] = dD;
  theta[4] = rR;  theta[5] = dR;  theta[6] = K;
  for (s in 1:N_segs) theta[6 + s] = t_segs[s];

  // t0 = -1.0 so t_obs[1] = 0 is strictly after the initial time
  array[N_obs] vector[2] cells =
    ode_rk45(dpop_dt, y0, -1.0, t_obs, theta, D_packed);

  vector[N_obs] tumour_pred;
  for (i in 1:N_obs)
    tumour_pred[i] = fmax(1e-6, cells[i][1] + cells[i][2]);
}

model {
  cost     ~ beta(1, 1);
  turnover ~ beta(2, 5);
  n0       ~ beta(1, 1);    // flat — data determines initial TumourSize
  rFrac    ~ beta(1, 30);
  sigma    ~ lognormal(-2, 1);

  psa_norm ~ normal(tumour_pred, sigma);
}

generated quantities {
  vector[N_obs] psa_pred;
  vector[N_obs] S_pred;
  vector[N_obs] R_pred;
  for (i in 1:N_obs) {
    psa_pred[i] = normal_rng(tumour_pred[i], sigma);
    S_pred[i]   = cells[i][1];
    R_pred[i]   = cells[i][2];
  }
}
