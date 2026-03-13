library(rstan)
library(tidyverse)
library(tidybayes)
rstan_options(auto_write = TRUE)
options(mc.cores = 4)

# --- Constants ---
rS <- 0.027; K <- 1.0; dD <- 1.5; k <- 0.2
N_PRIOR <- 20000
CACHE_DIR <- "fits_cache_p12_clinical_seq"
dir.create(CACHE_DIR, showWarnings = FALSE)

# --- Helper functions ---
calc_double_time_fixed_k <- function(dR, k, K, rR, R0) {
  (1 / (dR - (1 - k) * rR)) *
    (log(1 + (rR * R0) / (dR * K - (1 - k) * K * rR + rR * R0)) - log(2))
}

calc_delta_at <- function(dR, k, K, rR, R0) {
  delta_k <- calc_double_time_fixed_k(dR, k, K, rR, R0)
  delta_0 <- calc_double_time_fixed_k(dR, 0, K, rR, R0)
  (delta_k - delta_0) / delta_0
}

# KL divergence from two sample vectors via kernel density on shared grid
kl_div <- function(p_samp, q_samp, n_grid = 1024, eps = 1e-10) {
  p_samp <- p_samp[is.finite(p_samp)]
  q_samp <- q_samp[is.finite(q_samp)]
  if (length(p_samp) < 20 || length(q_samp) < 20) return(NA_real_)
  xr <- quantile(c(p_samp, q_samp), c(0.001, 0.999))
  if (diff(xr) < 1e-10) return(NA_real_)
  p_d <- pmax(density(p_samp, from = xr[1], to = xr[2], n = n_grid)$y, eps)
  q_d <- pmax(density(q_samp, from = xr[1], to = xr[2], n = n_grid)$y, eps)
  p_d <- p_d / sum(p_d)
  q_d <- q_d / sum(q_d)
  sum(p_d * log(p_d / q_d))
}

# --- Parse data ---
raw <- read_csv(
  "Bruchovsky_et_al/patient012.txt",
  col_names = c("id", "date", "dose_mg", "testo_thresh",
                "psa", "testosterone", "cycle", "on_treatment",
                "calendar_day", "day"),
  show_col_types = FALSE
) %>% arrange(day)

psa_data <- raw %>%
  filter(!is.na(psa)) %>%
  select(day, psa, on_treatment) %>%
  mutate(psa_norm = psa / psa[1])

day0 <- psa_data$day[1]

all_segs <- raw %>%
  select(day, on_treatment) %>%
  mutate(on_treatment = as.integer(on_treatment)) %>%
  distinct(day, .keep_all = TRUE) %>%
  mutate(prev_trt = lag(on_treatment, default = 1L)) %>%
  filter(on_treatment != prev_trt | row_number() == 1) %>%
  bind_rows(tibble(day = 0L, on_treatment = 1L), .) %>%
  distinct(day, .keep_all = TRUE) %>%
  arrange(day) %>%
  mutate(t_seg = as.double(day) - day0)

t_obs_all <- as.double(psa_data$day) - day0
N_obs <- nrow(psa_data)

trt_abs <- tibble(
  day = 0:max(psa_data$day),
  trt = all_segs$on_treatment[findInterval(as.double(0:max(psa_data$day)),
                                            all_segs$day)]
) %>%
  mutate(trt = ifelse(is.na(trt), 1L, trt),
         block = cumsum(c(TRUE, diff(trt) != 0))) %>%
  filter(trt == 1) %>%
  group_by(block) %>%
  summarise(xmin = min(day) - 0.5, xmax = max(day) + 0.5, .groups = "drop")

cat("Patient 12:", N_obs, "PSA observations\n")
cat("Days (shifted):", paste(round(t_obs_all), collapse=", "), "\n")

# --- Sequential Stan fits ---
for (t in 2:N_obs) {
  cache_path <- file.path(CACHE_DIR, sprintf("p12_t%02d.rds", t))
  if (file.exists(cache_path)) {
    fit <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(fit) && length(rstan::extract(fit)$cost) > 0) next
  }

  psa_sub  <- psa_data[1:t, ]
  t_obs_t  <- t_obs_all[1:t]
  t_max    <- max(t_obs_t)

  segs_t   <- all_segs %>% filter(t_seg <= t_max)
  if (nrow(segs_t) == 0) segs_t <- all_segs[1, ]

  data_list <- list(
    N_obs    = t,
    t_obs    = as.array(t_obs_t),
    psa_norm = as.array(psa_sub$psa_norm),
    N_segs   = nrow(segs_t),
    t_segs   = as.array(segs_t$t_seg),
    trt_segs = as.array(segs_t$on_treatment),
    rS = rS, K = K, dD = dD
  )

  fit <- tryCatch(
    stan("cancer_bruchovsky_clinical.stan",
         data    = data_list,
         chains  = 4, cores  = 4,
         warmup  = 1500, iter  = 3000,
         seed    = 42,
         control = list(adapt_delta = 0.95),
         refresh = 0),
    error = function(e) NULL
  )
  if (!is.null(fit)) saveRDS(fit, cache_path)
}
