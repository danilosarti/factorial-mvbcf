#' Simulate a two-treatment multi-environment trial with known effects
#'
#' Generates a synthetic MET with TWO crossed binary treatments (e.g.
#' \code{Z1} = drought stress, \code{Z2} = a biological inoculant / N top-dressing)
#' applied in a \eqn{2\times2} factorial within each genotype x environment cell.
#' The multivariate response (yield, thousand-grain weight, height) is generated
#' from known prognostic and effect surfaces, including a genuine
#' \strong{interaction} (the inoculant buffers part of the drought loss). All
#' ground-truth surfaces are returned (columns \code{true_*}) for validation; the
#' model never sees them. Assignment is (optionally) confounded through the
#' covariates so the propensity adjustment matters.
#'
#' @param n_gen,n_env,n_check Genotypes, environments, and always-present checks.
#' @param confounded If TRUE, treatment assignment depends on covariates.
#' @param seed Random seed.
#' @return list(\code{data}, \code{genotype}, \code{environment}).
#' @examples
#' met <- simulate_met2(n_gen = 20, n_env = 5, seed = 1)
#' nrow(met$data)
#' @export
simulate_met2 <- function(n_gen = 40, n_env = 10, n_check = 6,
                          confounded = TRUE, seed = 2024) {
  set.seed(seed)
  gen <- sprintf("G%03d", seq_len(n_gen))
  geno <- data.frame(gen = gen, tolerance = stats::rbeta(n_gen, 2, 2),
                     responsiveness = stats::rbeta(n_gen, 2, 2),
                     g_yield = stats::rnorm(n_gen, 0, 500),
                     g_tgw = stats::rnorm(n_gen, 0, 3),
                     g_hgt = stats::rnorm(n_gen, 0, 6), stringsAsFactors = FALSE)
  locs <- c("Viterbo","Santaella","Cordoba","Sevilla","Marchouch","TelHadya")
  years <- c(2020, 2021, 2022)
  combos <- expand.grid(loc = locs, year = years, stringsAsFactors = FALSE)
  combos <- combos[sample(nrow(combos), n_env), ]
  E <- nrow(combos)
  env <- data.frame(env = paste(combos$loc, combos$year, sep = "_"),
                    loc = combos$loc, year = combos$year,
                    temp_mean = round(stats::runif(E, 16, 27), 2),
                    precip_total = round(stats::runif(E, 120, 520), 1),
                    vpd = round(stats::runif(E, 0.8, 2.6), 2),
                    radiation = round(stats::runif(E, 15, 26), 2),
                    rel_humidity = round(stats::runif(E, 35, 75), 1),
                    stringsAsFactors = FALSE)
  clay <- stats::runif(E, 8, 45); sand <- stats::runif(E, 15, 60)
  env$soil_clay <- round(clay, 1); env$soil_sand <- round(sand, 1)
  env$soil_silt <- round(pmax(100 - clay - sand, 5), 1)
  env$soil_pH <- round(stats::runif(E, 5.5, 8.3), 2)
  env$soil_OM <- round(stats::runif(E, 0.6, 3.5), 2)
  env$soil_N  <- round(stats::runif(E, 0.04, 0.22), 3)
  env$soil_CEC <- round(stats::runif(E, 8, 35), 1)
  sw <- 0.012 * env$soil_clay + 0.10 * env$soil_OM + stats::rnorm(E, 0, 0.05)
  env$soil_water_cap <- round((sw - min(sw)) / (max(sw) - min(sw)) * 0.8 + 0.1, 3)
  z <- function(x) (x - mean(x)) / stats::sd(x)
  sev_lin <- -0.9*z(env$precip_total) + 0.8*z(env$vpd) + 0.5*z(env$temp_mean) - 0.3*z(env$rel_humidity)
  env$severity <- round(0.2 + 0.7 * stats::plogis(sev_lin), 3)
  env$e_yield <- round(stats::rnorm(E, 0, 800), 0)
  env$e_tgw <- round(stats::rnorm(E, 0, 2.5), 2)
  env$e_hgt <- round(stats::rnorm(E, 0, 8), 2)

  checks <- gen[seq_len(n_check)]; others <- setdiff(gen, checks)
  mem <- do.call(rbind, lapply(seq_len(E), function(i) {
    present <- c(checks, sample(others, sample(max(8, n_gen%/%3):(n_gen%/%2), 1)))
    data.frame(gen = present, env = env$env[i], stringsAsFactors = FALSE)
  }))
  gmap <- geno; rownames(gmap) <- geno$gen
  emap <- env;  rownames(emap) <- env$env
  rows <- list(); r <- 0
  for (m in seq_len(nrow(mem))) {
    rg <- gmap[mem$gen[m], ]; re <- emap[mem$env[m], ]
    # ---- ground-truth effect surfaces ----
    impact <- re$severity * (1 - 0.7 * rg$tolerance) * (1 - 0.4 * re$soil_water_cap)  # drought
    boost  <- rg$responsiveness * (0.3 + 0.5 * re$soil_N * 4)                          # inoculant/N gain
    synergy<- 0.6 * impact * rg$responsiveness                                         # inoculant buffers drought
    mu_y <- 6200 + rg$g_yield + re$e_yield + 25 * re$soil_OM
    mu_t <- 42 + rg$g_tgw + re$e_tgw; mu_h <- 92 + rg$g_hgt + re$e_hgt
    t1y <- -2600*impact; t1t <- -9*impact;  t1h <- -16*impact          # drought (negative)
    t2y <-  900*boost;   t2t <-  2.5*boost; t2h <-  6*boost            # inoculant (positive)
    t12y<-  1400*synergy;t12t<-  4*synergy; t12h<- 5*synergy           # interaction (buffering)
    # ---- (optionally confounded) 2x2 assignment ----
    if (confounded) {
      pz1 <- stats::plogis(-0.4 + 1.2*(re$severity - 0.5))            # drought likelier in harsh env
      pz2 <- stats::plogis(-0.2 + 1.5*(rg$responsiveness - 0.5))      # inoculant likelier on responsive geno
    } else { pz1 <- 0.5; pz2 <- 0.5 }
    for (Z1 in c(0,1)) for (Z2 in c(0,1)) {
      # keep the design roughly factorial but let propensity tilt inclusion
      if (confounded && stats::runif(1) > (Z1*pz1 + (1-Z1)*(1-pz1)) *
                                          (Z2*pz2 + (1-Z2)*(1-pz2)) * 1.6) next
      r <- r + 1
      rows[[r]] <- data.frame(
        gen = mem$gen[m], env = mem$env[m], loc = re$loc, year = re$year,
        Z1 = Z1, Z2 = Z2,
        temp_mean = re$temp_mean, precip_total = re$precip_total, vpd = re$vpd,
        radiation = re$radiation, rel_humidity = re$rel_humidity,
        soil_clay = re$soil_clay, soil_sand = re$soil_sand, soil_silt = re$soil_silt,
        soil_pH = re$soil_pH, soil_OM = re$soil_OM, soil_N = re$soil_N,
        soil_CEC = re$soil_CEC, soil_water_cap = re$soil_water_cap,
        yield_kg_ha = round(mu_y + Z1*t1y + Z2*t2y + Z1*Z2*t12y + stats::rnorm(1,0,220), 1),
        tgw_g       = round(mu_t + Z1*t1t + Z2*t2t + Z1*Z2*t12t + stats::rnorm(1,0,1.1), 2),
        plant_height_cm = round(mu_h + Z1*t1h + Z2*t2h + Z1*Z2*t12h + stats::rnorm(1,0,3.0), 1),
        true_tau1_yield = round(t1y,1), true_tau2_yield = round(t2y,1), true_tau12_yield = round(t12y,1),
        true_tau1_tgw = round(t1t,3),   true_tau2_tgw = round(t2t,3),   true_tau12_tgw = round(t12t,3),
        true_tau1_height = round(t1h,2),true_tau2_height = round(t2h,2),true_tau12_height = round(t12h,2),
        severity = re$severity, tolerance = round(rg$tolerance,4),
        responsiveness = round(rg$responsiveness,4), stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows); df <- cbind(row_id = seq_len(nrow(df)), df)
  list(data = df, genotype = geno, environment = env)
}

#' Simulate a generic T-treatment causal dataset with known effects
#'
#' A compact, domain-neutral generator for validating the general-order engine:
#' \code{Ti} crossed binary treatments, \code{p} covariates, \code{q} outcomes,
#' known heterogeneous main effects and (selected) pairwise interactions.
#'
#' @param n,p,q,Ti Units, covariates, outcomes, treatments.
#' @param confounded Confound the assignment through the covariates.
#' @param seed Seed.
#' @return list(\code{data}, \code{truth}) where \code{truth} holds the true
#'   per-row component surfaces keyed by subset label.
#' @examples
#' sim <- simulate_multi(n = 100, p = 5, q = 2, Ti = 2, seed = 1)
#' str(sim$data[, 1:6])
#' @export
simulate_multi <- function(n = 1500, p = 6, q = 2, Ti = 3,
                           confounded = TRUE, seed = 7) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p); colnames(X) <- paste0("x", seq_len(p))
  # assignment (optionally confounded on x1,x2)
  Z <- sapply(seq_len(Ti), function(t) {
    lin <- if (confounded) 0.6 * X[, 1] * (t %% 2 == 1) - 0.5 * X[, 2] * (t %% 2 == 0) else 0
    stats::rbinom(n, 1, stats::plogis(lin - 0.1 * t))
  })
  colnames(Z) <- paste0("Z", seq_len(Ti))
  mu <- 2 + 1.5 * X[, 1] - X[, 2]^2 * 0.4 + 0.8 * X[, 3]
  # main effects: heterogeneous in different covariates
  tau_main <- lapply(seq_len(Ti), function(t)
    (1 + t) * 0.8 + 1.2 * X[, ((t - 1) %% p) + 1] * (t %% 2 == 1) +
      0.9 * (X[, ((t) %% p) + 1] > 0))
  # real pairwise interactions (only those whose treatments exist for this Ti)
  pair_true <- list()
  if (Ti >= 2) pair_true[["1:2"]] <- 1.4 * (X[, 1] > 0) * (X[, 2] > 0)
  if (Ti >= 3) pair_true[["2:3"]] <- -1.1 * X[, 3]
  q_scale <- seq(1, by = 0.6, length.out = q)   # outcome-specific scaling
  Ey <- sapply(seq_len(q), function(k) {
    s <- q_scale[k]
    val <- s * mu
    for (t in seq_len(Ti)) val <- val + Z[, t] * s * tau_main[[t]]
    if (Ti >= 2) val <- val + Z[, 1] * Z[, 2] * s * pair_true[["1:2"]]
    if (Ti >= 3) val <- val + Z[, 2] * Z[, 3] * s * pair_true[["2:3"]]
    val
  })
  Sig <- 0.5 * diag(q) + 0.3; E <- matrix(stats::rnorm(n * q), n) %*% chol(Sig)
  y <- Ey + E; colnames(y) <- paste0("y", seq_len(q))
  data <- data.frame(row_id = seq_len(n), y, Z, X)
  truth <- list(mu = outer(mu, q_scale),
                main = lapply(seq_len(Ti), function(t) outer(tau_main[[t]], q_scale)),
                pair = lapply(pair_true, function(v) outer(v, q_scale)))
  list(data = data, truth = truth,
       responses = colnames(y), treatments = colnames(Z), covariates = colnames(X))
}
