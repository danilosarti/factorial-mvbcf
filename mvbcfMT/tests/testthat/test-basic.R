test_that("build_components enumerates the factorial lattice", {
  Z <- matrix(c(0,1,1,0, 0,0,1,1), ncol = 2,
              dimnames = list(NULL, c("Z1", "Z2")))
  comp <- build_components(Z, c("Z1", "Z2"), order = 2L)
  expect_equal(comp$labels, c("mu", "tau_Z1", "tau_Z2", "tau_Z1:Z2"))
  expect_true(comp$is_prognostic[1])
  # interaction indicator is the product Z1 * Z2
  expect_equal(comp$indicator[, 4], Z[, 1] * Z[, 2])
  expect_equal(ncol(comp$indicator), length(comp$labels))
})

test_that("simulate_multi returns a coherent design", {
  sim <- simulate_multi(n = 120, p = 5, q = 2, Ti = 2, seed = 1)
  expect_true(all(sim$responses %in% names(sim$data)))
  expect_true(all(sim$treatments %in% names(sim$data)))
  expect_true(all(sim$covariates %in% names(sim$data)))
  expect_equal(nrow(sim$data), 120L)
  Z <- as.matrix(sim$data[, sim$treatments])
  expect_true(all(Z %in% c(0, 1)))
})

test_that("sensitivity_effect flags interval sign correctly", {
  s1 <- sensitivity_effect(estimate = -0.45, lo = -0.58, hi = -0.34, outcome_sd = 1)
  expect_true(s1$excludes_zero)
  expect_gt(s1$robustness_sd, 0)
  s2 <- sensitivity_effect(estimate = 0.02, lo = -0.05, hi = 0.09)
  expect_false(s2$excludes_zero)
  expect_equal(s2$bias_to_reach_CI_null, 0)
})

test_that("end-to-end fit, effects and diagnostics run", {
  skip_on_cran()
  sim <- simulate_multi(n = 160, p = 5, q = 2, Ti = 2, seed = 3)
  fit <- fit_mvbcf_multi(
    sim$data, responses = sim$responses, treatments = sim$treatments,
    covariates = sim$covariates,
    n_iter = 80, n_burn = 40, n_tree = 25, n_tree_tau = 15
  )
  expect_s3_class(fit, "mvbcfMT")
  eff <- component_effect(fit, "tau_Z1", test = FALSE)
  expect_equal(nrow(eff$mean), nrow(sim$data))
  expect_equal(ncol(eff$mean), length(sim$responses))

  ov <- overlap_diagnostic(fit)
  expect_s3_class(ov, "mvbcf_overlap")
  expect_true(all(c("p_min", "p_max") %in% names(ov$treatments)))

  db <- debiased_ate(fit, propensity = "glm", folds = 3)
  expect_s3_class(db, "mvbcf_debiased")
  expect_true(all(c("raw_plugin", "debiased", "se", "lo", "hi") %in% names(db)))
})
