#`dog
# this is to test for accuracy against local multi threadding pure cpp version 
# vs the offiical library version
###

library(testthat)
library(dplyr, warn.conflicts = FALSE)
library(survival)

test_that("tsesimp (weibull AFT): local tsesimp_mt matches CRAN tsesimp", {
  # need callr to run the official one in a clean R process
  if (!requireNamespace("callr", quietly = TRUE)) {
    skip("callr not installed")
  }
  
  # make sure your dev wrapper exists (exported or not)
  if (!exists("tsesimp_mt", where = asNamespace("trtswitch"), inherits = FALSE)) {
    skip("trtswitch::tsesimp_mt not found; build/export your wrapper first")
  }
  
  # data used by the package’s own tests / vignettes
  data("shilong", package = "trtswitch")
  
  # === build shilong3 exactly as in the example ===
  shilong1 <- shilong %>%
    arrange(bras.f, id, tstop) %>%
    group_by(bras.f, id) %>%
    slice(n()) %>%
    select(-c("ps", "ttc", "tran"))
  
  shilong2 <- shilong %>%
    filter(pd == 0 | tstart <= dpd) %>%
    arrange(bras.f, id, tstop) %>%
    group_by(bras.f, id) %>%
    slice(n()) %>%
    select(bras.f, id, ps, ttc, tran)
  
  shilong3 <- shilong1 %>%
    left_join(shilong2, by = c("bras.f", "id"))
  
  # common args for both
  common_args <- list(
    data       = as.data.frame(shilong3),   # plain data.frame to avoid tibble quirks
    time       = "tstop",
    event      = "event",
    treat      = "bras.f",
    censor_time= "dcut",
    pd         = "pd",
    pd_time    = "dpd",
    swtrt      = "co",
    swtrt_time = "dco",
    base_cov   = c("agerand","sex.f","tt_Lnum","rmh_alea.c","pathway.f"),
    base2_cov  = c("agerand","sex.f","tt_Lnum","rmh_alea.c","pathway.f","ps","ttc","tran"),
    aft_dist   = "weibull",
    recensor   = TRUE,
    swtrt_control_only = FALSE,
    offset     = 1,
    boot       = TRUE,
    n_boot     = 500,
    seed = 1L
  )
  
  # official result from a clean R session using the CRAN version
  official <- callr::r(function(args) {
    library(trtswitch)
    do.call(tsesimp, args)
  }, list(common_args))
  
  # local result from your dev wrapper
  local <- do.call(trtswitch::tsesimp_mt, common_args)
  
  # compare (tight tolerance first)
  tol <- 1e-4
  expect_equal(as.numeric(local$hr),     as.numeric(official$hr),     tolerance = tol)
  expect_equal(as.numeric(local$hr_CI),  as.numeric(official$hr_CI),  tolerance = tol)
  expect_equal(as.numeric(local$psi),    as.numeric(official$psi),    tolerance = tol)
  expect_equal(as.numeric(local$psi_CI), as.numeric(official$psi_CI), tolerance = tol)
})
