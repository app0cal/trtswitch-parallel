#`dog
# this is to test for accuracy against local multi threadding pure cpp version 
# vs the offiical library version
###

library(testthat)
library(dplyr, warn.conflicts = FALSE)
library(survival)
official_meta <- callr::r(function() {
  library(trtswitch)
  c(
    version = as.character(packageVersion("trtswitch")),
    path    = find.package("trtswitch")
  )
})
print(official_meta)


test_that("threads = 0 uses automatic detectoin", {
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

  local_args <- common_args
  local_args$threads <- 0L
  # local result from local tsesimp version
  local <- do.call(trtswitch::tsesimp_mt, local_args)
  
  expect_true(local$threads_used >= 1L)
})
