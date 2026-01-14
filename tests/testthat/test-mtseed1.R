test_that("tsesimp_mt is repeatable for same seed + same threads", {
  data("shilong", package = "trtswitch")
  
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
  
  args <- list(
    data       = as.data.frame(shilong3),
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
    n_boot     = 100,
    seed       = 1L,
    threads    = 8L
  )
  
  out1 <- do.call(trtswitch::tsesimp_mt, args)
  out2 <- do.call(trtswitch::tsesimp_mt, args)
  
  # Edit these values to maybe be a bit more strict if we care about repeatability:
  tol_point <- 1e-6          
  tol_ci    <- 1e-6          
  
  expect_equal(as.numeric(out1$hr),  as.numeric(out2$hr),  tolerance = tol_point)
  expect_equal(as.numeric(out1$psi), as.numeric(out2$psi), tolerance = tol_point)
  
  expect_equal(as.numeric(out1$hr_CI),  as.numeric(out2$hr_CI),  tolerance = tol_ci)
  expect_equal(as.numeric(out1$psi_CI), as.numeric(out2$psi_CI), tolerance = tol_ci)
})
