test_that("pure llik/score/info match reference lifereg", {
  data("shilong", package = "trtswitch")
  
  # Build shilong3 exactly like your main test
  shilong1 <- shilong |>
    dplyr::arrange(bras.f, id, tstop) |>
    dplyr::group_by(bras.f, id) |>
    dplyr::slice(dplyr::n()) |>
    dplyr::select(-c("ps", "ttc", "tran"))
  
  shilong2 <- shilong |>
    dplyr::filter(pd == 0 | tstart <= dpd) |>
    dplyr::arrange(bras.f, id, tstop) |>
    dplyr::group_by(bras.f, id) |>
    dplyr::slice(dplyr::n()) |>
    dplyr::select(bras.f, id, ps, ttc, tran)
  
  shilong3 <- dplyr::left_join(shilong1, shilong2, by = c("bras.f","id"))
  
  # Which arm has at least 1 event among PD patients?
  pick_arm <- NA_integer_
  for (h in 0:1) {
    n_ev <- with(shilong3, sum(bras.f == h & pd == 1 & event == 1, na.rm = TRUE))
    if (n_ev > 0L) { pick_arm <- h; break }
  }
  if (is.na(pick_arm)) skip("No events post-progression in either arm for this dataset")
  
  sub <- subset(shilong3, bras.f == pick_arm & pd == 1)
  
  # AFT dataset used inside tsesimp: pps, event, swtrt plus covariates_aft
  covariates_aft <- c(
    "swtrt",                         # must be first
    "agerand","sex.f","tt_Lnum","rmh_alea.c","pathway.f","ps","ttc","tran"
  )
  
  data1 <- data.frame(
    pps   = sub$tstop - sub$dpd + 1, # offset = 1
    event = as.integer(sub$event),
    swtrt = as.integer(sub$co),
    agerand = sub$agerand,
    sex.f = sub$sex.f,
    tt_Lnum = sub$tt_Lnum,
    rmh_alea.c = sub$rmh_alea.c,
    pathway.f = sub$pathway.f,
    ps = sub$ps,
    ttc = sub$ttc,
    tran = sub$tran
  )
  
  # sanity: at least one event
  expect_gt(sum(data1$event == 1L, na.rm = TRUE), 0L)
  
  out <- trtswitch:::debug_aft_llik_score_info(
    data1 = data1,
    covariates_aft = covariates_aft,
    dist = "weibull",
    alpha = 0.05
  )
  
  # Compare log-likelihood at the same parameter vector
  expect_equal(out$llik_pure, out$llik_ref, tolerance = 1e-6)
  
  # Score at MLE should be near 0
  expect_lt(max(abs(unlist(out$score_pure))), 1e-4)
})