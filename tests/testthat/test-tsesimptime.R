suppressPackageStartupMessages({
  library(dplyr, warn.conflicts = FALSE)
  library(callr)
})

# --- build shilong3 once in the parent ---
data("shilong", package = "trtswitch")

shilong1 <- shilong %>%
  arrange(bras.f, id, tstop) %>%
  group_by(bras.f, id) %>%
  slice(n()) %>%
  select(-c("ps","ttc","tran"))

shilong2 <- shilong %>%
  filter(pd == 0 | tstart <= dpd) %>%
  arrange(bras.f, id, tstop) %>%
  group_by(bras.f, id) %>%
  slice(n()) %>%
  select(bras.f, id, ps, ttc, tran)

shilong3 <- shilong1 %>% left_join(shilong2, by = c("bras.f","id"))

common_args <- list(
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
  n_boot     = 5000,
  seed       = 1L
)

# Make threads deterministic in both children
omp_env <- c(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  OMP_PLACES = "",
  OMP_PROC_BIND = ""
)

## 1) Time the CRAN/installed version (needs tsesimp exported there)
official_median <- callr::r(
  function(args, env) {
    do.call(Sys.setenv, as.list(env))
    suppressPackageStartupMessages(library(trtswitch))
    stopifnot("tsesimp" %in% getNamespaceExports("trtswitch"))
    invisible(do.call(trtswitch::tsesimp, args))  # warm-up
    times <- numeric(3L)
    for (i in seq_along(times)) {
      gc()
      t0 <- proc.time()[["elapsed"]]
      invisible(do.call(trtswitch::tsesimp, args))
      times[i] <- proc.time()[["elapsed"]] - t0
    }
    median(times)
  },
  args = list(common_args, omp_env)
)

cat(sprintf("Official tsesimp (1 thread) median: %.3f s\n", official_median))

## 2) Time your dev build (loaded via pkgload in an isolated process)
dev_path <- "C:/Users/marco/OneDrive/Desktop/trtswitch-test"

local_medians <- vapply(c(1L,2L,4L,8L), function(tn) {
  callr::r(function(args, path, env, tn) {
    env$OMP_NUM_THREADS <- as.character(tn)
    do.call(Sys.setenv, as.list(env))
    suppressPackageStartupMessages(library(pkgload))
    pkgload::load_all(path, quiet = TRUE)               # loads your dev namespace
    stopifnot("tsesimp_mt" %in% getNamespaceExports("trtswitch"))
    invisible(do.call(trtswitch::tsesimp_mt, args))     # warm-up
    times <- numeric(5L)
    for (i in seq_along(times)) {
      gc()
      t0 <- proc.time()[["elapsed"]]
      invisible(do.call(trtswitch::tsesimp_mt, args))
      times[i] <- proc.time()[["elapsed"]] - t0
    }
    median(times)
  }, args = list(common_args, dev_path, omp_env, tn))
}, numeric(1))

res <- data.frame(
  threads = c(1L,2L,4L,8L),
  median_sec = as.numeric(local_medians),
  speedup_vs_official_1t = as.numeric(official_median) / as.numeric(local_medians)
)
print(res, row.names = FALSE)
