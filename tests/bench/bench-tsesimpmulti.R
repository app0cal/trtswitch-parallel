library(devtools)
library(dplyr, warn.conflicts = FALSE)

devtools::load_all(".")  # make sure your tsesimp_mt is loaded

cat("DEV trtswitch version:", as.character(packageVersion("trtswitch")), "\n")

# --- build shilong3 once (data prep cost is small) ---
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
  left_join(shilong2, by = c("bras.f", "id")) %>%
  as.data.frame()

common_args <- list(
  data       = shilong3,
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
  seed       = 1L
)

time_one <- function(expr) {
  gc()
  t <- system.time(force(expr))
  unname(t[["elapsed"]])
}

# --- Local timings (your fork / dev) ---
# warm-up
tmp <- common_args; tmp$threads <- 1L
invisible(do.call(trtswitch::tsesimp_mt, tmp))

t_mt1 <- time_one(do.call(trtswitch::tsesimp_mt, tmp))

tmp8 <- common_args; tmp8$threads <- 8L
t_mt8 <- time_one(do.call(trtswitch::tsesimp_mt, tmp8))

# --- CRAN official timing in a clean process ---
if (!requireNamespace("callr", quietly = TRUE)) {
  stop("Install callr to benchmark CRAN tsesimp: install.packages('callr')")
}

t_official <- callr::r(function(args) {
  library(trtswitch)
  gc()
  t <- system.time(do.call(tsesimp, args))[["elapsed"]]
  unname(t)
}, list(common_args))

cat(sprintf("elapsed CRAN tsesimp:       %.3fs\n", t_official))
cat(sprintf("elapsed DEV tsesimp_mt (1): %.3fs\n", t_mt1))
cat(sprintf("elapsed DEV tsesimp_mt (8): %.3fs\n", t_mt8))
cat(sprintf("speedup CRAN -> mt8: %.2fx\n", t_official / t_mt8))
cat(sprintf("speedup mt1 -> mt8:  %.2fx\n", t_mt1 / t_mt8))