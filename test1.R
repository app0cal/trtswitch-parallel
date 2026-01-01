# 0) Only if you added/renamed any // [[Rcpp::export]] functions:
Rcpp::compileAttributes()

# 1) Clean old compiled objects
devtools::clean_dll()

# 2) Rebuild & load (this compiles the C++)
devtools::load_all(recompile = TRUE)

library(dplyr, warn.conflicts = FALSE)
library(survival)

data("shilong", package = "trtswitch")
# build shilong3 exactly as you did before...

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

res_cpp <- tsesimpcpp_mt(data = shilong3, time = "tstop", event = "event",
                  treat = "bras.f", censor_time = "dcut", pd = "pd",
                  pd_time = "dpd", swtrt = "co", swtrt_time = "dco",
                  base_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                               "pathway.f"),
                  base2_cov = c("agerand", "sex.f", "tt_Lnum", "rmh_alea.c",
                                "pathway.f", "ps", "ttc", "tran"),
                  aft_dist = "weibull", alpha = 0.05,
                  recensor = TRUE, swtrt_control_only = FALSE, offset = 1,
                  boot = TRUE, n_boot = 500, seed = 1)

# sanity checks
res_cpp$hr
res_cpp$hr_CI
res_cpp$hr_CI_type
sd(log(res_cpp$hr_boots[!res_cpp$fail_boots]))