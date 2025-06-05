library(devtools)
devtools::load_all(".", recompile = TRUE, reset = TRUE)

data <- read.csv("immdef_test.csv")

# serial
system.time(
  res_serial <- rpsftm(
    data,
    id          = "id",
    time        = "progyrs",
    event       = "prog",
    treat       = "imm",
    rx          = "rx",
    censor_time = "censyrs",
    boot        = TRUE,
    n_boot      = 100,      # small while testing
    parallel    = FALSE
  )
)

# parallel
system.time(
  res_parallel <- rpsftm(
    data,
    id          = "id",
    time        = "progyrs",
    event       = "prog",
    treat       = "imm",
    rx          = "rx",
    censor_time = "censyrs",
    boot        = TRUE,
    n_boot      = 100,
    parallel    = TRUE
  )
)

all.equal(res_serial$hr, res_parallel$hr, tolerance = 1e-2)
