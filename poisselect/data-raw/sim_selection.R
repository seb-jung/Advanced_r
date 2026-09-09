# This script creates the example data set 'sim_selection' that ships with
# the package (data/sim_selection.rda). It only needs to be run again if the
# simulator or the true parameters change. Run it from the package root with
# the package loaded (devtools::load_all()), otherwise simulate_poisselect()
# is not found.

set.seed(20240501)

sim_selection <- simulate_poisselect(
  n = 800L,
  beta = c(0.5, 0.8, -0.4),
  gamma = c(0.3, 0.5, 0.7),
  sigma = 0.6,
  rho = 0.5
)

save(
  sim_selection,
  file = "data/sim_selection.rda",
  version = 3,
  compress = "xz"
)
