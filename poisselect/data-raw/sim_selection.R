# Generates the bundled example data set 'sim_selection'.
# Run with the package loaded, for example via devtools::load_all().

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
