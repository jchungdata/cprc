# Reproducible generation script for inst/extdata/synth_spm.csv
# Replicates Simulation Setting 2 of Wang et al. (2018), Human Brain Mapping, 39(12), 4811-4825.
# 121 voxels (11 x 11 slice) x 20 subjects:
#   - Voxels 1..43: active region; Voxels 44..121: null background
#   - Subjects 1..17: consistent moderate activation
#   - Subjects 18..20: outlying subjects with inflated variance across voxels 9..43

set.seed(2018)
M <- 121L
N <- 20L
X <- matrix(rnorm(M * N, mean = 0, sd = 1), nrow = M, ncol = N)

# Moderate consensus activation in voxels 1..43 for subjects 1..17
X[1:43, 1:17] <- matrix(rnorm(43 * 17, mean = 2.0, sd = 0.8), nrow = 43, ncol = 17)

# Outlying subjects 18..20 with inflated variance across voxels 9..43
X[1:8, 18:20]  <- matrix(rnorm(8 * 3, mean = 2.0, sd = 1.0), nrow = 8, ncol = 3)
X[9:43, 18:20] <- matrix(rnorm(35 * 3, mean = 2.0, sd = 8.0), nrow = 35, ncol = 3)

colnames(X) <- paste0("Sub", seq_len(N))
# write.csv(X, file = "inst/extdata/synth_spm.csv", row.names = FALSE)
