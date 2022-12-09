library(dplyr)
library(tibble)

args = commandArgs(trailingOnly=TRUE)
print(args)
expr_fname = args[1]
expr_noise_out = args[2]
gauss_sd = args[3]
n_rep = args[4]

# Extract dimensions of expression matrix
expr <- readRDS(expr_fname)
n <- dim(expr)[[1]]
m <- dim(expr)[[2]]

# Generate and add gaussian noise
set.seed(42)
gauss_noise <- matrix(rnorm(n*m, mean = 0, sd = as.numeric(gauss_sd)), n, m)
expr_noise <- expr + gauss_noise

# Shift expression so it's non-negative
if (min(expr_noise) < 0) {
	expr_noise <- expr_noise + abs(min(expr_noise)) + 1e-16
}

# Save 
saveRDS(expr_noise, expr_noise_out)