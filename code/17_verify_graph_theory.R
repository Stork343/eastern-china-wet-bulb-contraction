###############################################################################
# Numerical checks of the graph-dispersion identities used in the manuscript.
###############################################################################

library(data.table)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript")
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(script_path)
source(file.path(code_dir, "esh_utils.R"))
source(file.path(code_dir, "graph_esh_utils.R"))

set.seed(20260871L)
n <- 25L
sites <- data.table(site_id = seq_len(n), lon = runif(n, 105, 125),
                    lat = runif(n, 20, 42))
op <- make_graph_operators(sites, 0.5)[[1]]
y <- rnorm(n)
z <- y - mean(y)

# Laplacian quadratic form equals the weighted pairwise-difference sum.
pairwise_sum <- 0
for (i in seq_len(n - 1L)) for (j in (i + 1L):n) {
  pairwise_sum <- pairwise_sum + op$W[i, j] * (y[i] - y[j])^2
}
matrix_sum <- as.numeric(crossprod(y, op$L %*% y))
stopifnot(all.equal(pairwise_sum, matrix_sum, tolerance = 1e-10))

# Translation invariance and degree-two scale equivariance.
D <- graph_dispersion(y, op)
stopifnot(all.equal(D, graph_dispersion(y + 17, op), tolerance = 1e-10))
stopifnot(all.equal(9 * D, graph_dispersion(3 * y, op), tolerance = 1e-10))

# Spectral Poincare bounds for centered signals.
lower <- op$lambda_2 * sum(z^2) / (2 * op$weight_sum)
upper <- op$lambda_max * sum(z^2) / (2 * op$weight_sum)
stopifnot(D >= lower - 1e-10, D <= upper + 1e-10)

# Complete equal-weight graph reduces to n/(n-1) times population variance.
W_complete <- matrix(1, n, n); diag(W_complete) <- 0
L_complete <- diag(rowSums(W_complete)) - W_complete
op_complete <- list(L = L_complete, weight_sum = n * (n - 1) / 2)
D_complete <- graph_dispersion(y, op_complete)
V_population <- mean((y - mean(y))^2)
stopifnot(all.equal(D_complete, n / (n - 1) * V_population,
                    tolerance = 1e-10))

cat("All graph-dispersion identities and spectral bounds passed.\n")

