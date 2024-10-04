library(spam)

#### generate gaussian processes in a [0,2] by [0,2] windows with resolution 100 by 100
#### time-comsuming (half an hour), need large ram (16 G)

set.seed(1)
nrows <- 100
ncols <- 100
x <- seq(0, 2, length.out = ncols)
y <- seq(0, 2, length.out = nrows)
grid <- expand.grid(x = x, y = y)
range_param <- 0.05
sill_param <- 1
covariance_function <- function(d) sill_param * exp(-d / range_param)
distances <- as.matrix(dist(grid))
cov_matrix <- covariance_function(distances)
simulated_field_1_large <- rmvnorm(n = 1, mean = rep(0, nrow(grid)), sigma = cov_matrix)
simulated_field_2_large <- rmvnorm(n = 1, mean = rep(0, nrow(grid)), sigma = cov_matrix)

write.csv(simulated_field_1_large, file = "simulated_field_1_large.csv", row.names = FALSE)
write.csv(simulated_field_2_large, file = "simulated_field_2_large.csv", row.names = FALSE)