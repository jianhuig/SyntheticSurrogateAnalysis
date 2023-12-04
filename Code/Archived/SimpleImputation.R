source("SynSurrogateSim.R")

library(dplyr)
library(doParallel)
library(ggplot2)

n.sim <- 10^4 # number of replicates
n.sample <- 10^3 # total sample size

# calculate effect size
var.y <- 1 / (1 - 0.005 - 0.2) # 1/var.y = proportion of variance explained by res
beta.g <- sqrt(0.005 * var.y / (2 * 0.25 * 0.75)) # (var(G) * beta^2)/var.y = 0.005
beta.x <- rep(sqrt(0.2 * var.y / 2.5), 2) # beta.1^2 + beta.2^2 + 0.5*beta.1*beta.2 = 0.2 * var.y

# generate covariates
X <- mvnfast::rmvn(
  n = n.sample, mu = c(0, 0),
  sigma = matrix(c(1, 0.5, 0.5, 1), ncol = 2)
)

#-------------- Known Model ------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
result <- foreach(
  i = 1:n.sim, .combine = rbind, .packages = c("doParallel"),
  .errorhandling = "pass"
) %dopar% {
  # genotype
  g <- rbinom(n.sample, size = 2, prob = 0.25)
  # y complete
  y <- beta.g * g + X %*% beta.x + rnorm(n.sample)
  # oracle procedure
  assoc.oracle <- lm(y ~ X + g - 1)
  beta.oracle <- summary(assoc.oracle)$coefficients["g", 1]
  # missing_index
  missing.id.25 <- sample(1:n.sample, size = ceiling(n.sample * 0.25))
  missing.id.50 <- sample(1:n.sample, size = ceiling(n.sample * 0.50))
  missing.id.75 <- sample(1:n.sample, size = ceiling(n.sample * 0.75))
  # imputed y
  y.hat <- beta.g * g + X %*% beta.x
  y.imputed.25 <- y.imputed.50 <- y.imputed.75 <- y
  y.imputed.25[missing.id.25] <- y.hat[missing.id.25]
  y.imputed.50[missing.id.50] <- y.hat[missing.id.50]
  y.imputed.75[missing.id.75] <- y.hat[missing.id.75]
  # regression with imputed value
  assoc.imputed.25 <- lm(y.imputed.25 ~ X + g - 1)
  beta.imputed.25 <- summary(assoc.imputed.25)$coefficients["g", 1]
  assoc.imputed.50 <- lm(y.imputed.50 ~ X + g - 1)
  beta.imputed.50 <- summary(assoc.imputed.50)$coefficients["g", 1]
  assoc.imputed.75 <- lm(y.imputed.75 ~ X + g - 1)
  beta.imputed.75 <- summary(assoc.imputed.75)$coefficients["g", 1]

  out <- c(beta.oracle, beta.imputed.25, beta.imputed.50, beta.imputed.75)
  names(out) <- c("oracle", "25% missing", "50% missing", "75% missing")
  out
}
stopCluster(cl)

result %>%
  data.frame() %>%
  tidyr::gather() %>%
  ggplot(aes(x = key, y = value, color = key)) +
  geom_boxplot() +
  geom_hline(yintercept = beta.g, color = "red")

#-----Multiple Imputation, Incorrectly specified ------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
result <- foreach(
  i = 1:n.sim, .combine = rbind, .packages = c("doParallel"),
  .errorhandling = "pass"
) %dopar% {
  # genotype
  g <- rbinom(n.sample, size = 2, prob = 0.25)
  # y complete
  y <- beta.g * g + X %*% beta.x + rnorm(n.sample)
  # oracle procedure
  assoc.oracle <- lm(y ~ X + g)
  beta.oracle <- summary(assoc.oracle)$coefficients["g", 1]
  # missing_index
  missing.id.25 <- sample(1:n.sample, size = ceiling(n.sample * 0.25))
  missing.id.50 <- sample(1:n.sample, size = ceiling(n.sample * 0.50))
  missing.id.75 <- sample(1:n.sample, size = ceiling(n.sample * 0.75))
  # imputed y
  y.hat <- beta.g * g + X %*% beta.x
  y.imputed.25 <- y.imputed.50 <- y.imputed.75 <- y
  y.imputed.25[missing.id.25] <- y.hat[missing.id.25]
  y.imputed.50[missing.id.50] <- y.hat[missing.id.50]
  y.imputed.75[missing.id.75] <- y.hat[missing.id.75]
  # regression with imputed value
  assoc.imputed.25 <- lm(y.imputed.25 ~ X + g)
  beta.imputed.25 <- summary(assoc.imputed.25)$coefficients["g", 1]
  assoc.imputed.50 <- lm(y.imputed.50 ~ X + g)
  beta.imputed.50 <- summary(assoc.imputed.50)$coefficients["g", 1]
  assoc.imputed.75 <- lm(y.imputed.75 ~ X + g)
  beta.imputed.75 <- summary(assoc.imputed.75)$coefficients["g", 1]

  out <- c(beta.oracle, beta.imputed.25, beta.imputed.50, beta.imputed.75)
  names(out) <- c("oracle", "25% missing", "50% missing", "75% missing")
  out
}
stopCluster(cl)

result %>%
  data.frame() %>%
  tidyr::gather() %>%
  ggplot(aes(x = key, y = value, color = key)) +
  geom_boxplot() +
  geom_hline(yintercept = beta.g, color = "red")

#-------------- Changing the target outcome ------
cl <- makeCluster(detectCores())
registerDoParallel(cl)
result <- foreach(
  i = 1:n.sim, .combine = rbind, .packages = c("doParallel", "dplyr"),
  .errorhandling = "pass"
) %dopar% {
  # genotype
  g <- rbinom(n.sample, size = 2, prob = 0.25)
  # y complete
  y <- beta.g * g + X %*% beta.x + rnorm(n.sample)
  # oracle procedure
  assoc.oracle <- lm(y ~ X + g)
  beta.oracle <- summary(assoc.oracle)$coefficients["g", 1]
  # missing_index
  missing.id.25 <- sample(1:n.sample, size = ceiling(n.sample * 0.25))
  missing.id.50 <- sample(1:n.sample, size = ceiling(n.sample * 0.50))
  missing.id.75 <- sample(1:n.sample, size = ceiling(n.sample * 0.75))
  # missing y
  y.imputed.25 <- y.imputed.50 <- y.imputed.75 <- y
  y.imputed.25[missing.id.25] <- NA
  y.imputed.50[missing.id.50] <- NA
  y.imputed.75[missing.id.75] <- NA
  # regression with imputed value
  assoc.imputed.25 <- lm(y.imputed.25 ~ X + g)
  y.imputed.25 <- predict(assoc.imputed.25, cbind(X, g) %>% data.frame())
  assoc.imputed.25 <- lm(y.imputed.25 ~ X + g)
  beta.imputed.25 <- summary(assoc.imputed.25)$coefficients["g", 1]
  assoc.imputed.50 <- lm(y.imputed.50 ~ X + g)
  y.imputed.50 <- predict(assoc.imputed.50, cbind(X, g) %>% data.frame())
  assoc.imputed.50 <- lm(y.imputed.50 ~ X + g)
  beta.imputed.50 <- summary(assoc.imputed.50)$coefficients["g", 1]
  assoc.imputed.75 <- lm(y.imputed.75 ~ X + g)
  y.imputed.75 <- predict(assoc.imputed.75, cbind(X, g) %>% data.frame())
  assoc.imputed.75 <- lm(y.imputed.75 ~ X + g)
  beta.imputed.75 <- summary(assoc.imputed.75)$coefficients["g", 1]

  out <- c(beta.oracle, beta.imputed.25, beta.imputed.50, beta.imputed.75)
  names(out) <- c("oracle", "25% missing", "50% missing", "75% missing")
  out
}
stopCluster(cl)

result %>%
  data.frame() %>%
  tidyr::gather() %>%
  ggplot(aes(x = key, y = value, color = key)) +
  geom_boxplot() +
  geom_hline(yintercept = beta.g, color = "red")
