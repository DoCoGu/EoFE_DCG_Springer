library(readxl)
library(writexl)
library(ggplot2)

# Load data
Ret <- read_excel('../Data/Returns.xlsx')
Factors <- read_excel('../Data/FF_Factors.xlsx')

# Extract required data
Rf <- as.matrix(Factors[-1, 5] / 1200)
R <- as.matrix(Ret[, -1])
N <- ncol(R)

# Calculate variables
z <- colMeans(R)
sig <- apply(R, 2, sd)
V <- cov(R)
V1 <- solve(V)

H <- t(z - mean(Rf)) %*% V1 %*% (z - mean(Rf))
sqrtH <- sqrt(H)

mu_p <- seq(mean(Rf), 0.009, length.out = 100)

A <- t(z) %*% V1 %*% z
B <- t(z) %*% V1 %*% rep(1, N)
C <- sum(rep(1, N) %*% V1 %*% rep(1, N))
D <- A * C - B^2

# Variance and standard deviation
sig2_p <- 1 / H * (mu_p - mean(Rf))^2
sig_p <- sqrt(sig2_p)

mu_t <- mean(Rf) + H / (rep(1, N) %*% V1 %*% (z - rep(mean(Rf), N)))
sig_t <- sqrtH / (rep(1, N) %*% V1 %*% (z - rep(mean(Rf), N)))

A <- t(z) %*% V1 %*% z
B <- t(z) %*% V1 %*% rep(1, N)
C <- sum(rep(1, N) %*% V1 %*% rep(1, N))
D <- A * C - B^2
mu_p <- seq(0.001, 0.008, length.out = 100)

# Variance and standard deviation
sig2_p <- (1 / D) * (C * mu_p^2 - 2 * B * mu_p + A)
sig_pp <- sqrt(sig2_p)
sig_pp[mu_p < 0.005] <- NaN
x <- cbind(sig_pp, mu_p)

# Plot

p <- ggplot() +
  geom_line(aes(x = sig_p, y = mean(Rf) + sqrtH * sig_p), color = 'blue', size = 1) +
  geom_path(aes(x = sig_pp, y = mu_p), color = 'black', size = 1) +
  geom_point(aes(x = sig_t, y = mu_t), fill = 'black', size = 3) +
  labs(
    title = 'Capital market line',
    x = 'Portfolio Risk',
    y = 'Portfolio Expected Return'
  ) +
  theme_minimal() +
  theme(legend.position = 'northwest') +
  xlim(0, 0.05) +
  geom_hline(yintercept = mu_t, linetype = 'dashed', color = 'black', size = 1) +
  geom_line(
    aes(x = sig_p, y = mean(Rf) - sqrtH * sig_p),
    color = 'blue',
    size = 1,
    linetype = 'solid'
  ) +
  geom_vline(xintercept = sig_t, linetype = 'dashed', color = 'black', size = 1)
  
  
ggsave('MV_rf.png', device = 'png')

# Optimal weights
mup <- seq(0.001, 0.05, length.out = 10)
wp <- V1 %*% (z - rep(mean(Rf), N)) %*% as.numeric((mup - mean(Rf)) / H)

# Tangency portfolio weights
wT <- (V1 %*% (z - rep(mean(Rf), N))) / sum(V1 %*% (z - rep(mean(Rf), N)))

# Weight on the risk-free asset
rf <- 1 - colSums(cbind(wT, wp))

# Create Weights table
library(tibble)


Weights <- data.frame(rbind(cbind(wT,wp),rf))

colnames(Weights) <- c("GMVP", paste0(round(mup * 100, 2), "%"))
rownames(Weights) <- c(colnames(Ret)[-1],"rf")


# Save Weights to Excel
write_xlsx(Weights, "Weights_MV_rf.xlsx")
