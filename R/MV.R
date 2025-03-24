library(readxl)
library(writexl)
library(ggplot2)

# Read data from Excel file
Ret <- read_xlsx("../Data/Returns.xlsx", col_names = TRUE)

# Extract returns data
R <- as.matrix(Ret[, -1])
N <- ncol(R)

z <- colMeans(R)
sig <- apply(R, 2, sd)
V <- cov(R)
V1 <- solve(V)

A <- t(z) %*% V1 %*% z
B <- t(z) %*% V1 %*% rep(1, N)
C <- rep(1, N) %*% V1 %*% rep(1, N)
D <- A * C - B^2

mu_p <- seq(0, 0.015, by = 0.0001)

sig2_p <- (1 / D) * (C * mu_p^2 - 2 * B * mu_p + A)
sig_p <- sqrt(sig2_p)

# Plot

eff_frontier <- data.frame(sig_p, mu_p)
colnames(eff_frontier) <- c("Portfolio Risk", "Portfolio Expected Return")

eff_frontier_point<-data.frame(`Portfolio Risk` = 1 / sqrt(C), `Portfolio Expected Return` = B / C)
colnames(eff_frontier_point) <- c("Portfolio Risk", "Portfolio Expected Return")

p <- ggplot(eff_frontier, aes(x = `Portfolio Risk`, y = `Portfolio Expected Return`)) +
  geom_path(color = "black", size = 1) +
  geom_hline(yintercept = B / C, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = 1 / sqrt(C), linetype = "dashed", color = "black", size = 1) +
  geom_point(data = data.frame(sig, z), aes(x = sig, y = z), color = c("orange","red","blue","green","pink"), size = 3, shape = 21, fill =c("orange","red","blue","green","pink")) +
  geom_point(data = eff_frontier_point, aes(x = `Portfolio Risk`, y = `Portfolio Expected Return`), color = "black", size = 3, shape = 21, fill = "black") +
  xlim(0, 0.1) +
  labs(
    title = "Efficient Frontier",
    x = "Portfolio Risk",
    y = "Portfolio Expected Return"
  ) +
  theme_minimal()

ggsave("MV.png", plot = p, width = 10, height = 6)

# Optimal weights
g <- as.numeric((1 / D)) * (as.numeric(A) * (V1 %*% rep(1, N)) - as.numeric(B) * (V1) %*% z)
h <- as.numeric(1 / D) * (as.numeric(C) * (V1 %*% z) - as.numeric(B) * (V1 %*% rep(1, N)))

mup <- seq(0.001, 0.05, length.out = 10)
wp <- rep(g,N) + h %*% mup

# Global minimum variance portfolio weights
w_mvp <- g + h * as.numeric(B / C)

Weights <- data.frame(cbind(w_mvp,wp)
)

colnames(Weights) <- c("GMVP", paste0(round(mup * 100, 2), "%"))
rownames(Weights) <- colnames(Ret)[-1]

# Save weights to Excel
write_xlsx(Weights, "Weights_MV.xlsx")
