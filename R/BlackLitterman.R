# The Black-Litterman Model - Application
library(readxl)
library(dplyr)

source("getuncef.R")


# Import and prepare data
ret <- read_excel("../Data/Returns.xlsx")  # 5 assets
ff <- read_excel("../Data/FF_Factors.xlsx")  # Fama-French factors

ret$`Mkt-RF` <- ff$`Mkt-RF`[2:length(ff$`Mkt-RF`)] / 100
ret$rf <- ff$RF[2:length(ff$RF)] / 100

# Excess returns
stock_returns <- as.matrix(ret[, 2:(ncol(ret) - 2)])  # matrix of size 275 x 5
rf <- as.matrix(ret[, ncol(ret)])  # vector of size 275 x 1

exret <- stock_returns - rep(rf, times = ncol(stock_returns))
exmkt <- as.matrix(ret[, ncol(ret) - 1])  # excess returns on the market

# Compute the equilibrium returns
T <- nrow(exret)
n <- ncol(exret)
beta <- matrix(NA, 2, n)

for (i in 1:n) {
  mdl <- lm(exret[, i] ~ exmkt)
  beta[, i] <- coef(mdl)[2]
}

muExret <- colMeans(exret)
Sigma <- cov(exret)  # Covariance matrix of the excess returns
Pi <- beta[2, ] * muExret  # Market equilibrium returns
tau <- 1 / T  # Uncertainty factor (as defined by Litterman and He, 1999)

# Define the Q and P matrices
Q <- c(0.04 / 12,  # Divided by 12 to express it monthly
       0.02 / 12, 
       0.10 / 12)

P <- matrix(c(1, 0, -1, 0, 0, 
              0, -1, 0, 1, 0, 
              0, 0, 0, 0, 1), ncol = 5, byrow = TRUE)

# Omega, the covariance matrix of the views
Omega <- P %*% (tau * Sigma) %*% t(P)

# Blend the equilibrium returns Pi with the views
SigmaBL <- solve(solve(tau * Sigma) + t(P) %*% solve(Omega) %*% P)
muBL <- SigmaBL %*% (solve(tau * Sigma) %*% Pi + t(P) %*% solve(Omega) %*% Q)


muR <- seq(0, 0.015, by = 0.00001)


# Table of muExret, Pi, muBL
muTable <- data.frame(rbind(muExret, Pi, t(muBL)))
rownames(muTable) <- c("exret", "Pi", "muBL")
write.csv(muTable, "muTable.csv")

# Efficient frontiers
muR <- seq(0, 0.015, by = 0.00001)

# Function to calculate the efficient frontier (assuming you have it defined elsewhere)
result_MV <- getuncef(muExret, Sigma, muR)
result_BL  <- getuncef(muBL, Sigma + SigmaBL, muR)

ef_table <- data.frame(
  r = muR,
  sigma = sqrt(result_MV$OptSigma),
  sigmaBL = sqrt(result_BL$OptSigma)
)

# Plot the efficient frontiers
library(ggplot2)

EF_plot <- ggplot(ef_table, aes(x = sigma, y = r)) +
  geom_path(color = 'blue', size = 1) +
  geom_path(aes(x = sigmaBL), color = 'red', size = 1) +
  ggtitle('Efficient Frontier') +
  xlim(0, 0.1) +
  labs(x = 'Portfolio Risk', y = 'Portfolio Expected Return') +
  theme_minimal() +
  theme(plot.title = element_text(size = 16),
        axis.title = element_text(size = 16)) +
  scale_color_manual(values = c("blue", "red")) +
  theme(legend.position = "northwest")



# Optimal weights
w <- result_MV$w
wBL <- result_BL$w

weights <- data.frame(matrix(ncol = n, nrow = 2))
colnames(weights) <- colnames(ret)[1:5]
rownames(weights) <- c("M-V", "B-L")

weights[1, ] <- w[, which(ef_table$r == 0.0060)]  # M-V weights at target return
weights[2, ] <- wBL[, which(ef_table$r == 0.0060)]  # B-L weights at target return

write.csv(weights, "Weights_BL.xlsx")
