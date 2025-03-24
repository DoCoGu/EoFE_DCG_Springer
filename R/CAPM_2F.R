rm(list = ls())


library(readxl)

# Import Stock Returns
Ret <- read_excel('../Data/Returns.xlsx')
Factors <- read_excel('../Data/FF_Factors.xlsx')
Rf <- as.matrix(Factors[, 5] / 100)
Mkt <- as.matrix(Factors[2:276, 2] / 100)

ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]

# Import Uncertainty data
Uncertainty <- read_excel('../Data/Uncertainty.xlsx')

# Calculate VIX and EPU
VIX <- as.matrix((Uncertainty[-1, 2]) / (Uncertainty[-nrow(Uncertainty), 2]) - 1)
EPU <- as.matrix((Uncertainty[-1, 3]) / (Uncertainty[-nrow(Uncertainty), 3]) - 1)
# Testing the CAPM (One-Factor Model)


# VIX
mS <- ncol(ExRet)
alpha_VIX <- matrix(NA, nrow = 2, ncol = mS)
beta1_VIX <- matrix(NA, nrow = 2, ncol = mS)
beta2_VIX <- matrix(NA, nrow = 2, ncol = mS)
R2_VIX <- numeric(mS)

for (i in 1:mS) {
  mdl <- lm(ExRet[, i] ~ Mkt + VIX)
  alpha_VIX[1, i] <- coef(mdl)[1]
  alpha_VIX[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1]
  beta1_VIX[1, i] <- coef(mdl)[2]
  beta1_VIX[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2]
  beta2_VIX[1, i] <- coef(mdl)[3]
  beta2_VIX[2, i] <- summary(mdl)[["coefficients"]][, "t value"][3]
  R2_VIX[i] <- summary(mdl)$adj.r.squared
}

# EPU
alpha_EPU <- matrix(NA, nrow = 2, ncol = mS)
beta1_EPU <- matrix(NA, nrow = 2, ncol = mS)
beta2_EPU <- matrix(NA, nrow = 2, ncol = mS)
R2_EPU <- numeric(mS)

for (i in 1:mS) {
  mdl <- lm(ExRet[, i] ~ Mkt + EPU)
  alpha_EPU[1, i] <- coef(mdl)[1]
  alpha_EPU[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1]
  beta1_EPU[1, i] <- coef(mdl)[2]
  beta1_EPU[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2]
  beta2_EPU[1, i] <- coef(mdl)[3]
  beta2_EPU[2, i] <- summary(mdl)[["coefficients"]][, "t value"][3]
  R2_EPU[i] <- summary(mdl)$adj.r.squared
}

# Display Results
Stock_VIX <- data.frame(
  alpha = c(alpha_VIX[1, ]),
  `(alpha t-stat)` = c(alpha_VIX[2, ]),
  beta1 = c(beta1_VIX[1, ]),
  `(beta t-stat)` = c(beta1_VIX[2, ]),
  beta2 = c(beta2_VIX[1, ]),
  `(beta t-stat)` = c(beta2_VIX[2, ]),
  `Adj. R2` = R2_VIX
)


Stock_VIX <- t(Stock_VIX)

rownames(Stock_VIX) <- c("alpha", "(alpha t-stat)", "Mkt", "(Mkt t-stat)", "VIX", "(VIX t-stat)", "Adj. R2")
colnames(Stock_VIX) <- colnames(ExRet)

# Display Results
Stock_EPU <- data.frame(
  alpha = c(alpha_EPU[1, ]),
  `(alpha t-stat)` = c(alpha_EPU[2, ]),
  beta1 = c(beta1_EPU[1, ]),
  `(beta t-stat)` = c(beta1_EPU[2, ]),
  beta2 = c(beta2_EPU[1, ]),
  `(beta t-stat)` = c(beta2_EPU[2, ]),
  `Adj. R2` = R2_EPU
)


Stock_EPU <- t(Stock_EPU)

rownames(Stock_EPU) <- c("alpha", "(alpha t-stat)", "Mkt", "(Mkt t-stat)", "EPU", "(EPU t-stat)", "Adj. R2")
colnames(Stock_EPU) <- colnames(ExRet)

# Write the results to a CSV file
write.csv(Stock_VIX, file = "VIX_Stock.csv", row.names = TRUE)
write.csv(Stock_EPU, file = "EPU_Stock.csv", row.names = TRUE)
