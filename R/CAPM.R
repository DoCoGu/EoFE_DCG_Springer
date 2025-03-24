library(readxl)


# Import Stock Returns
Ret <- read_excel('../Data/Returns.xlsx')
Factors <- read_excel('../Data/FF_Factors.xlsx')
Rf <- Factors[, 5] / 100
Mkt <- Factors[2:276, 2] / 100



ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]
Mkt <- as.matrix(Mkt)
# Testing the CAPM (One-Factor Model)

# Stocks
mS <- ncol(ExRet)
alpha_S <- matrix(NA, nrow = 2, ncol = mS)
beta_S <- matrix(NA, nrow = 2, ncol = mS)
R2_S <- numeric(mS)

for (i in 1:mS) {
  mdl <- lm(ExRet[, i] ~ Mkt)
  alpha_S[1, i] <- coef(mdl)[1]
  alpha_S[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1]
  beta_S[1, i] <- coef(mdl)[2]
  beta_S[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2]
  R2_S[i] <- summary(mdl)$adj.r.squared
}


# Display Results
Results <- data.frame(
  alpha = c(alpha_S[1, ]),
  `(alpha t-stat)` = c(alpha_S[2, ]),
  beta = c(beta_S[1, ]),
  `(beta t-stat)` = c(beta_S[2, ]),
  `Adj. R2` = R2_S
)

Results <- t(Results)

rownames(Results) <- c("alpha", "(alpha t-stat)", "beta", "(beta t-stat)", "Adj. R2")
colnames(Results) <- colnames(ExRet)

# Write the results to a CSV file
write.csv(Results, file = "CAPM_Stock.csv", row.names = TRUE)
