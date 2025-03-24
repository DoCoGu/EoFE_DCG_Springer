rm(list=ls())

library(readxl)

# Import Stock Returns
Ret <- read_excel('../Data/Returns.xlsx')
Factors <- read_excel('../Data/FF_Factors.xlsx')

FF <- Factors[2:276,2:4]/100
Rf <- Factors[, 5] / 100
Mkt <- Factors[2:276, 2] / 100



ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]
Mkt <- as.matrix(Mkt[2:276,])
FF <- as.matrix(FF)
# Testing the CAPM (One-Factor Model)

# Stocks
mS <- ncol(ExRet)
alpha <- matrix(NA, nrow = 2, ncol = mS)
beta1 <- matrix(NA, nrow = 2, ncol = mS)
beta2 <- matrix(NA, nrow = 2, ncol = mS)
beta3 <- matrix(NA, nrow = 2, ncol = mS)
R2_S <- numeric(mS)

for (i in 1:mS) {
  mdl <- lm(ExRet[, i] ~ FF)
  alpha[1, i] <- coef(mdl)[1]
  alpha[2, i] <- summary(mdl)[["coefficients"]][, "t value"][1]
  beta1[1, i] <- coef(mdl)[2]
  beta1[2, i] <- summary(mdl)[["coefficients"]][, "t value"][2]
  beta2[1, i] <- coef(mdl)[2]
  beta2[2, i] <- summary(mdl)[["coefficients"]][, "t value"][3]
  beta3[1, i] <- coef(mdl)[2]
  beta3[2, i] <- summary(mdl)[["coefficients"]][, "t value"][4]
  
  
  R2_S[i] <- summary(mdl)$adj.r.squared
}


# Display Results
Stocks <- data.frame(
  alpha = c(alpha[1, ]),
  `(alpha t-stat)` = c(alpha[2, ]),
  beta1 = c(beta1[1, ]),
  `(beta t-stat)` = c(beta1[2, ]),
  beta2 = c(beta1[1, ]),
  `(beta t-stat)` = c(beta2[2, ]),
  beta3 = c(beta1[1, ]),
  `(beta t-stat)` = c(beta3[2, ]),
  
    `Adj. R2` = R2_S
)

Stocks <- t(Stocks)

rownames(Stocks) <- c("alpha", "(alpha t-stat)", "mkt", "(mkt t-stat)",
                             "smb", "(smb t-stat)", "hml", "(hml t-stat)", "Adj. R2")
colnames(Stocks) <- colnames(ExRet)

# Write the results to a CSV file
write.csv(Stocks, file = "CAPM_3F.csv", row.names = TRUE)
