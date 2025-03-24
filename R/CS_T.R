# Load required libraries
library(readxl)
library(sandwich)
library(lmtest)

setwd(getSrcDirectory(function(){})[1])

# Data
Ret <- read_xlsx('../Data/PortfoliosLong.xlsx')
FFactors <- read_xlsx('../Data/FF_FactorsLong.xlsx')
Temperature <- read_xlsx('../Data/Temperature.xlsx')
Consumption <- read_xlsx('../Data/ConsumptionLong.xlsx')


Temp <- as.matrix(diff(log(as.matrix(Temperature[, 2]))))
C <- as.matrix(Consumption[-1, 2] / Consumption[-nrow(Consumption), 2] - 1)

Mkt <- as.matrix(FFactors[-1, 2] / 100)
Rf <- as.matrix(FFactors[-1, 5] / 100)

ExRet <- as.matrix(Ret[-1, -1] / 100 - Rf)

Factors <- cbind(Mkt, C, Temp)

# First Stage Regression
n1 <- nrow(ExRet)
n2 <- ncol(ExRet)
nF <- ncol(cbind(Mkt, C, Temp))

CoefAll <- matrix(NA, nrow = nF, ncol = n2)
Res <- matrix(NA, nrow = n1, ncol = n2)

for (i in 1:n2) {
  lm_result <- lm(ExRet[, i] ~ Mkt + C + Temp)
  CoefAll[, i] <- coef(lm_result)[-1]
  Res[, i] <- residuals(lm_result)
}

VarCovErr <- cov(Res)

# Second-Stage Regression
MeanRet <- colMeans(ExRet)
Betas <- t(CoefAll)

# Fit linear model for MeanRet against Betas
mdl <- lm(MeanRet ~ Betas - 1)  # -1 removes intercept
SE <- summary(mdl)$coefficients[,2]
Lambda <- summary(mdl)$coefficients[,1]
Tstat <- Lambda / SE

# Shanken correction
Sigma_f <- cov(Factors)
VarLam <- (solve(t(Betas) %*% Betas) %*% t(Betas) %*% VarCovErr %*% Betas %*% (solve(t(Betas) %*% Betas)) *
             as.numeric(1 + t(Lambda) %*% solve(Sigma_f) %*% Lambda) + Sigma_f) / n1
SE_Shanken <- sqrt(diag(VarLam))
Tstat_Shanken <- Lambda / SE_Shanken

# T cross-section regressions
LambdaFull <- matrix(NA, n1, nF)
for (j in 1:n1) {
  MeanRet <- ExRet[j,]
  mdl <- lm(MeanRet ~ Betas - 1)
  LambdaFull[j,] <- coef(mdl)
}

LambdaMean <- colMeans(LambdaFull)

# HAC-corrected standard-errors
X <- rep(1, n1)
hac_cov <- numeric(nF)
for (k in 1:nF) {
  y <- LambdaFull[, k]
  mdl <- lm(y ~ X - 1)
  hac_cov[k] <- NeweyWest(mdl, lag = 1, sandwich = T, prewhite = F, adjust = T)
}

SE_NW <- sqrt(hac_cov)
Tstat_NW <- LambdaMean / SE_NW

# Results
NamePort <- colnames(Ret)[-1]
FirstStageReg <- data.frame(Mkt = Betas[, 1], dC = Betas[, 2], T = Betas[, 3])
rownames(FirstStageReg) <- NamePort

SecondStage <- data.frame(Lambda = Lambda,
                          Tstat = Tstat,
                          Tstat_HAC = Tstat_NW,
                          Tstat_Shanken = Tstat_Shanken)
SecondStage <- t(SecondStage)
rownames(SecondStage) <- c('Lambda', 'tstat', 't-stat HAC', 't-stat Shanken')
colnames(SecondStage) <- c('Mkt', 'C', 'T')

write.csv(FirstStageReg, 'FirstStage_EPU.csv', row.names = TRUE)
write.csv(SecondStage, 'SecondStage_EPU.csv', row.names = TRUE)
