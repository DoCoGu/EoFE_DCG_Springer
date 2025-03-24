Sys.setlocale("LC_TIME", "en_US")

# Load required libraries
library(readxl)
library(dplyr)
library(lubridate)
library(sandwich)
library(lmtest)
# Setting
maxlag <- 3
cutoff <- 150
fD <- 3
cons <- 1
America <- 1

# Import Data
Ret <- read_excel("../Data/IndustryPortfolios.xlsx", col_names = TRUE)
DatesReturn <- ymd(Ret$...1)

Events <- read_excel("../Data/Events.xlsx", col_names = TRUE)

if (America == 0) {
  Events <- Events[which(Events$`Total Deaths` > cutoff),4]
} else {
  Events <- Events[which(Events$`Total Deaths` > cutoff & Events$Zone==1),4]
}

ExRet <- Ret$Trans

# Convert dates to Date format
Events <- dmy(Events$`Date start`)

# Preallocate vectors of weekday dummies
dow <- matrix(NA, nrow = length(DatesReturn), ncol = 4)

# Create dummy for weekdays
days <- c("Monday", "Tuesday", "Wednesday", "Thursday")
for (d in 1:4) {
  for (i in 1:length(DatesReturn)) {
    if (weekdays(DatesReturn[i]) == days[d]) {
      dow[i, d] <- 1
    } else {
      dow[i, d] <- 0
    }
  }
}

# Create events dummy
Match <- DatesReturn %in% Events

EDummy <- matrix(0, nrow = length(DatesReturn)+3, ncol = 3)

for (i in 1:length(DatesReturn)) {
  if (Match[i]) {
    EDummy[i + 1, 1] <- 1
    EDummy[i + 2, 2] <- 1
    EDummy[i + 3, 3] <- 1
  } else {
    EDummy[i + 1, 1] <- 0
    EDummy[i + 2, 2] <- 0
    EDummy[i + 3, 3] <- 0
  }
}

EDummy <- EDummy[-c(1:3), ]

TaxDay <- rep(0, length(Ret))
dday <- day(DatesReturn)
dmonth <- month(DatesReturn)

TaxDay[dmonth == 1 & dday >= 1 & dday <= 5] <- 1
TaxDay[is.na(TaxDay)] = 0
Res <- matrix(NA, nrow = 2, ncol = 1)

RetLag <-cbind(lag(ExRet, n = 1), lag(ExRet, n = 2), lag(ExRet, n = 3)) 

# Regression
X <- cbind(RetLag, dow, TaxDay, EDummy[,1:fD])

if (cons == 1) {
  mdl <- lm(ExRet ~ X)
  coeff <- round(coef(mdl), 5)
  se <- sqrt(diag(vcov(mdl)))
} else {
  mdl <- lm(ExRet ~ X - 1)
  coeff <- round(coef(mdl), 5)
  se <- sqrt(diag(vcov(mdl)))
}

Sig <- coeff / se
res <- matrix(rbind(coeff, Sig), nrow = 2)
varnames <- c("Mon", "Tue", "Wed", "Thu", "TaxDays")
rString <- paste("R_t-", 1:maxlag, sep = "")
eString <- paste("Event +", 1:fD)
var <- c("Constant", rString, varnames, eString)

tab <- data.frame(matrix(NA, ncol = length(var), nrow = 2))
colnames(tab) <- var
rownames(tab) <- c("Coefficient", "t-Statistic")
tab[] <- res

if (America == 0) {
  write.table(tab, "RegressionAnalysis_All.csv", sep = ",", col.names = TRUE, row.names = T)
} else {
  write.table(tab, "RegressionAnalysis_America.csv", sep = ",", col.names = TRUE, row.names = T)
}
