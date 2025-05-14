# Load required libraries
library(readxl)
library(dplyr)
library(lubridate)
library(stats)
library(ggplot2)
# Settings
America <- 1  # Use only aviation disasters in America
cutoff <- 150  # Cutoff for the number of casualties
start_CAR <- -1
end_CAR <- 5

start_CAPM <- -250
end_CAPM <- -50

# Import Data
Portfolios <- read_excel("../Data/IndustryPortfolios.xlsx")
Factors <- read_excel("../Data/FF_FactorsDaily.xlsx")
Events <- read_excel("../Data/Events.xlsx")

if (America == 0) {
  Events <- Events[which(Events$`Total Deaths` > cutoff),4]
} else {
  Events <- Events[which(Events$`Total Deaths` > cutoff & Events$Zone==1),4]
}

Events <- Events$`Date start`

MktRet <- Factors$`Mkt-RF` / 100
Rf <- Factors$RF / 100

Ret <- Portfolios$Trans / 100
ExRet <- as.matrix(Ret - Rf)

# Convert dates to Date format
Events <- dmy(Events)
DatesReturn <- ymd(Factors$Date)

n <- nrow(ExRet)
nS <- ncol(ExRet)

if (is.null(nS)){
  nS = 1
}

# Get event's date in returns dates
Dates <- match(Events, DatesReturn)
for (i in 1:length(Events)) {
  while (is.na(Dates[i])) {  # If Dates[i] is not a trading day, try the next day
    Events[i] <- Events[i] + days(1)
    Dates[i] <- match(Events[i], DatesReturn)
  }
}

Dates <- Dates[Dates > abs(start_CAPM)]
Dates <- unique(Dates)
nE <- length(Dates)

# Event Study
# Preallocate vectors for CAPM coefficients
alpha <- matrix(NA, nrow=nE, ncol=nS)
beta <- matrix(NA, nrow=nE, ncol=nS)

# Compute the CAPM model for each stock and around each event
for (i in 1:nE) {
  for (k in 1:nS) {
    onefctmdl <- lm(ExRet[(Dates[i] + start_CAPM) : (Dates[i] + end_CAPM)] ~ MktRet[(Dates[i] + start_CAPM):(Dates[i] + end_CAPM)])
    alpha[i, k] <- coef(onefctmdl)[1]
    beta[i, k] <- coef(onefctmdl)[2]
  }
}

# Preallocate a vector of CAPM predicted returns
PredRet <- array(NA, c(abs(start_CAR) + abs(end_CAR) + 1, nE, nS))

# Get predicted returns for each stock and around each event, for the CARs period
for (t in 1:(abs(start_CAR) + abs(end_CAR) + 1)) {
  for (i in 1:nE) {
    for (k in 1:nS) {
      PredRet[t, i, k] <- alpha[i, k] + beta[i, k] * (MktRet[Dates[i] + start_CAR - 1 + t])
    }
  }
}

# Preallocate a vector of observed returns
ObsRet_agg <- array(NA, c(abs(start_CAR) + abs(end_CAR) + 1, nE, nS))

# Get observed returns for each stock and around each event, for the CARs period
for (i in 1:nE) {
  for (t in 1:(abs(start_CAR) + abs(end_CAR) + 1)) {
    for (k in 1:nS) {
      ObsRet_agg[t, i, k] <- ExRet[Dates[i] + start_CAR - 1 + t, k]
    }
  }
}

# Get abnormal returns (Observed ret. - Predicted ret.)
AbnRet <- ObsRet_agg - PredRet

# Get cumulative abnormal returns
CAR <- apply(AbnRet, 2, function(x) cumsum(x))
CAAR <- apply(CAR, 1, mean, na.rm = TRUE) * 100

# Plot
date <- seq(start_CAR, end_CAR, by = 1)
zero <- rep(0, length(date))

plot_data <- data.frame(date = date, CAAR = CAAR)

# Create the ggplot
CAR_plot <- ggplot(plot_data, aes(x = date, y = CAAR)) +
  geom_line(color = "blue") +
  labs(x = "Days Relative to Events", y = "CAR") +
  xlim(start_CAR, end_CAR) +
  ylim(min(CAAR) - 0.001, max(CAAR) + 0.001) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(start_CAR, end_CAR, 1), limits = c(start_CAR, end_CAR)) +
  scale_color_manual(values = c("CAR" = "blue")) +
  guides(color = guide_legend(title = "Legend", override.aes = list(linetype = 1, size = 2))) +
  labs(color = "Legend")+
  theme_minimal() +
  theme(legend.position = "top")
ggsave("CARs.eps", plot = CAR_plot, device = "eps")
