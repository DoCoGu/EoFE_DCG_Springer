
# Load required libraries
library(readxl)
library(lubridate)
library(ggplot2)

# Import Stock Returns
Ret <- read_excel('../Data/Returns.xlsx')
Factors <- read_excel('../Data/FF_Factors.xlsx')
Rf <- Factors[, 5] / 100
Mkt <- Factors[2:276, 2] / 100


ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]
Mkt <- as.matrix(Mkt)
# Rolling-Window Size
rw <- 40

# Stocks
nS <- nrow(ExRet)
mS <- ncol(ExRet)
alpha <- matrix(NA, nS - rw + 1, mS)
beta <- matrix(NA, nS - rw + 1, mS)

for (t in 1:(nS - rw + 1)) {
  for (i in 1:mS) {
    mdl <- lm(ExRet[t:(t + rw - 1), i] ~ Mkt[t:(t + rw - 1)])
    alpha[t, i] <- coef(mdl)[1]
    beta[t, i] <- coef(mdl)[2]
  }
}

# %% Plots

# Generate dates for plotting
Dates_for_plot <- seq(as.Date(Ret[rw, 1]$Date), as.Date(Ret[nS, 1]$Date), length.out = 20)
Dates_for_plot <- format(Dates_for_plot, "%b-%Y")

alpha_plot <- ggplot() +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 1]), color = "red", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 2]), color = "black", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 3]), color = "blue", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 4]), color = "green", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 5]), color = "magenta", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(1, nS - rw + 1, length.out = 20), labels = Dates_for_plot) +
  labs(title = "Alpha", x = "", y = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 

ggsave("alpha_S.eps", plot = alpha_plot, device = "eps")

# Beta plot

beta_plot <- ggplot() +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 1]), color = "red", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 2]), color = "black", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 3]), color = "blue", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 4]), color = "green", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 5]), color = "magenta", size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = seq(1, nS - rw + 1, length.out = 20), labels = Dates_for_plot) +
  labs(title = "Beta", x = "", y = "") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) 
ggsave("beta_S.eps", plot = beta_plot, device = "eps")

