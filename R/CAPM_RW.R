# CAPM_RW.R
# This script implements a rolling-window regression to estimate time-varying
# CAPM alpha and beta for a set of assets and plots the results in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Load required libraries
library(readxl)   # For reading Excel files
library(lubridate) # For working with dates
library(ggplot2)  # For creating plots

# --- Data Loading and Preparation ---
# Import Stock Returns from an Excel file
Ret <- read_excel('../Data/Returns.xlsx')
# Import Fama-French Factors from an Excel file
Factors <- read_excel('../Data/FF_Factors.xlsx')

# Extract Risk-Free Rate (Rf) and Market Excess Return (Mkt)
# Select the 5th column (index 5) for Rf, convert from percentage to decimal
Rf <- Factors[, 5] / 100
# Select the 2nd column (index 2) for Mkt-RF, convert from percentage to decimal
# Note: The original code uses rows 2 to 276 for Mkt, assuming a specific data range.
# Using rows 2 to end for consistency with other scripts if possible, or keep original range if intentional.
# Keeping the original range 2:276 as provided.
Mkt <- Factors[2:276, 2] / 100

# Calculate Excess Returns for assets (Asset Return - Risk-Free Rate)
# Convert asset return columns (all except the first column, which is date) to a matrix
# Subtract the Risk-Free Rate (Rf) from row 2 to 276 to match the length of ExRet
ExRet <- as.matrix(Ret[, -1]) - Rf[2:276,]
# Convert Mkt to a matrix
Mkt <- as.matrix(Mkt)

# --- Rolling Regression Setup ---
# Define the size of the rolling window
rw <- 40 # Number of observations in each window

# Get the size of the Excess Returns matrix
nS <- nrow(ExRet) # nS = total number of observations in ExRet
mS <- ncol(ExRet) # mS = number of assets

# Initialize matrices to store rolling regression results
# alpha and beta will store the coefficients for each window and asset
alpha <- matrix(NA, nS - rw + 1, mS) # Rows = number of windows, Columns = number of assets
beta <- matrix(NA, nS - rw + 1, mS)  # Rows = number of windows, Columns = number of assets

# --- Perform Rolling Window Regression ---
# Outer loop iterates through the starting point (t) of each window
# The loop runs from 1 to the total number of possible windows (nS - rw + 1)
for (t in 1:(nS - rw + 1)) {
  # Inner loop iterates through each asset (i)
  for (i in 1:mS) {
    # Fit the linear regression model for the current window and asset
    # Model: ExRet_i = alpha_i + beta_i * Mkt + epsilon_i
    # Use the lm() function for linear modeling
    # Dependent variable: ExRet for the current window (rows t to t+rw-1), for asset i
    # Independent variable: Mkt for the current window (rows t to t+rw-1)
    mdl <- lm(ExRet[t:(t + rw - 1), i] ~ Mkt[t:(t + rw - 1)])

    # Extract and store the regression results (coefficients) for the current window
    # coef(mdl)[1] is the intercept (alpha)
    # coef(mdl)[2] is the coefficient for Mkt (beta)
    alpha[t, i] <- coef(mdl)[1]
    beta[t, i] <- coef(mdl)[2]
  }
}

# --- Plots ---

# Generate dates for the x-axis of the plots
# Select dates corresponding to the end of each rolling window (from row rw to nS of Ret)
# Create a sequence of 20 equally spaced dates within this range
Dates_for_plot <- seq(as.Date(Ret[rw, 1]$Date), as.Date(Ret[nS, 1]$Date), length.out = 20)
# Format the dates as "Month-Year" strings (e.g., "Jan-2000")
Dates_for_plot <- format(Dates_for_plot, "%b-%Y")

# Alpha plot
# Create a ggplot object for the alpha plot
alpha_plot <- ggplot() +
  # Add lines for the rolling alpha of each asset
  # aes() maps data variables to plot aesthetics (x-axis position, y-axis position)
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 1]), color = "red", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 2]), color = "black", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 3]), color = "blue", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 4]), color = "green", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = alpha[, 5]), color = "magenta", size = 1) +
  # Add a horizontal dashed line at zero for reference
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  # Configure the x-axis ticks and labels
  # breaks defines the positions of the ticks, labels provides the text for the ticks
  scale_x_continuous(breaks = seq(1, nS - rw + 1, length.out = 20), labels = Dates_for_plot) +
  # Set plot title and axis labels
  labs(title = "Rolling Alpha (CAPM)", x = "Time", y = "Alpha") + # Added meaningful axis labels
  # Apply a minimal theme for cleaner appearance
  theme_minimal() +
  # Customize the appearance of x-axis text (angle and alignment for readability)
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# Save the alpha plot as an EPS file
ggsave("alpha_S.eps", plot = alpha_plot, device = "eps")

# Beta plot
# Create a ggplot object for the beta plot
beta_plot <- ggplot() +
  # Add lines for the rolling beta of each asset
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 1]), color = "red", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 2]), color = "black", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 3]), color = "blue", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 4]), color = "green", size = 1) +
  geom_line(aes(x = 1:(nS - rw + 1), y = beta[, 5]), color = "magenta", size = 1) +
  # Add a horizontal dashed line at one for reference (CAPM beta is expected to be 1 for the market)
  geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
  # Configure the x-axis ticks and labels (same as alpha plot)
  scale_x_continuous(breaks = seq(1, nS - rw + 1, length.out = 20), labels = Dates_for_plot) +
  # Set plot title and axis labels
  labs(title = "Rolling Beta (CAPM)", x = "Time", y = "Beta") + # Added meaningful axis labels
  # Apply a minimal theme for cleaner appearance
  theme_minimal() +
  # Customize the appearance of x-axis text
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# Save the beta plot as an EPS file
ggsave("beta_S.eps", plot = beta_plot, device = "eps")


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
