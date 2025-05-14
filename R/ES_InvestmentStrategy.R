# ES_InvestmentStrategy.R
# This script analyzes different investment strategies, including buy-and-hold
# and an event-driven strategy based on aviation disasters, and plots the
# cumulative returns in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Load required libraries
library(readxl)    # For reading Excel files
library(dplyr)     # For data manipulation (used for lag, select)
library(lubridate) # For working with dates (ymd, day, month, weekdays)
library(ggplot2)   # For creating plots
library(scales)    # For formatting plot labels (e.g., percent_format)

# Clear the workspace (optional, but good practice)
rm(list = ls())

# --- Setting Parameters ---
cutoff <- 150     # Cutoff value for filtering events (based on number of casualties)
America <- 0      # Flag to filter events by Zone (1 for America, 0 for All)
# Note: fD is not explicitly used in the EDummy creation loop, but the EDummy
# matrix is initialized with ncol = 3, implying fD = 3 future days.

# --- Data Loading ---
# Import industry portfolio returns (assuming it includes Transport industry)
# Use col_names = TRUE to ensure the first row is treated as column names
Ret <- read_excel("../Data/IndustryPortfolios.xlsx", col_names = TRUE)
# Extract dates from the first column and convert to Date objects using ymd()
DatesReturn <- ymd(Ret$Date) # Assuming the date column is named "...1" after import

# Import daily Fama-French Factors (specifically for Risk-Free Rate)
# Use col_names = TRUE to ensure the first row is treated as column names
Rf_data <- read_excel("../Data/FF_FactorsDaily.xlsx", col_names = TRUE)
# Extract the Risk-Free Rate (Rf) column
Rf <- Rf_data$RF

# Import Events data
# Use col_names = TRUE to ensure the first row is treated as column names
Events <- read_excel("../Data/Events.xlsx", col_names = TRUE)

# --- Filter Events ---
# Filter events based on the America flag and casualty cutoff
# Assuming 'Total Deaths' is the column name for casualties and 'Zone' for the zone
if (America == 0) {
  # Filter for all events where 'Total Deaths' is greater than the cutoff
  # Select the 4th column (assumed to be the event date column based on original code)
  Events <- Events[which(Events$`Total Deaths` > cutoff), 4]
} else {
  # Filter for events where 'Total Deaths' is greater than the cutoff AND 'Zone' is 1
  # Select the 4th column (assumed to be the event date column)
  Events <- Events[which(Events$`Total Deaths` > cutoff & Events$Zone == 1), 4]
}

# Extract the Transport industry return series (assuming it's in the column named "Trans")
Ret_asset <- Ret$Trans # Renamed to avoid conflict with the Ret data frame

# Convert filtered event dates to Date objects using dmy() (assuming format "dd-mm-yyyy")
# The column name after filtering is likely `Date start` based on the original code snippet
Events <- dmy(Events$`Date start`)

# --- Create Event Dummies ---
# Create event dummy matrix. Initialize with zeros.
# The matrix will have dimensions: number of observations + 3 (for lags/future days) x 3 (for fD=3 future days)
# Note: The original code initializes with length(DatesReturn) rows, but then adds 3 rows in the loop.
# Let's initialize with length(DatesReturn) + 3 rows to match the loop's indexing.
EDummy <- matrix(0, nrow = length(DatesReturn) + 3, ncol = 3)

# Find the row index in DatesReturn that corresponds to each event date
# `DatesReturn %in% Events` creates a logical vector indicating if each date in DatesReturn is also in Events
Match <- DatesReturn %in% Events

# Create dummies for event days and subsequent days (up to 3 days)
# The dummy EDummy[i+d, d] is 1 if day i is an event day and we are looking at day i+d.
for (i in 1:length(DatesReturn)) { # Loop through indices of DatesReturn
  if (Match[i]) { # If the current date (index i) is an event day
    # Set dummies for the next 3 days (indices 1, 2, 3 relative to event day)
    # Ensure the index (i + d) is within the bounds of the EDummy matrix
    if (i + 1 <= nrow(EDummy)) { # Check for day + 1
      EDummy[i + 1, 1] <- 1
    }
    if (i + 2 <= nrow(EDummy)) { # Check for day + 2
      EDummy[i + 2, 2] <- 1
    }
    if (i + 3 <= nrow(EDummy)) { # Check for day + 3
      EDummy[i + 3, 3] <- 1
    }
  }
  # Note: The original code did not have an else block setting dummies to 0, which is correct
  # since the matrix is initialized with zeros.
}

# Remove the first 3 rows of EDummy as they cannot have event dummies referring to previous events
EDummy <- EDummy[-c(1:3), ] # Adjusting index for 1-based R arrays

# --- Tax Day Dummy (Jan 1–5) ---
# Note: This dummy is created in the original code but not used in the investment strategies below.
# Keeping it for completeness based on the original code.
TaxDay <- rep(0, length(Ret_asset)) # Use length(Ret_asset) for correct size
dday <- day(DatesReturn)
dmonth <- month(DatesReturn)
# Set dummy to 1 if the date is in January (month 1) and between the 1st and 5th day
TaxDay[dmonth == 1 & dday >= 1 & dday <= 5] <- 1

# --- Investment Strategies Analysis ---

# Buy and Hold strategy: Cumulative excess return over time
# Calculate excess returns (Asset Return - Risk-Free Rate)
# Note: Ensure Rf is aligned correctly with Ret_asset. Assuming Rf has the same number of observations.
BuyandHold <- cumsum(Ret_asset - Rf)

# Aviation disasters event-driven strategy
# This strategy seems to modify the excess return on the day after an aviation disaster
# (where EDummy[, 1] == 1).
# The original code calculates `Rf - Ret` on these days.
# A common interpretation of switching to the risk-free asset is that the *excess* return becomes zero on that day.
# If the goal is to show the cumulative *total* return difference, then the return on event+1 day is Rf.
# If the goal is cumulative *excess* return, setting excess return to 0 on those days might be more standard.
# Keeping the original calculation as provided in the user's code.
AviationStrategy <- Ret_asset - Rf # Start with excess returns
# Identify the indices where the dummy for the day after an event is 1 (EDummy[, 1])
aviation_idx <- which(EDummy[, 1] == 1)
# On these days, replace the excess return with `Rf - Ret`.
AviationStrategy[aviation_idx] <- Rf[aviation_idx] - Ret_asset[aviation_idx]

# Calculate the cumulative return for the aviation strategy
AviationStrategy <- cumsum(AviationStrategy)

# --- Plotting the Strategies ---
# Create a data frame for the plot, combining DatesReturn, BuyandHold, and AviationStrategy
# Note: EDummy and AviationStrategy start 3 days after DatesReturn starts due to the EDummy construction.
# We need to align the dates for plotting.
# The EDummy matrix has length(DatesReturn) rows after trimming.
# BuyandHold and AviationStrategy also have length(DatesReturn) rows.
# So, the dates should align.
plot_data <- data.frame(
  DatesReturn = DatesReturn,
  BuyandHold = BuyandHold,
  AviationStrategy = AviationStrategy
)

# Create the ggplot object
Inv_plot <- ggplot(plot_data, aes(x = DatesReturn)) +
  # Add the line for the Buy and Hold strategy
  geom_line(aes(y = BuyandHold, color = "Buy and hold"), linetype = "solid", linewidth = 1, show.legend = TRUE) +
  # Add the line for the Aviation Strategy
  geom_line(aes(y = AviationStrategy, color = "Aviation strategy"), linetype = "dashed", linewidth = 1, show.legend = TRUE) +
  # Set axis labels
  labs(x = "Date", y = "Cumulative return") +
  # Configure x-axis dates (labels every 5 years)
  scale_x_date(date_labels = "%Y-%m-%d", breaks = "5 years") +
  # Apply a minimal theme
  theme_minimal() +
  # Customize x-axis text appearance (angle and alignment)
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Customize legend position and justification
  theme(legend.position = "top", legend.justification = "left") +
  # Format y-axis labels as percentages (scaling by 1 assumes the data is already in decimal)
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  # Define custom colors for the lines in the legend
  scale_color_manual(values = c("Buy and hold" = "black", "Aviation strategy" = "blue")) +
  # Customize legend title
  guides(color = guide_legend(title = "Strategy")) + # Changed legend title to "Strategy"
  # Set plot title
  ggtitle("Investment Strategy")

# Save the plot as a PNG file
ggsave("ES_InvStrategy.jpg", Inv_plot, width = 10, height = 6) # Specify plot object and dimensions

# Print the plot to the console (optional, for display in R environment)
print(Inv_plot)


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
