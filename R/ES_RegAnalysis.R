# ES_RegAnalysis.R
# This script performs regression analysis on industry portfolio returns,
# including event studies and calendar effects (weekday and tax days) in R.
#
# Developed for: Essentials of Financial Economics
# Authors: Michael Donadelli, Michele Costola, Ivan Gufler
# Date: May 8, 2025

# Set locale for consistent weekday naming, especially important for `weekdays()`
Sys.setlocale("LC_TIME", "en_US")

# Load required libraries
library(readxl)      # For reading Excel files
library(dplyr)       # For data manipulation (e.g., lag, select)
library(lubridate)   # For working with dates (ymd, day, month, weekdays)
library(sandwich)    # For robust standard errors (like Newey-West, though not explicitly used in the provided code's SE calculation)
library(lmtest)      # For testing linear regression models (e.g., coeftest, though not explicitly used)
# Note: The provided code calculates standard errors manually after fitting lm.

# --- Settings ---
maxlag <- 3       # Maximum lag for lagged returns in the regression
cutoff <- 150     # Cutoff value for filtering events (based on number of casualties)
fD <- 3           # Number of future days to create event dummies for
cons <- 1         # Flag: 1 to include a constant (intercept) in the regression, 0 to exclude
America <- 1      # Flag: 1 to filter events by Zone (America), 0 for all events

# --- Import Data ---

# Read industry portfolio returns and dates from an Excel file
# Use col_names = TRUE to ensure the first row is treated as column names
Ret <- read_excel("../Data/IndustryPortfolios.xlsx", col_names = TRUE)
# Extract dates from the first column and convert to Date objects using ymd()
DatesReturn <- ymd(Ret$Date) # Assuming the date column is named "...1" after import

# Read events data from an Excel file
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
ExRet <- Ret$Trans

# Convert filtered event dates to Date objects using dmy() (assuming format "dd-mm-yyyy")
# The column name after filtering is likely `Date start` based on the original code snippet
Events <- dmy(Events$`Date start`)

# --- Weekday Dummies (Mon to Thu) ---
# Preallocate matrix for weekday dummies (Monday=1, Tuesday=2, ..., Thursday=4)
# The matrix will have dimensions: number of observations x 4
dow <- matrix(NA, nrow = length(DatesReturn), ncol = 4)

# Create dummy variables for Monday (weekday 2 in R), Tuesday (3), Wednesday (4), Thursday (5)
# Note: The original code comments mention weekday 2=Mon, ..., 5=Thu, but the loop uses d=1:4.
# R's weekdays() function returns "Monday", "Tuesday", etc. Let's match the string names.
days <- c("Monday", "Tuesday", "Wednesday", "Thursday")
for (d in 1:4) { # Loop for Monday (index 1 in 'days') to Thursday (index 4 in 'days')
  for (i in 1:length(DatesReturn)) { # Loop through each observation date
    # Check if the weekday of the current date matches the current dummy day (days[d])
    if (weekdays(DatesReturn[i]) == days[d]) {
      dow[i, d] <- 1 # Set dummy to 1 if it's that weekday
    } else {
      dow[i, d] <- 0 # Set dummy to 0 otherwise
    }
  }
}

# --- Event Dummies ---
# Create event dummy matrix. Initialize with zeros.
# The matrix will have dimensions: number of observations + 3 (for lags) x 3 (for fD=3 future days)
EDummy <- matrix(0, nrow = length(DatesReturn) + 3, ncol = 3)

# Find the row index in DatesReturn that corresponds to each event date
# `DatesReturn %in% Events` creates a logical vector indicating if each date in DatesReturn is also in Events
Match <- DatesReturn %in% Events

# Create dummies for event days and subsequent days (up to fD)
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
  # Note: The original code also had else blocks setting dummies to 0, which is redundant
  # since the matrix is initialized with zeros. Removed the redundant else blocks.
}

# Remove the first 3 rows of EDummy as they cannot have event dummies referring to previous events
EDummy <- EDummy[-c(1:3), ]

# --- Tax Day Dummy (Jan 1–5) ---
# Initialize array for tax day dummy with zeros
TaxDay <- rep(0, length(DatesReturn)) # Use length(DatesReturn) for correct size
# Extract day and month from return dates
dday <- day(DatesReturn)
dmonth <- month(DatesReturn)
# Set dummy to 1 if the date is in January (month 1) and between the 1st and 5th day
TaxDay[dmonth == 1 & dday >= 1 & dday <= 5] <- 1
# Replace any NA values in TaxDay with 0 (though lubridate functions usually don't produce NA here)
TaxDay[is.na(TaxDay)] = 0

# --- Lagged Returns ---
# Create lagged returns matrix using dplyr's lag() function
# lag() introduces NA values at the beginning for the lagged periods
RetLag <- cbind(lag(ExRet, n = 1), lag(ExRet, n = 2), lag(ExRet, n = 3))

# Handle initial NA values in RetLag (introduced by lag)
# The original code sets these to 0. An alternative is to remove the corresponding rows later.
# Keeping the original approach of setting to 0 for consistency.
RetLag[1:maxlag, ] <- 0

# --- Design Matrix ---
# Combine all independent variables into a single matrix X
# Include lagged returns, weekday dummies, tax day dummy, and event dummies (up to fD)
X <- cbind(RetLag, dow, TaxDay, EDummy[, 1:fD])

# --- Regression ---
# Fit the linear regression model using lm()
# The dependent variable is ExRet. The independent variables are in matrix X.
# The 'cons' flag determines whether to include an intercept.
if (cons == 1) {
  # Include intercept: ExRet ~ X + intercept
  mdl <- lm(ExRet ~ X)
} else {
  # Exclude intercept: ExRet ~ X - 1
  mdl <- lm(ExRet ~ X - 1)
}

# Extract coefficients and standard errors from the regression model
coeff <- round(coef(mdl), 5) # Coefficients, rounded to 5 decimal places
se <- sqrt(diag(vcov(mdl))) # Standard errors (using vcov to get variance-covariance matrix)

# Calculate t-statistics (Coefficient / Standard Error)
Sig <- coeff / se

# Combine coefficients and t-statistics into a single matrix for the output table
# Use rbind to stack the coefficient and t-statistic vectors as rows
res <- matrix(rbind(coeff, Sig), nrow = 2)

# --- Variable Names for Output Table ---
# Define variable names for the output table columns
# Names for lagged returns
rString <- paste("R_t-", 1:maxlag, sep = "")
# Names for weekday dummies
varnames_dow <- c("Mon", "Tue", "Wed", "Thu")
# Name for tax day dummy
varnames_tax <- c("TaxDays")
# Names for event dummies
eString <- paste("Event +", 1:fD)
# Combine all predictor variable names
predictor_vars <- c(rString, varnames_dow, varnames_tax, eString)

# Add "Constant" to the variable names if the intercept was included
if (cons == 1) {
  var <- c("Constant", predictor_vars)
} else {
  var <- predictor_vars # No constant if excluded
}

# --- Output Table ---
# Create a data frame to display the regression results
# The data frame will have 2 rows ('Coefficient', 't-Statistic') and columns for each variable
tab <- data.frame(matrix(NA, ncol = length(var), nrow = 2))
# Set column names to the defined variable names
colnames(tab) <- var
# Set row names to identify the statistics
rownames(tab) <- c("Coefficient", "t-Statistic")
# Populate the data frame with the results matrix (res)
tab[] <- res # Assign the results matrix to the data frame

# --- Save Output ---
# Define the output file name based on the America flag
output_file <- if (America == 0) "RegressionAnalysis_All.csv" else "RegressionAnalysis_America.csv"
# Write the results table to a CSV file
# Use write.table for more control over separators and row/column names
write.table(tab, output_file, sep = ",", col.names = TRUE, row.names = TRUE)


# This R script contains code examples and material from the book:
# "Essentials of Financial Economics" by Michael Donadelli, Michele Costola, and Ivan Gufler.
# You can find more information and download the book at: link.springer.com/book/9783031861895
