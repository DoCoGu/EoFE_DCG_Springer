# Load required libraries
library(readxl)
library(dplyr)
library(lubridate)
library(ggplot2)

# Setting
cutoff <- 150
America <- 1

# Import Data
Ret <- read_excel("../Data/IndustryPortfolios.xlsx", col_names = TRUE)
Rf <- read_excel("../Data/FF_FactorsDaily.xlsx", col_names = TRUE)
Rf <- Rf$RF
DatesReturn <- ymd(Ret$...1)

Events <- read_excel("../Data/Events.xlsx", col_names = TRUE)

if (America == 0) {
  Events <- Events[which(Events$`Total Deaths` > cutoff),4]
} else {
  Events <- Events[which(Events$`Total Deaths` > cutoff & Events$Zone==1),4]
}

Ret <- Ret$Trans

# Convert dates to Date format
Events <- dmy(Events$`Date start`)

# Create events dummy
Match <- DatesReturn %in% Events

EDummy <- matrix(0, nrow = length(DatesReturn), ncol = 3)

for (i in 1:length(DatesReturn)) {
  if (Match[i]) {
    EDummy[i + 1, 1] <- 1
    EDummy[i + 2, 2] <- 1
    EDummy[i + 3, 3] <- 1
  }
}

EDummy <- EDummy[-c(1:3), ]

T <- rep(0, length(Ret))
dday <- day(DatesReturn)
dmonth <- month(DatesReturn)

T[dmonth == 1 & dday >= 1 & dday <= 5] <- 1

# Investment Strategy
# Buy and Hold strategy
BuyandHold <- cumsum(Ret - Rf)

# Aviation disasters
AviationStrategy <- Ret - Rf
AviationStrategy[EDummy[, 1] == 1] <- Rf[EDummy[, 1] == 1] - Ret[EDummy[, 1] == 1]
AviationStrategy <- cumsum(AviationStrategy)

# Create a data frame for the plot
plot_data <- data.frame(DatesReturn = DatesReturn, BuyandHold = BuyandHold, AviationStrategy = AviationStrategy)

# Create the ggplot
Inv_plot <- ggplot(plot_data, aes(x = DatesReturn)) +
  geom_line(aes(y = BuyandHold, color = "Buy and hold"), linetype = "solid", linewidth = 1, show.legend = TRUE) +
  geom_line(aes(y = AviationStrategy, color = "Aviation strategy"), linetype = "dashed", linewidth = 1, show.legend = TRUE) +
  labs(x = "Date", y = "Cumulative return") +
  scale_x_date(date_labels = "%Y-%m-%d", breaks="5 years") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="top",legend.justification="left")+ 
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  scale_color_manual(values = c("Buy and hold" = "black", "Aviation strategy" = "blue")) +
  guides(color = guide_legend(title = "Legend")) +
  ggtitle("Investment Strategy")

# Set plot size
ggsave("ES_InvStrategy.png", Inv_plot, width = 10, height = 6)

# Print the plot
print(Inv_plot)

