# =============================================================================

# Script: US Inflation Data Preparation and Analysis
# Description: 
#   This script prepares monthly US inflation data by merging seasonally 
#   adjusted CPI data from FRED with non-seasonally adjusted labor data. 
#   It calculates multiple inflation measures, tests for seasonality, 
#   produces summary statistics and time series plots, and saves the 
#   processed dataset for further use in forecasting models.
#
# Input:
#   - "CPI US labour dataset.xlsx": Monthly CPI (non-seasonally adjusted)
#   - "FRED.csv": Seasonally adjusted macroeconomic variables including CPI
#
# Output:
#   - inflation_df_monthly.RData: Cleaned and merged monthly inflation dataset
# =============================================================================

# Load required libraries
library(tidyverse)     # For data wrangling and plotting
library(e1071)         # For skewness and kurtosis
library(ggplot2)       # For additional plotting features (already included in tidyverse)
library(readxl)        # For reading Excel files
library(rstudioapi)    # To set working directory to current script location
library(dplyr)
library(tseries)

# Set working directory to the location of the currently opened R script
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the non-seasonally adjusted CPI data from Excel
CPI_US_labour_dataset <- read_excel("CPI US labour dataset.xlsx", 
                                    range = "A12:M124")

# Reshape from wide to long format and create a proper date column
CPI_US_labour_long <- CPI_US_labour_dataset %>%
  # Assume first column is Year, pivot the remaining month columns to long format
  pivot_longer(
    cols = -1,
    names_to = "Month",
    values_to = "CPInonSA"
  ) %>%
  rename(Year = 1) %>%
  filter(Year >= 1959) %>%
  mutate(
    Month_num = match(Month, month.abb),                        # Convert month abbreviation to number
    Date = as.Date(paste(Year, Month_num, "01", sep = "-"))     # Create proper date format
  ) %>%
  arrange(Date)

# Compute month-over-month inflation (non-seasonally adjusted) as log difference in CPI
CPI_US_labour_long <- CPI_US_labour_long %>%
  arrange(Date) %>%
  mutate(inflationNonSA = log(CPInonSA / lag(CPInonSA)) * 100)

# Load the FRED data (seasonally adjusted CPI and other indicators)
data <- read.csv("FRED.csv")

# Omit the first two rows as we only require raw data
data <- data[-c(1, 2),]

inflation_df <- data %>%
  select(c("sasdate", "CPIAUCSL", "UNRATE", "IPFINAL", "CUMFNS", "RPI", "RETAILx", "VIXCLSx","GS1","USGOVT","INDPRO")) %>%
  mutate(sasdate = as.Date(sasdate, "%m/%d/%Y")) %>%
  
  # Calculate inflation by taking the logs of the CPI divided by its lag
  mutate(inflationSA = log(CPIAUCSL/lag(CPIAUCSL))*100) %>%
  mutate(
    ldGS1 = log(GS1) - log(lag(GS1)),
    dCUMFNS = c(NA,diff(CUMFNS)),
    dIPFINAL = c(NA,diff(IPFINAL)),
    dUNRATE = c(NA, diff(UNRATE)),
    dINDPRO = log(INDPRO) - log(lag(INDPRO)),
    dUSGOVT = log(USGOVT) - log(lag(USGOVT)),
    dRETAIL = log(RETAILx) - log(lag(RETAILx)),
    dRPI = log(RPI) - log(lag(RPI))
  ) %>%
  filter(sasdate >= as.Date("1959-06-01"))

# Join with non-seasonally adjusted CPI and inflation data
inflation_df <- inflation_df %>%
  left_join(
    CPI_US_labour_long %>% 
      filter(Date >= as.Date("1959-06-01")) %>% 
      select(Date, CPInonSA, inflationNonSA),
    by = c("sasdate" = "Date")
  )

# Ensure only historical data is included (exclude future observations)
inflation_df <- inflation_df %>%
  filter(sasdate < as.Date("2025-01-01"))

# Generate descriptive statistics (mean, median, sd, etc.) for all columns except date and VIX
info_df <- inflation_df %>%
  select(-sasdate, -VIXCLSx) %>%  # Exclude date and VIXCLSx columns if not needed
  summarise_all(list(
    Mean = mean,
    Median = median,
    SD = sd,
    Min = min,
    Max = max,
    Skewness = ~ skewness(., na.rm = TRUE),
    Kurtosis = ~ kurtosis(., na.rm = TRUE)
  )) %>%
  pivot_longer(
    cols = everything(), 
    names_to = c("Statistic", "Variable"), 
    names_sep = "_"
  ) %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  as.data.frame()

# Generate VIX summary statistics (data available only from July 1962)
VIX_summary <- inflation_df %>%
  select(c("sasdate", "VIXCLSx")) %>%
  filter(sasdate >= as.Date("1962-07-01"))

vix_summary_df <- VIX_summary[-1] %>%
  summarise_all(list(
    Mean = mean,
    Median = median,
    SD = sd,
    Min = min,
    Max = max,
    Skewness = ~ skewness(., na.rm = TRUE),
    Kurtosis = ~ kurtosis(., na.rm = TRUE)
  )) %>%
  pivot_longer(cols = everything(), names_to = c("Statistic", "Variable"), names_sep = "_") %>%
  pivot_wider(names_from = "Variable", values_from = "value") %>%
  as.data.frame()

# -----------------------------------------------------------------------------
# Plot 1: Seasonally Adjusted CPI and Inflation over time
# -----------------------------------------------------------------------------
inflation_df %>%
  ggplot(aes(x = sasdate)) +
  
  # Volcker Disinflation (1979–1987)
  geom_rect(aes(xmin = as.Date("1979-01-01"), xmax = as.Date("1987-12-01"),
                ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
  
  # 2008 Financial Crisis (Sept 2008–June 2009)
  geom_rect(aes(xmin = as.Date("2008-09-01"), xmax = as.Date("2009-06-01"),
                ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
  
  # COVID/Post-COVID (March 2020–Dec 2023)
  geom_rect(aes(xmin = as.Date("2020-03-01"), xmax = as.Date("2023-12-01"),
                ymin = -Inf, ymax = Inf),
            fill = "grey70", alpha = 0.2, inherit.aes = FALSE) +
  
  # Time series lines
  geom_line(aes(y = CPIAUCSL, color = "CPI"), size = 1, linetype = "dashed") +
  geom_line(aes(y = inflationSA * 100, color = "Inflation"), size = 0.7) +
  
  # Axis and labels
  scale_y_continuous(
    name = "CPI",
    sec.axis = sec_axis(~ . / 100, name = "Inflation (%)")
  ) +
  labs(title = "",
       x = "Date",
       color = "Variable") +
  scale_color_manual(values = c("CPI" = "blue", "Inflation" = "red")) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(5, "pt")
  )


# -----------------------------------------------------------------------------
# Plot 2: Non-Seasonally Adjusted CPI and Inflation over time
# -----------------------------------------------------------------------------
inflation_df %>%
  ggplot(aes(x = sasdate)) +
  geom_line(aes(y = CPInonSA, color = "CPI"), size = 1, linetype = "dashed") +  # CPI (dashed)
  geom_line(aes(y = inflationNonSA * 100, color = "Inflation"), size = 0.7) +   # Inflation (scaled again)
  scale_y_continuous(
    name = "CPI",
    sec.axis = sec_axis(~ . / 100, name = "Inflation (%)")  # Adjust for visual scale
  ) +
  labs(title = "",
       x = "Date",
       color = "Variable") +
  scale_color_manual(values = c("CPI" = "blue", "Inflation" = "red")) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(5, "pt") 
  )

# -----------------------------------------------------------------------------
# Seasonality Tests
# -----------------------------------------------------------------------------

# Test for seasonal patterns in monthly inflation (non-seasonally adjusted)
seastests::kw(inflation_df$inflationNonSA, freq = 12)      # Kruskal-Wallis test
seastests::seasdum(inflation_df$inflationNonSA, freq = 12) # Seasonal dummies test

# Test for seasonal patterns in monthly inflation (seasonally adjusted)
seastests::kw(inflation_df$inflationSA, freq = 12)
seastests::seasdum(inflation_df$inflationSA, freq = 12)

# Test for seasonality in CPI indices
seastests::kw(inflation_df$CPInonSA, freq = 12)
seastests::seasdum(inflation_df$CPInonSA, freq = 12)
seastests::kw(inflation_df$CPIAUCSL, freq = 12)
seastests::seasdum(inflation_df$CPIAUCSL, freq = 12)

# -----------------------------------------------------------------------------
# Add Year-over-Year Inflation
# -----------------------------------------------------------------------------

inflation_df <- inflation_df %>%
  arrange(sasdate) %>%
  mutate(
    # Year-over-year seasonal adjusted inflation and previous year's seasonal adjusted CPI
    inflation_yoy_SA = log(CPIAUCSL / lag(CPIAUCSL, 12)) * 100,
    
    # Year-over-year non-seasonally adjusted inflation and previous year's non-seasonally adjusted CPI
    inflation_yoy_nonSA = log(CPInonSA / lag(CPInonSA, 12)) * 100
  )

# -----------------------------------------------------------------------------
# Save Final Monthly Dataset
# -----------------------------------------------------------------------------

# List of variables of interest
vars <- c("CPIAUCSL", "UNRATE", "IPFINAL", "CUMFNS", "RPI", "RETAILx", "VIXCLSx", "inflationNonSA")

# Create inflation column and filter CPI-based start date
inflation_df_monthly <- data %>%
  select(sasdate, all_of(vars[-length(vars)])) %>%
  mutate(
    sasdate = as.Date(sasdate, "%m/%d/%Y"),
    inflation = log(CPIAUCSL / lag(CPIAUCSL)) * 100
  ) %>%
  filter(sasdate >= as.Date("1959-06-01"))

# Filter VIX from its available start date
inflation_df <- inflation_df %>%
  mutate(VIXCLSx = ifelse(sasdate < as.Date("1962-07-01"), NA, VIXCLSx))

# Compute summary stats
summary_stats <- inflation_df %>%
  select(all_of(vars)) %>%
  summarise(across(
    everything(),
    list(
      Mean = ~ mean(., na.rm = TRUE),
      Median = ~ median(., na.rm = TRUE),
      SD = ~ sd(., na.rm = TRUE),
      Skewness = ~ skewness(., na.rm = TRUE),
      Kurtosis = ~ kurtosis(., na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  )) %>%
  pivot_longer(everything(), names_to = c("Variable", "Statistic"), names_sep = "_") %>%
  pivot_wider(names_from = Statistic, values_from = value)

summary_stats_rounded <- summary_stats %>%
  mutate(across(where(is.numeric), ~ sprintf("%.2f", .)))

print(summary_stats_rounded)

library(ggplot2)

# Compute ACF with lag.max = 28
acf_obj <- acf(inflation_df$inflationNonSA, lag.max = 28, plot = FALSE, na.action = na.pass)

# Create data frame for plotting (excluding lag 0)
acf_df <- data.frame(
  Lag = 1:28,
  ACF = acf_obj$acf[2:29]
)

# Confidence interval (95%)
conf_level <- 1.96 / sqrt(length(na.omit(inflation_df$inflationNonSA)))

# Plot with x-axis ticks every 4 lags
ggplot(acf_df, aes(x = Lag, y = ACF)) +
  geom_col(fill = "steelblue", width = 0.7) +
  geom_hline(yintercept = 0, color = "black") +
  geom_hline(yintercept = c(-conf_level, conf_level), linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = seq(4, 28, 4)) +
  labs(title = "",
       x = "Lag (Months)",
       y = "Autocorrelation") +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 10)
  )

# Ljung-Box test for autocorrelation up to lag 28
Box.test(
  inflation_df$inflationNonSA,
  lag = 28,
  type = "Ljung-Box"
)

# ADF test on non-seasonally adjusted inflation
adf.test(inflation_df$inflationNonSA, alternative = "stationary")

# -----------------------------------------------------------------------------
# Test for correlation 
# -----------------------------------------------------------------------------


# Load the FRED data (seasonally adjusted CPI and other indicators)
data <- read.csv("FRED.csv")

# Omit the first two rows as we only require raw data
data <- data[-c(1, 2),]

inflation_df_cor <- data %>%
  mutate(sasdate = as.Date(sasdate, "%m/%d/%Y")) %>%
  
  # Calculate inflation by taking the logs of the CPI divided by its lag
  mutate(inflationSA = log(CPIAUCSL/lag(CPIAUCSL))*100) %>%
  filter(sasdate >= as.Date("1959-06-01"))

# Join with non-seasonally adjusted CPI and inflation data
inflation_df_cor <- inflation_df_cor %>%
  left_join(
    CPI_US_labour_long %>% 
      filter(Date >= as.Date("1959-06-01")) %>% 
      select(Date, CPInonSA, inflationNonSA),
    by = c("sasdate" = "Date")
  )

# Ensure only historical data is included (exclude future observations)
inflation_df_cor <- inflation_df_cor %>%
  filter(sasdate < as.Date("2025-01-01"))


# Compute correlation matrix for numeric variables (excluding date)
cor_matrix <- inflation_df_cor %>%
  select(c("inflationNonSA","UNRATE","IPFINAL","CUMFNS","RPI","RETAILx","GS1","USGOVT","INDPRO")) %>%
  select_if(is.numeric) %>%
  cor(use = "pairwise.complete.obs")

# Extract correlations with inflationNonSA
cor_with_inflationNonSA <- cor_matrix["inflationNonSA", ]

# Display all correlations with inflationNonSA
cat("All correlations with inflationNonSA:\n")
print(round(cor_with_inflationNonSA, 3))

# Filter for absolute correlation > 0.85 (excluding inflationNonSA itself)
high_corr_vars <- cor_with_inflationNonSA[abs(cor_with_inflationNonSA) > abs(0.170) & names(cor_with_inflationNonSA) != "inflationNonSA"]

# Display strong correlations
cat("\nCorrelations with inflationNonSA above 0.85 or below -0.85:\n")
print(round(high_corr_vars, 3))

# Rename the final dataframe for use in forecasting and save
inflation_df_monthly <- inflation_df
save(inflation_df_monthly, file = "inflation_df_monthly.RData")


cor(inflation_df_monthly$ldGS1, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dCUMFNS, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dIPFINAL, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dUNRATE, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dINDPRO, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dRETAIL, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dUSGOVT, inflation_df_monthly$inflationNonSA)
cor(inflation_df_monthly$dRPI, inflation_df_monthly$inflationNonSA)

# -----------------------------------------------------------------------------
# Plot the distribution of non-SA inflation
# -----------------------------------------------------------------------------

# Compute mean and sd for the normal curve
mean_val <- mean(inflation_df_monthly$inflationNonSA, na.rm = TRUE)
sd_val <- sd(inflation_df_monthly$inflationNonSA, na.rm = TRUE)

# Plot using ggplot2
ggplot(df, aes(x = inflation)) +
  geom_histogram(aes(y = ..density..), 
                 bins = 30, 
                 fill = "lightblue", 
                 color = "white") +
  stat_function(fun = dnorm, 
                args = list(mean = mean_val, sd = sd_val), 
                color = "red", 
                size = 1) +
  labs(title = "Histogram of Non-SA Inflation with Normal Distribution",
       x = "Inflation (Non-SA)",
       y = "Density") +
  theme_minimal()

# For each threshold value from 0.0
