# =============================================================================
# Script: US Inflation Data Preparation and Analysis (Quarterly Version)
# Description: 
#   This script prepares US inflation data by merging seasonally adjusted CPI 
#   data from FRED with non-seasonally adjusted labor data. It calculates 
#   multiple inflation measures, converts to quarterly frequency with 
#   end-of-quarter dates, and saves the processed dataset for further use.
#
# Input:
#   - "CPI US labour dataset.xlsx": Monthly CPI (non-seasonally adjusted)
#   - "FRED.csv": Seasonally adjusted macroeconomic variables including CPI
#
# Output:
#   - inflation_df_quarterly.RData: Cleaned and merged quarterly dataset
# =============================================================================

# Load required libraries
library(tidyverse)
library(e1071)
library(ggplot2)
library(readxl)
library(rstudioapi)
library(dplyr)
library(tseries)
library(lubridate)

# Load the non-seasonally adjusted CPI data from Excel
CPI_US_labour_dataset <- read_excel("CPI US labour dataset.xlsx", 
                                    range = "A12:M124")

# Reshape from wide to long format and create a proper date column
CPI_US_labour_long <- CPI_US_labour_dataset %>%
  pivot_longer(
    cols = -1,
    names_to = "Month",
    values_to = "CPInonSA"
  ) %>%
  rename(Year = 1) %>%
  filter(Year >= 1959) %>%
  mutate(
    Month_num = match(Month, month.abb),
    Date = as.Date(paste(Year, Month_num, "01", sep = "-"))
  ) %>%
  arrange(Date)

# Load the FRED data (seasonally adjusted CPI and other indicators)
data <- read.csv("FRED.csv")
data <- data[-c(1, 2),]

inflation_df <- data %>%
  select(c("sasdate", "CPIAUCSL", "UNRATE", "IPFINAL", "CUMFNS", "RPI", 
           "RETAILx", "VIXCLSx", "GS1", "USGOVT", "INDPRO")) %>%
  mutate(sasdate = as.Date(sasdate, "%m/%d/%Y")) %>%
  mutate(
    ldGS1 = log(GS1) - log(lag(GS1)),
    dCUMFNS = c(NA, diff(CUMFNS)),
    dIPFINAL = c(NA, diff(IPFINAL)),
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
      select(Date, CPInonSA),
    by = c("sasdate" = "Date")
  ) %>%
  filter(sasdate < as.Date("2025-01-01"))

# Add quarter label
inflation_df <- inflation_df %>%
  mutate(Quarter = quarter(sasdate),
         Year = year(sasdate),
         Quarter_end = ceiling_date(sasdate, unit = "quarter") - days(1))

# Calculate quarterly inflation based on CPI log difference (start vs end of quarter)
# *** MODIFIED ***
inflation_df_quarterly <- inflation_df %>%
  group_by(Year, Quarter) %>%
  arrange(sasdate) %>%
  summarise(
    Quarter = max(Quarter_end),
    CPIAUCSL_start = first(CPIAUCSL),
    CPIAUCSL_end = last(CPIAUCSL),
    CPInonSA_start = first(CPInonSA),
    CPInonSA_end = last(CPInonSA),
    
    # Quarterly inflation: log diff of CPI end vs. start
    inflationSA = 100 * log(CPIAUCSL_end / CPIAUCSL_start),
    inflationNonSA = 100 * log(CPInonSA_end / CPInonSA_start),
    
    # End-of-quarter and averaged values
    UNRATE = mean(UNRATE, na.rm = TRUE),
    IPFINAL = mean(IPFINAL, na.rm = TRUE),
    CUMFNS = mean(CUMFNS, na.rm = TRUE),
    RPI = mean(RPI, na.rm = TRUE),
    RETAILx = mean(RETAILx, na.rm = TRUE),
    VIXCLSx = mean(VIXCLSx, na.rm = TRUE),
    GS1 = mean(GS1, na.rm = TRUE),
    USGOVT = mean(USGOVT, na.rm = TRUE),
    INDPRO = mean(INDPRO, na.rm = TRUE),
    
    # Quarterly changes
    ldGS1 = sum(ldGS1, na.rm = TRUE),
    dCUMFNS = sum(dCUMFNS, na.rm = TRUE),
    dIPFINAL = sum(dIPFINAL, na.rm = TRUE),
    dUNRATE = sum(dUNRATE, na.rm = TRUE),
    dINDPRO = sum(dINDPRO, na.rm = TRUE),
    dUSGOVT = sum(dUSGOVT, na.rm = TRUE),
    dRETAIL = sum(dRETAIL, na.rm = TRUE),
    dRPI = sum(dRPI, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(Quarter)

# Optional: Save the dataset
# save(inflation_df_quarterly, file = "inflation_df_quarterly.RData")

# =============================================================================
# Script: MAR Model Estimation and Forecast Evaluation
# Description:
#   This script estimates and compares Mixed Causal-Noncausal Autoregressive (MAR)
#   models for US inflation using non-seasonally adjusted data. It includes:
#     - Model selection via information criteria
#     - Residual diagnostics for AR models
#     - Forecasting with various model configurations
#     - Out-of-sample forecast performance evaluation using RMSE
# =============================================================================

# Load required libraries
library(MASS)          # For statistical distributions and matrix functions
source("MARX_functions.R")  # Custom functions for MAR model estimation
source("MART.R")            # MART model training and forecasting routines (includes information criteria calculations)
library(forecast)      # For ARIMA modeling and forecast tools
library(pbmcapply)     # For parallel processing with progress bar
library(stats)

#basic MAR
marx(inflation_df_quarterly$inflationNonSA, NULL, p_max = 6, sig_level = 0.1)

# Define maximum lag orders for causal (p_C) and noncausal (p_NC) components
p_C_max <- 6
p_NC_max <- 6

# Initialize matrices to store information criteria and model diagnostics
LL <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
AIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
BIC <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
HQ <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)
N <- matrix(NA, nrow = p_C_max + 1, ncol = p_NC_max + 1)

# Loop over all possible (p_C, p_NC) combinations and compute info criteria
for (i in 0:p_C_max) {
  for (j in 0:p_NC_max) {
    marx_loop <- marx.t(inflation_df_quarterly$inflationNonSA, NULL, p_C = i, p_NC = j)
    information <- information.criteria(type = "MARX", marx_loop)
    LL[(i+1),(j+1)] <- information$loglikelihood
    AIC[(i+1),(j+1)] <- information$aic
    BIC[(i+1),(j+1)] <- information$bic
    HQ[(i+1),(j+1)] <- information$hq
    N[(i+1),(j+1)] <- information$n
  }
}