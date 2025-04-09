# Import necessary libraries
library(tidyverse)
library(e1071)
library(ggplot2)
library(readxl)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
    Month_num = match(Month, month.abb),
    Date = as.Date(paste(Year, Month_num, "01", sep = "-"))
  ) %>%
  arrange(Date)

# Add an inflation column computed as log(CPI(t)/CPI(t-1)) * 100
CPI_US_labour_long <- CPI_US_labour_long %>%
  arrange(Date) %>%
  mutate(inflationNonSA = log(CPInonSA / lag(CPInonSA)) * 100)

# Load the FRED data
FRED_data <- read.csv("FRED.csv")

# Omit the first two rows as we only require raw data
FRED_data <- FRED_data[-c(1, 2),]

inflation_df <- FRED_data %>%
  select(c("sasdate", "CPIAUCSL", "UNRATE", "IPFINAL", "CUMFNS", "RPI", "RETAILx", "VIXCLSx")) %>%
  mutate(sasdate = as.Date(sasdate, "%m/%d/%Y")) %>%
  
  # Calculate inflation by taking the logs of the CPI divided by its lag
  mutate(inflationSA = log(CPIAUCSL/lag(CPIAUCSL))*100) %>%
  filter(sasdate >= as.Date("1959-06-01"))

inflation_df <- inflation_df %>%
  left_join(
    CPI_US_labour_long %>% 
      filter(Date >= as.Date("1959-06-01")) %>% 
      select(Date, CPInonSA, inflationNonSA),
    by = c("sasdate" = "Date")
  )

#ensure all data is available
inflation_df <- inflation_df %>%
  filter(sasdate < as.Date("2025-01-01"))

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

VIX_summary <- inflation_df %>%
  select(c("sasdate", "VIXCLSx")) %>%
  filter(sasdate >= as.Date("1962-07-01")) #vix is only available from july 1962

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

# Plot the SA inflation data with SA CPI on one y axis and inflation on the other
inflation_df %>%
  ggplot(aes(x = sasdate)) +
  geom_line(aes(y = CPIAUCSL, color = "CPI"), size = 1, linetype = "dashed") +  # Dashed line for CPI
  geom_line(aes(y = inflationSA * 100, color = "Inflation"), size = 0.7) +  # Scale inflation
  scale_y_continuous(
    name = "CPI (Seasonally Adjusted)",
    sec.axis = sec_axis(~ . / 100, name = "Inflation (%) (Seasonally Adjusted)")  # Scale back for display
  ) +
  labs(title = "",
       x = "Date",
       color = "Variable") +
  scale_color_manual(values = c("CPI" = "blue", "Inflation" = "red")) +  # Define colors properly
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(5, "pt"), 
  )

# Plot the non SA inflation data with CPi on one y axis and inflation on the other
inflation_df %>%
  ggplot(aes(x = sasdate)) +
  geom_line(aes(y = CPInonSA, color = "CPI"), size = 1, linetype = "dashed") +  # Dashed line for CPI
  geom_line(aes(y = inflationNonSA * 100, color = "Inflation"), size = 0.7) +  # Scale inflation
  scale_y_continuous(
    name = "CPI (raw)",
    sec.axis = sec_axis(~ . / 100, name = "Inflation (%) (raw)")  # Scale back for display
  ) +
  labs(title = "",
       x = "Date",
       color = "Variable") +
  scale_color_manual(values = c("CPI" = "blue", "Inflation" = "red")) +  # Define colors properly
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black"), 
    axis.ticks.length = unit(5, "pt"), 
  )


# Test for seasonality in inflation
seastests::kw(inflation_df$inflationNonSA, freq = 12)
seastests::seasdum(inflation_df$inflationNonSA, freq = 12)
seastests::kw(inflation_df$inflationSA, freq = 12)
seastests::seasdum(inflation_df$inflationSA, freq = 12)

# Test for seasonality in CPI
seastests::kw(inflation_df$CPInonSA, freq = 12)
seastests::seasdum(inflation_df$CPInonSA, freq = 12)
seastests::kw(inflation_df$CPIAUCSL, freq = 12)
seastests::seasdum(inflation_df$CPIAUCSL, freq = 12)

#add inflation with regard to year previously
inflation_df <- inflation_df %>%
  arrange(sasdate) %>%
  mutate(
    # Year-over-year seasonal adjusted inflation and previous year's seasonal adjusted CPI
    inflation_yoy_SA = log(CPIAUCSL / lag(CPIAUCSL, 12)) * 100,
    
    # Year-over-year non-seasonally adjusted inflation and previous year's non-seasonally adjusted CPI
    inflation_yoy_nonSA = log(CPInonSA / lag(CPInonSA, 12)) * 100,
  )

# Rename the df for use in other scripts
inflation_df_monthly <- inflation_df
save(inflation_df_monthly, file = "inflation_df_monthly.RData")



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
