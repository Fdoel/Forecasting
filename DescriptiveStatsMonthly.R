# Import necessary libraries
library(tidyverse)
library(e1071)
library(ggplot2)

# Load the data
data <- read.csv("FRED.csv")

# Omit the first two rows as we only require raw data
data <- data[-c(1, 2),]

inflation_df <- data %>%
  select(c("sasdate", "CPIAUCSL", "UNRATE", "IPFINAL", "CUMFNS", "RPI", "RETAILx", "VIXCLSx")) %>%
  mutate(sasdate = as.Date(sasdate, "%m/%d/%Y")) %>%
  # Calculate inflation by taking the logs of the CPI divided by its lag
  mutate(inflation = log(CPIAUCSL/lag(CPIAUCSL))*100) %>%
  filter(sasdate >= as.Date("1959-06-01"))
  
# Get summary statistics on all columns except date including Kurtosis, Max, Min, Mean, Median, Skewness, and Standard Deviation
info_df <- inflation_df[-1] %>%
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

# Plot the inflation data with CPi on one y axis and inflation on the other
inflation_df %>%
  ggplot(aes(x = sasdate)) +
  geom_line(aes(y = CPIAUCSL, color = "CPI"), size = 1, linetype = "dashed") +  # Dashed line for CPI
  geom_line(aes(y = inflation * 100, color = "Inflation"), size = 0.7) +  # Scale inflation
  scale_y_continuous(
    name = "CPI",
    sec.axis = sec_axis(~ . / 100, name = "Inflation (%)")  # Scale back for display
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
seastests::kw(inflation_df$inflation, freq = 12)
seastests::seasdum(inflation_df$inflation, freq = 12)
# Do not reject no seasonality at the 5% level

# Test for seasonality in inflation
seastests::kw(inflation_df$CPIAUCSL, freq = 12)
seastests::seasdum(inflation_df$CPIAUCSL, freq = 12)
# Do not reject no seasonality at the 5% level

# Rename the df for use in other scripts
inflation_df_monthly <- inflation_df
save(inflation_df_monthly, file = "inflation_df_monthly.RData")

