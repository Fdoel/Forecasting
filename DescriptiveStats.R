# Import necessary libraries
library(tidyverse)
library(e1071)

# Load the data
data <- read.csv("FREDQ.csv")

# Omit the first two rows as we only require raw data
data <- data[-c(1, 2),]

inflation_df <- data %>%
  select(c("sasdate", "CPIAUCSL")) %>%
  mutate(sasdate = as.Date(sasdate, "%m/%d/%Y")) %>%
  # Calculate inflation by taking the logs of the CPI divided by its lag
  mutate(inflation = log(CPIAUCSL/lag(CPIAUCSL))*100) %>%
  filter(sasdate >= "1980-1-1")
  
# Get summary statistics on inflation including Kurtosis, Max, Min, Mean, Median, Skewness, and Standard Deviation
info_df <- rbind(Kurtosis = kurtosis(inflation_df$inflation, type = 2),
              Max = max(inflation_df$inflation),
              Min = min(inflation_df$inflation),
              Mean = mean(inflation_df$inflation),
              Median = median(inflation_df$inflation),
              Skewness = skewness(inflation_df$inflation),
              Standard_Deviation = sd(inflation_df$inflation))

colnames(info_df) <- "Inflation"

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

# Test for seasonality in 
# Test for seasonality
seastests::kw(inflation_df$inflation, freq = 4)
seastests::seasdum(inflation_df$inflation, freq = 4)
# Do not reject no seasonality at the 5% level
