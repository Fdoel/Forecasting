# Load required library
library(MASS)
source("MARX_functions.R")


#gwn interface volgen
marx(inflation_df_monthly$inflationNonSA, NULL, p_max=6, sig_level = 0.1)

marx_test1 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=0, p_NC=4)
marx_test2 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=1, p_NC=3)
marx_test3 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=2, p_NC=2)
marx_test4 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=3, p_NC=1)
marx_test5 <- mixed(inflation_df_monthly$inflationNonSA, NULL, p_C=4, p_NC=0)

forecast <- forecast.marx(y=inflation_df_monthly$inflationNonSA, X=NULL, p_C=1, p_NC=3, h=1, M=50, N=10000)$Y.for








