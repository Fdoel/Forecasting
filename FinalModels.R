# Load data and functions
load("inflation_df_monthly.RData")
source("MART.R")
source("MARX_functions.R")

# Define c for pseudo models
c_med <- median(inflation_df_monthly$inflationNonSA)

exo <- cbind(inflation_df_monthly$ldGS1, inflation_df_monthly$dRPI, inflation_df_monthly$dRETAIL)

# AR models
AR <- arx.ls(inflation_df_monthly$inflationNonSA,NULL,12)
AR_x <- arx.ls(inflation_df_monthly$inflationNonSA, exo, 12)

# SETAR models
SETAR <- MART(inflation_df_monthly$inflationNonSA, NULL, 2, 0, 0.6,1) # SETAR using grid search
SETAR_x <- MART(inflation_df_monthly$inflationNonSA, exo, 2, 0, 0.3 ,1) # SETAR-X using grid search

# STAR models
STAR <- SMART(inflation_df_monthly$inflationNonSA, NULL, 2,0, 0.6,1,1) # STAR using grid search
STAR_x <- SMART(inflation_df_monthly$inflationNonSA, exo, 2,0, 0.3,1,1) # STAR-X using grid search

# All MAR models
MAR_grid_pseudo <- marx.t(inflation_df_monthly$inflationNonSA, NULL, 1,12) # MAR is the same for grid search and pseudo
MAR_x <- marx.t(inflation_df_monthly$inflationNonSA, exo, 1,12) # MAR-X using grid search DEZE MOET ER NOG IN

# All MART models
MART_pseudo <- MART(inflation_df_monthly$inflationNonSA, NULL, 1,1, c_med, d=3)
MART_grid <- MART(inflation_df_monthly$inflationNonSA, NULL, 1,1, c=0.6, d=4)
MART_x_grid <- MART(inflation_df_monthly$inflationNonSA, exo, 1,3, c=0.3, d=1)
MART_x_pseudo <- MART(inflation_df_monthly$inflationNonSA, exo, 1,4, c=c_med, d=1)

# SMART models
SMART_grid <- SMART(inflation_df_monthly$inflationNonSA, NULL, 1,1, c=0.6, gamma = 3, d=3)
SMART_pseudo <- SMART(inflation_df_monthly$inflationNonSA, NULL, 1,1, c=c_med, gamma = 2, d=3)
