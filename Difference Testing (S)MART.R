## -----------------------------------------------------------------------------
# Testing if the coefficients differ significantly in the regimes F test
# -----------------------------------------------------------------------------
model.MART <- MART(inflation_df_monthly$inflationNonSA,NULL,1,1,median(inflation_df_monthly$inflationNonSA),3)
model.ART <- MART(inflation_df_monthly$inflationNonSA,NULL,2,0,median(inflation_df_monthly$inflationNonSA),3)

mart_resid_var <- var(model.MART$residuals)
art_resid_var <- var(model.ART$residuals)

t <- length(mart_resid_var)

F_c <- t * (art_resid_var - mart_resid_var) / mart_resid_var

p_value <- 1 - pchisq(F_c, df = 2 + 1)
p_value



model.SMART <- SMART(inflation_df_monthly$inflationNonSA,NULL,1,1,median(inflation_df_monthly$inflationNonSA),3,3)
model.STAR <- SMART(inflation_df_monthly$inflationNonSA,NULL,2,0,0.6,1,1)

