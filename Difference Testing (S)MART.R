## -----------------------------------------------------------------------------
# Testing if the coefficients differ significantly in the regimes F test
# -----------------------------------------------------------------------------
# Compare AR to ART model
model.ar <- arx.ls(inflation_df_monthly$inflationNonSA, NULL, 12)
model.art <- MART(inflation_df_monthly$inflationNonSA, NULL, 12, 0, 0.6,1)

model.mart <- MART(inflation_df_monthly$inflationNonSA, NULL,1,1, 0.6, d=4)
model.mar <- marx.t(inflation_df_monthly$inflationNonSA, NULL, 1,1)
mart_resid_var <- var(model.mart$residuals)
mar_resid_var <- var(model.mar$residuals)

t <- length(mart_resid_var)
F_c <- t * (mar_resid_var - mart_resid_var) / mart_resid_var
p_value_ar <- 1 - pchisq(F_c, df = 2 + 1)


# Compare MART and MAR models
model.mart <- MART(inflation_df_monthly$inflationNonSA, NULL,1,11, 0.6, d=4)
model.mar <- marx.t(inflation_df_monthly$inflationNonSA, NULL, 1,11)
mart_resid_var <- var(model.mart$residuals)
mar_resid_var <- var(model.mar$residuals)

t <- length(mart_resid_var)

F_c <- t * (mar_resid_var - mart_resid_var) / mart_resid_var
p_value_art <- 1 - pchisq(F_c, df = 2 + 1)

print(paste("P-value for AR vs SETAR:", p_value_ar))
print(paste("P-value for MART vs MAR:", p_value_art))


