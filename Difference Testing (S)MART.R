## -----------------------------------------------------------------------------
# Testing if the residual variance differ significantly in the regimes F test
# -----------------------------------------------------------------------------
# Compare AR to ART model
model.ar <- marx.t(inflation_df_monthly$inflationNonSA, NULL, 12,0)
model.art <- MART(inflation_df_monthly$inflationNonSA, NULL, 12, 0, 0.6,1)

art_resid_var <- var(model.art$residuals)
ar_resid_var <- var(model.ar$residuals)

t <- length(mart_resid_var)
F_c <- t * (ar_resid_var - art_resid_var) / art_resid_var
p_value_ar <- 1 - pchisq(F_c, df = 2 + 1)


# Compare MART and MAR models
model.mart <- MART(inflation_df_monthly$inflationNonSA, NULL,1,11, 0.6, d=4)
model.mar <- marx.t(inflation_df_monthly$inflationNonSA, NULL, 1,11)
mart_resid_var <- var(model.mart$residuals)
mar_resid_var <- var(model.mar$residuals)

t <- length(mart_resid_var)

F_c <- t * (mar_resid_var - mart_resid_var) / mart_resid_var
p_value_art <- 1 - pchisq(F_c, df = 2 + 1)



# MAR vs. SMART
model.smart <- SMART(inflation_df_monthly$inflationNonSA, NULL, 1,11, c=0.6, gamma = 3, d=3)
smart_resid_var <- var(model.smart$residuals)

t <- length(smart_resid_var)

F_c <- t * (mar_resid_var - smart_resid_var) / smart_resid_var
p_value_smart <- 1 - pchisq(F_c, df = 2 + 1)

# STAR vs AR
model.star <- SMART(inflation_df_monthly$inflationNonSA, NULL, 12,0, 0.6,1,1) # STAR using grid search
star_resid_var <- var(model.star$residuals)
t <- length(star_resid_var)
F_c <- t * (ar_resid_var - star_resid_var) / star_resid_var
p_value_star <- 1 - pchisq(F_c, df = 2 + 1)

print(paste("P-value for AR vs SETAR:", p_value_ar))
print(paste("P-value for MART vs MAR:", p_value_art))
print(paste("P-value for MAR vs SMART:", p_value_smart))
print(paste("P-value for STAR vs AR:", p_value_star))
