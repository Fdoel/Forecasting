


#use inflation_df_monthly

first_try <- marx(inflation_df_monthly$inflation, NULL, p_max = 10, sig_level = 0.1, p_C=0, p_NC = 6)
second_try <- marx(inflation_df_monthly$inflation, NULL, p_max = 10, sig_level = 0.05, p_C = 3, p_NC = 3)

ff_proberen <- marx(inflation_df_monthly$inflation, NULL, p_max = 15, sig_level = 0.05, p_C=3)
