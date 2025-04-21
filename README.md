## Licensing and Credits

This project is based in part on the R package MARX originally developed by Sean Telg [aut, cre, cph], Alain Hecq [ctb], Lenard Lieb [ctb] and licensed under the [GNU General Public License v2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).

Modifications have been made to suit the needs of this project. All derived work remains under GPL-2 as per the license terms.

Please see the `LICENSE` file for full licensing details.

## Explanation of the repository

# Marx_functions.R: 

The original functions of the MARX package, with some minor tweaks made, e.g. the line breaks before some else statements were removed to prevent errors with newer versions of R

# MART.R:

The estimation, and forecast function of the MART and SMART models. This contains the main usability of MART and SMART for other use than to validate research results.

# GridSearch.R 

The code used grid searching the models used in the paper. The code needs to be adjusted for MART and SMART models according to the variables needed by the model.

# Descriptive statistics (..).R

Contains the preprocessing and some descriptive statistics of the various datasets used in the paper. A seperate version exists using Seasonally adjusted data, one using quarterly data and one with the main non-seasonally adjusted dataset used in the paper.

# Pseudo MART and SMART

These files contains the functionality of estimating pseudo causal MART and SMART models, as well as the code for the tests on the residuals of these models

# Inflation_df_monthly.RData

The main preprocessed dataset used in the paper

# Finalmodels.R 

Contains the estimation of most of the models used in the results section of the paper

# RMSFE.R

An attempt to centralise all previously ran root mean squared forecast errors. Also plots the forecasts.

# Basic 

Initial file used for running a grid search and residual tests for the AR model.

# Simulation/

The code used to simulate the various models

# Forecasts/

The main scripts used to run the forecasts. Due to the large amount of forecasts used in the paper, some scripts may need changing the parameters of the model to match that in the papaper.

# Robustness/
Scripts used for robustness that are not included elsewhere

# man/

The R markdown files of the original MARX package.

# Grid searches/
The results of the grid searches.

# Forecasting results/

.RData files of forecasts

# Simulation results/
.RData files containing results from the simulation

# FRED.CSV
The FRED MD dataset.

# CPOI US labour dataset
The raw unseasonaly adjusted CPI dataset.

