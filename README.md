# DiagnosR
**Quick Diagnostic Panel for Linear Models**
`diagnosR` provides a one-call solution for generating diagnostic plots and tests for linear regression models (`lm`). It helps identify common issues such as multicollinearity, autocorrelation, influential observations, and more.
# a install.packages("remotes") # if needed, that is, check if you have package already installed and install if not
remotes::install_github("Ian-123/diagnosR")
# or install.packages(???devtools???) # installs devtools
# devtools::install_github(???Ian-123/diagnosR???) #local install
#load library(diagnosR)
## b( Features
- **Diagnostic Plots**: Residuals vs fitted, Cook's distance, leverage plots, Q-Q plots.
- **Multicollinearity Checks**: VIF/GVIF and condition number.
- **Autocorrelation Tests**: Durbin???Watson and Breusch???Godfrey (with Ljung???Box fallback).
- **Compact Scorecard**: One-page summary of key diagnostics.Autocorrelation Tests.

**Example**
m <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris) 
- diag_plots(m) # base R diagnostics plots for linearity, normality, homoskedasticity, influence,independence
- diag_scorecard(m) # returns a report of basic and added diagnostics in one page
- diag_multicollinearity_tests(m) # Multicollinearity diagnostics
- diag_autocorr_tests(m) # Autocorrelation tests
