# DiagnosR

**Quick Diagnostic Panel for Linear Models**

`diagnosR` provides a one-call solution for generating diagnostic plots and tests for linear regression models (`lm`).  
It helps identify common issues such as multicollinearity, autocorrelation, influential observations, and more.

---

# How to install the package

## Option 1: Using remotes

install.packages("remotes")   # if not already installed

library(remotes) # Load remotes

remotes::install_github("Ian-123/diagnosR") # Install  package from GitHub

library(diagnosR) # Load package

## Option 2: Using devtool

install.packages("devtools") # Install devtools (only if not installed)

library(devtools) # Load devtools

install_github("Ian-123/diagnosR") # Install your package from GitHub

library(diagnosR)  # Load package

---

#  Features
- **Diagnostic Plots**: Residuals vs fitted, Cook's distance, leverage plots, Q-Q plots.
- **Multicollinearity Checks**: VIF/GVIF and condition number.
- **Autocorrelation Tests**: Durbin Watson and Breusch Godfrey (with Ljung Box fallback).
- **Compact Scorecard**: One-page summary of key diagnostics.Autocorrelation Tests.

---

# Example

m <- lm(Sepal.Length ~ Sepal.Width + Petal.Length, data = iris)

diag_plots(m) **Base R diagnostic plots for linearity, normality, homoskedasticity, and influence**

diag_scorecard(m) **Compact one-page scorecard of key diagnostics**

diag_multicollinearity_tests(m) **Multicollinearity diagnostics**

diag_autocorr_tests(m) **Autocorrelation tests**
