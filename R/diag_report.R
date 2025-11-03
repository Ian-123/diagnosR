#' Diagnostic report for an lm model (text summary)
#'
#' Runs common assumption checks and prints a concise, readable report.
#' Returns a list with all metrics invisibly.
#'
#' Checks included:
#' - Normality (Shapiro-Wilk on residuals)
#' - Homoskedasticity (Breusch-Pagan)
#' - Independence (Durbin-Watson)
#' - Outliers (|studentized residual| > 3)
#' - Influence (Cook's distance > 4/n)
#'
#' @param model Fitted \code{lm} object.
#' @return Invisibly returns a list of test results and flags.
#' @export
#' @importFrom stats residuals rstudent cooks.distance
diag_report <- function(model) {
  stopifnot(inherits(model, "lm"))

  # --- core quantities
  res   <- stats::residuals(model)
  rstu  <- stats::rstudent(model)
  cook  <- stats::cooks.distance(model)
  n     <- length(res)

  # --- tests (require lmtest)
  bp <- tryCatch(lmtest::bptest(model), error = function(e) NULL)
  dw <- tryCatch(lmtest::dwtest(model),  error = function(e) NULL)

  # Shapiro-Wilk (skip if n > 5000 to avoid spurious power)
  shap <- if (n <= 5000) {
    tryCatch(stats::shapiro.test(res), error = function(e) NULL)
  } else NULL

  # Outliers / Influence thresholds
  outlier_idx   <- which(abs(rstu) > 3)
  cooks_thresh  <- 4 / n
  infl_idx      <- which(cook > cooks_thresh)

  # Flags (TRUE = looks OK)
  ok_normal   <- is.null(shap) || shap$p.value > 0.05
  ok_bp       <- is.null(bp)   || bp$p.value   > 0.05
  ok_dw       <- is.null(dw)   || (dw$p.value  > 0.05)  # H0: no autocorr
  ok_outliers <- length(outlier_idx) == 0
  ok_cooks    <- length(infl_idx)    == 0

  cat("\n================  diagnosR :: diag_report  ================\n")
  cat("Model: ", deparse(formula(model)), "\n")
  cat("n =", n, "  p =", length(coef(model)), "\n\n")

  pflag <- function(ok) if (isTRUE(ok)) "PASS" else "ATTN"

  cat(sprintf("[%-4s] Normality (Shapiro): %s\n",
              pflag(ok_normal),
              if (is.null(shap)) "skipped (n > 5000)" else
                sprintf("W=%.3f, p=%.3f", shap$statistic, shap$p.value)))

  cat(sprintf("[%-4s] Homoskedasticity (Breusch–Pagan): %s\n",
              pflag(ok_bp),
              if (is.null(bp)) "unavailable (lmtest missing?)" else
                sprintf("BP=%.2f, p=%.3f", unname(bp$statistic), bp$p.value)))

  cat(sprintf("[%-4s] Independence (Durbin–Watson): %s\n",
              pflag(ok_dw),
              if (is.null(dw)) "unavailable (lmtest missing?)" else
                sprintf("DW=%.3f, p=%.3f", unname(dw$statistic), dw$p.value)))

  cat(sprintf("[%-4s] Outliers (|rstudent| > 3): %s\n",
              pflag(ok_outliers),
              if (ok_outliers) "none" else paste(length(outlier_idx), "obs")))

  cat(sprintf("\n[%-4s] Influence (Cook's D > 4/n): %s\n",
              pflag(ok_cooks),
              if (ok_cooks) "none" else paste(length(infl_idx), "obs")))
  if (!ok_cooks) {
    top <- head(order(cook, decreasing = TRUE), 5)
    cat("       Top Cook's D idx:", paste(top, collapse = ", "),
        "  values:", paste(round(cook[top], 3), collapse = ", "), "\n")
  }
  cat("============================================================\n\n")

  invisible(list(
    shapiro = shap, breusch_pagan = bp, durbin_watson = dw,
    outlier_index = outlier_idx, cooks_index = infl_idx,
    cooks_threshold = cooks_thresh
  ))
}
