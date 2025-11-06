#' Autocorrelation tests for linear model residuals

#'

#' Runs DurbinWatson (lag-1) and BreuschGodfrey tests for given lags.

#' @param model Fitted `lm` object.

#' @param lags Integer vector of BG lags to test (default 1:4).

#' @return A data.frame with test names, lag, statistic and p-value (invisibly).

#' @importFrom stats Box.test residuals

#' @importFrom lmtest dwtest bgtest bptest

#' @export

# R/diag_autocorr_tests.R

diag_autocorr_tests <- function(model, lags = 1:4) {

  stopifnot(inherits(model, "lm"))

  res <- stats::residuals(model)

  res <- as.numeric(res)

  res <- res[is.finite(res)]

  if (length(res) < 5) stop("Not enough residuals to test autocorrelation.")



  have_lmtest <- requireNamespace("lmtest", quietly = TRUE)



  # --- DurbinWatson ---

  dw_stat <- sum(diff(res)^2) / sum(res^2)     # always available

  dw_p <- NA_real_

  if (have_lmtest) {

    dw_p <- tryCatch(lmtest::dwtest(model)$p.value, error = function(e) NA_real_)

  }

  # Report DW even if p-value is NA

  dw_row <- data.frame(test = "DurbinWatson", lag = 1L,

                       statistic = dw_stat, p.value = dw_p)



  # --- BreuschGodfrey; fallback to LjungBox if lmtest fails ---

  rows <- list(dw_row)

  for (L in lags) {

    if (have_lmtest) {

      bg <- tryCatch(lmtest::bgtest(model, order = L),

                     error = function(e) NULL)

      if (!is.null(bg)) {

        rows[[length(rows)+1]] <- data.frame(

          test = "BreuschGodfrey", lag = L,

          statistic = unname(bg$statistic), p.value = bg$p.value

        )

        next

      }

    }

    # Fallback: LjungBox on residuals

    lb <- stats::Box.test(res, lag = L, type = "Ljung-Box")

    rows[[length(rows)+1]] <- data.frame(

      test = "LjungBox (fallback)", lag = L,

      statistic = unname(lb$statistic), p.value = lb$p.value

    )

  }



  out <- do.call(rbind, rows)

  print(out, row.names = FALSE)

  invisible(out)

}


