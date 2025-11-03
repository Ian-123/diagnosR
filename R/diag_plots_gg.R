#' ggplot2 diagnostic panel for an lm model
#'
#' Produces a 2x2 ggplot panel:
#'  - Residuals vs Fitted (with smooth)
#'  - Normal Q-Q
#'  - Scale-Location (sqrt(|std resid|) vs Fitted)
#'  - Cook's Distance bars
#'
#' @param model Fitted \code{lm} object.
#' @return A grid object (invisibly). Also prints the panel.
#' @export
#' @import ggplot2
#' @importFrom stats rstandard residuals fitted cooks.distance
diag_plots_gg <- function(model) {
  stopifnot(inherits(model, "lm"))

  res   <- stats::residuals(model)
  fit   <- stats::fitted(model)
  rstd  <- stats::rstandard(model)
  cook  <- stats::cooks.distance(model)
  idx   <- seq_along(res)

  df <- data.frame(fit = fit, res = res, rstd = rstd, cook = cook, idx = idx)

  p1 <- ggplot(df, aes(fit, res)) +
    geom_point() + geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = "Residuals vs Fitted", x = "Fitted", y = "Residuals")

  p2 <- ggplot(df, aes(sample = res)) +
    stat_qq() + stat_qq_line() +
    labs(title = "Normal Q-Q", x = "Theoretical", y = "Sample")

  p3 <- ggplot(df, aes(fit, sqrt(abs(rstd)))) +
    geom_point() + geom_smooth(se = FALSE, method = "loess", formula = y ~ x) +
    labs(title = "Scale-Location", x = "Fitted", y = "sqrt(|Std Resid|)")

  thr <- 4 / nrow(df)
  p4 <- ggplot(df, aes(idx, cook)) +
    geom_col() +
    geom_hline(yintercept = thr, linetype = 2) +
    labs(title = "Cook's Distance", x = "Obs", y = "Cook's D")

  gridExtra::grid.arrange(p1, p2, p3, p4, nrow = 2)
  invisible(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}
