#Ian Njuguna 2025
#' Diagnostic plot panel for linear models
#'
#' Produces a 3x3 base-R panel that replicates the common
#' checks you listed: linearity (vs Y and vs chosen Xs),
#' outliers (studentized residuals), homoskedasticity,
#' normal Q-Q, and independence (residuals over order).
#'
#' @param model A fitted lm object.
#' @param data Optional data.frame used to fit the model
#'   (if not supplied, tries to extract via model.frame()).
#' @param y Name of response variable (character). If NULL,
#'   inferred from the model.
#' @param x_vars Character vector of up to 2 predictor names
#'   to show residual-vs-X plots. If NULL, uses the first two
#'   non-intercept terms in the model.
#' @param order_index Optional numeric/ordering index for the
#'   independence plot. If NULL, uses the row order in data/model.
#' @param save_path Optional file path (e.g., "diag_panel.png").
#'   If supplied, the panel will be saved (PNG) and device closed.
#'
#' @return Invisibly returns a list with components used for plotting.
#' @export
diag_plots <- function(model,
                       data = NULL,
                       y = NULL,
                       x_vars = NULL,
                       order_index = NULL,
                       save_path = NULL) {

  stopifnot(inherits(model, "lm"))

  # Pull model frame & pieces
  mf <- tryCatch(model.frame(model),
                 error = function(...) NULL)

  if (is.null(data)) {
    data <- mf
  }
  if (is.null(data)) {
    stop("Could not infer `data`. Please pass the original data used in the model.")
  }

  # Response name
  if (is.null(y)) {
    tt <- terms(model)
    y <- deparse(attr(tt, "variables")[[2]])
  }
  if (!y %in% names(data)) {
    stop("`y` not found in data: ", y)
  }

  # Pick 2 Xs (for your specific residual-vs-X plots)
  # Prefer non-intercept, unique base terms
  term_labels <- attr(terms(model), "term.labels")
  # Strip interaction expansions to base names when possible
  base_terms <- unique(gsub(":.*|\\*.*", "", term_labels))
  base_terms <- base_terms[base_terms != "(Intercept)"]

  if (is.null(x_vars)) {
    x_vars <- base_terms[seq_len(min(2L, length(base_terms)))]
  }
  x_vars <- x_vars[x_vars %in% names(data)]
  if (length(x_vars) == 0) x_vars <- character(0)

  # Residuals, fitted, studentized residuals
  res  <- residuals(model)
  fit  <- fitted(model)
  rstu <- rstudent(model)

  # Order index for independence plot
  if (is.null(order_index)) {
    order_index <- seq_along(res)
  }
  if (length(order_index) != length(res)) {
    stop("`order_index` must have same length as residuals.")
  }

  # Optional save to file
  if (!is.null(save_path)) {
    grDevices::png(filename = save_path, width = 1400, height = 1400, res = 150)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

  ## 1) Residuals vs Y (linearity check vs response)
  plot(data[[y]], res,
       xlab = y, ylab = "Residuals",
       main = "Residuals vs Response")
  abline(h = 0, col = "gray")

  ## 2) Residuals vs X1
  if (length(x_vars) >= 1) {
    x1 <- data[[x_vars[1]]]
    if (is.factor(x1)) {
      boxplot(res ~ x1, xlab = x_vars[1], ylab = "Residuals",
              main = "Residuals vs X1 (factor)")
      abline(h = 0, col = "gray")
    } else {
      plot(x1, res, xlab = x_vars[1], ylab = "Residuals",
           main = "Residuals vs X1")
      abline(h = 0, col = "gray")
    }
  } else {
    plot.new(); title("Residuals vs X1 (none)")
  }

  ## 3) Residuals vs X2
  if (length(x_vars) >= 2) {
    x2 <- data[[x_vars[2]]]
    if (is.factor(x2)) {
      boxplot(res ~ x2, xlab = x_vars[2], ylab = "Residuals",
              main = "Residuals vs X2 (factor)")
      abline(h = 0, col = "gray")
    } else {
      plot(x2, res, xlab = x_vars[2], ylab = "Residuals",
           main = "Residuals vs X2")
      abline(h = 0, col = "gray")
    }
  } else {
    plot.new(); title("Residuals vs X2 (none)")
  }

  ## 4) Outliers: Studentized residuals vs Fitted
  plot(fit, rstu,
       xlab = "Fitted", ylab = "Studentized Residuals",
       main = "Outliers (rstudent vs Fitted)")
  abline(h = 0, col = "gray")

  ## 5) Homoskedasticity: Residuals vs Fitted
  plot(fit, res,
       xlab = "Fitted", ylab = "Residuals",
       main = "Homoskedasticity (Res vs Fitted)")
  abline(h = 0, col = "gray")

  ## 6) Normality: Q-Q plot
  qqnorm(res, main = "Normal Q-Q (Residuals)")
  qqline(res, col = "red")

  ## 7) Independence: residuals over order
  plot(order_index, res, type = "b",
       xlab = "Observation Order", ylab = "Residuals",
       main = "Independence (Res vs Order)")
  abline(h = 0, col = "red")

  ## 8–9) Free slots: leverage/influence & Cook's distance
  # Leverage vs residuals
  hatv <- hatvalues(model)
  plot(hatv, abs(rstu),
       xlab = "Leverage (hat values)", ylab = "|Studentized Residual|",
       main = "Influence (Leverage vs |rstudent|)")
  # Cook's distance
  cd <- cooks.distance(model)
  plot(cd, type = "h", ylab = "Cook's Distance",
       main = "Cook's Distance (Influence)")
  abline(h = 4/length(cd), col = "gray", lty = 2)

  invisible(list(
    y = y,
    x_vars = x_vars,
    fitted = fit,
    residuals = res,
    rstudent = rstu,
    leverage = hatv,
    cooks_distance = cd,
    order_index = order_index
  ))
}
