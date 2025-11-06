#' Automatic diagnostic plots in base R style (multi-page) + leverage + labels
#'
#' @param model Fitted lm object.
#' @param data Optional data frame used to fit model.
#' @param per_page Plots per page after the standard 3x3 panel.
#' @param label_influential logical; if TRUE, label influential points on the
#'   Residuals vs Leverage plot. Default TRUE.
#' @param label_which one of c("both","cook","leverage"); which criterion to use
#'   for labeling. Default "both".
#' @param label_top integer or NULL; if set, label only the top-k by Cook's D
#'   among those exceeding the chosen thresholds (helps avoid clutter). Default NULL = label all exceeding.
#' @param cook_cut Cook's D cutoff; default 4/n (common rule-of-thumb).
#' @param lev_cut Leverage cutoff; default 2p/n.
#' @param labels optional character vector of point labels (length n). If NULL,
#'   uses rownames(data/model.frame) or indices.
#' @importFrom graphics plot.new title mtext abline boxplot par lines points text
#' @importFrom stats model.frame residuals fitted rstudent cooks.distance hatvalues coef na.omit qqnorm qqline rnorm
#' @importFrom utils head tail
#' @importFrom grDevices dev.new
#' @export
diag_plots <- function(model, data = NULL, per_page = 4,
                       label_influential = TRUE,
                       label_which = c("both","cook","leverage"),
                       label_top = NULL,
                       cook_cut = NULL,
                       lev_cut  = NULL,
                       labels = NULL) {
  stopifnot(inherits(model, "lm"))
  label_which <- match.arg(label_which)

  mf   <- if (is.null(data)) model.frame(model) else data
  ynam <- names(mf)[1]
  res   <- residuals(model)
  fit   <- fitted(model)
  rstud <- tryCatch(rstudent(model), error = function(e) rep(NA_real_, length(res)))
  cook  <- tryCatch(cooks.distance(model), error = function(e) rep(NA_real_, length(res)))
  hatv  <- tryCatch(hatvalues(model), error = function(e) rep(NA_real_, length(res)))
  p     <- length(coef(model))
  n     <- length(res)
  idx   <- seq_len(n)

  cook_thr <- if (is.null(cook_cut)) 4 / max(1, n) else cook_cut
  lev_thr  <- if (is.null(lev_cut))  2 * p / max(1, n) else lev_cut

  # default labels for points
  lab <- labels
  if (is.null(lab)) {
    rn <- rownames(mf)
    lab <- if (!is.null(rn)) rn else as.character(idx)
  }

  # helpers -------------------------------------------------------------
  .safe_scatter <- function(x, y, main, xlab, ylab) {
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) {
      plot.new(); title(main); mtext("No finite data to plot", side = 3, line = -1, col = "red"); return(invisible(FALSE))
    }
    plot(x[ok], y[ok], main = main, xlab = xlab, ylab = ylab)
    abline(h = 0, col = "red")
    invisible(TRUE)
  }
  .safe_boxplot <- function(y, f, main, xlab, ylab) {
    ok <- is.finite(y) & !is.na(f)
    if (!any(ok)) { plot.new(); title(main); mtext("No finite data to plot", side = 3, line = -1, col = "red"); return(invisible(FALSE)) }
    boxplot(y[ok] ~ droplevels(f[ok]), main = main, xlab = xlab, ylab = ylab)
    abline(h = 0, col = "red")
    invisible(TRUE)
  }
  cook_contours <- function(h, p, D) {
    r <- sqrt(pmax(0, D * p * (1 - h) / pmax(h, 1e-8)))
    list(upper = r, lower = -r)
  }

  # FIRST PAGE ----------------------------------------------------------
  oldpar <- par(mfrow = c(3, 3))
  on.exit(par(oldpar), add = TRUE)

  # 1) Residuals vs Response (coerce y if needed, but prefer numeric)
  yvals <- mf[[1]]
  if (!is.numeric(yvals)) {
    suppressWarnings(y_plot <- as.numeric(yvals))
  } else y_plot <- yvals
  .safe_scatter(y_plot, res, "Residuals vs Response", xlab = ynam, ylab = "Residuals")

  # 23) Residuals vs first two predictors
  x_names <- setdiff(names(mf), ynam)
  x_show <- head(x_names, 2)
  for (xn in x_show) {
    x <- mf[[xn]]
    if (is.factor(x)) {
      .safe_boxplot(res, x, main = paste("Residuals vs", xn), xlab = xn, ylab = "Residuals")
    } else {
      .safe_scatter(x, res, main = paste("Residuals vs", xn), xlab = xn, ylab = "Residuals")
    }
  }
  while (length(x_show) < 2) { plot.new(); length(x_show) <- 2 }

  # 4) Outliers (studentized)
  .safe_scatter(fit, rstud, "Outliers (rstudent)", xlab = "Fitted", ylab = "Studentized Residuals")

  # 5) Homoskedasticity
  .safe_scatter(fit, res, "Homoskedasticity (Res vs Fitted)", xlab = "Fitted", ylab = "Residuals")

  # 6) Normal Q-Q
  qqnorm(na.omit(res), main = "Normal Q-Q"); qqline(na.omit(res), col = "red")

  # 7) Independence
  .safe_scatter(idx, res, "Independence (Res vs Order)", xlab = "Observation Order", ylab = "Residuals")

  # 8) Cook's Distance
  okc <- is.finite(cook)
  if (any(okc)) {
    plot(idx[okc], cook[okc], type = "h", main = "Cook's Distance", xlab = "Observation", ylab = "Cook's D")
    abline(h = cook_thr, col = "red")
  } else { plot.new(); title("Cook's Distance"); mtext("No finite data to plot", side = 3, line = -1, col = "red") }

  # 9) Residuals vs Leverage
  okrl <- is.finite(hatv) & is.finite(rstud)
  if (any(okrl)) {
    plot(hatv[okrl], rstud[okrl], xlab = "Leverage (hat values)", ylab = "Studentized Residuals",
         main = "Residuals vs Leverage")
    abline(h = 0, col = "red")
    abline(v = lev_thr, lty = 2, col = "blue")

    # contours if we have finite leverage range
    hmin <- min(hatv[okrl]); hmax <- max(hatv[okrl])
    if (is.finite(hmin) && is.finite(hmax) && hmax > 0) {
      hseq <- seq(max(1e-6, hmin * 0.98), min(0.999, hmax * 1.02), length.out = 200)
      for (D in c(0.5, 1)) {
        cc <- cook_contours(hseq, p, D)
        lines(hseq, cc$upper, lty = 3)
        lines(hseq, cc$lower, lty = 3)
        text(tail(hseq, 1), tail(cc$upper, 1), labels = paste0("Cook's D=", D), pos = 4, cex = 0.7)
      }
    }

    # optional labels
    if (isTRUE(label_influential)) {
      flag_cook <- cook > cook_thr
      flag_lev  <- hatv > lev_thr
      flag <- switch(label_which,
                     "cook"     = flag_cook,
                     "leverage" = flag_lev,
                     "both"     = flag_cook | flag_lev
      )
      ids <- which(flag & okrl)
      if (length(ids)) {
        if (!is.null(label_top) && is.finite(label_top) && label_top > 0 && length(ids) > label_top) {
          ids <- ids[order(cook[ids], decreasing = TRUE)][seq_len(label_top)]
        }
        # light jitter to reduce overlap
        xr <- diff(range(hatv[okrl])); yr <- diff(range(rstud[okrl]))
        jx <- hatv[ids] + rnorm(length(ids), sd = 0.01 * ifelse(is.finite(xr), xr, 1))
        jy <- rstud[ids] + rnorm(length(ids), sd = 0.01 * ifelse(is.finite(yr), yr, 1))
        text(jx, jy, labels = lab[ids], cex = 0.75, pos = 3, xpd = NA)
        points(hatv[ids], rstud[ids], pch = 21, bg = "yellow", cex = 1.1)
      }
    }
  } else { plot.new(); title("Residuals vs Leverage"); mtext("No finite data to plot", side = 3, line = -1, col = "red") }

  # extra pages ---------------------------------------------------------
  x_left <- setdiff(x_names, x_show)
  if (length(x_left) > 0) {
    groups <- split(x_left, ceiling(seq_along(x_left) / per_page))
    for (g in groups) {
      dev.new()
      par(mfrow = c(2, 2))
      for (xn in g) {
        x <- mf[[xn]]
        if (is.factor(x)) {
          .safe_boxplot(res, x, main = paste("Residuals vs", xn), xlab = xn, ylab = "Residuals")
        } else {
          .safe_scatter(x, res, main = paste("Residuals vs", xn), xlab = xn, ylab = "Residuals")
        }
      }
    }
  }

  invisible(list(
    residuals = res, fitted = fit, rstudent = rstud,
    cooks = cook, leverage = hatv, cook_thr = cook_thr, lev_thr = lev_thr
  ))
}

