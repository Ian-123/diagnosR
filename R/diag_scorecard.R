#' Internal: multicollinearity metrics for diag_scorecard

#' @keywords internal

#' @noRd

.mc_metrics_internal  <- function(model, vif_threshold = 10, kappa_threshold = 30) {

  X <- stats::model.matrix(model)

  if ("(Intercept)" %in% colnames(X)) {

    X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]

  }



  # drop zero-variance / all-NA columns

  keep <- apply(X, 2, function(z) is.finite(var(z)) && var(z) > 0)

  if (any(!keep)) X <- X[, keep, drop = FALSE]



  # VIF via 1/(1 - R2) fallback

  vif_vec <- rep(NA_real_, ncol(X))

  nm <- colnames(X)

  if (ncol(X) >= 2) {

    for (j in seq_along(nm)) {

      yj <- X[, j]

      Z  <- X[, -j, drop = FALSE]

      fit <- tryCatch(stats::lm.fit(Z, yj), error = function(e) NULL)

      if (is.null(fit)) next

      rss <- sum(fit$residuals^2)

      tss <- sum((yj - mean(yj))^2)

      R2  <- if (tss > 0) 1 - rss/tss else NA_real_

      if (!is.na(R2) && R2 < 1) vif_vec[j] <- 1 / (1 - R2)

    }

    names(vif_vec) <- nm

  } else if (ncol(X) == 1) {

    vif_vec[1] <- 1; names(vif_vec) <- nm

  }



  # condition number

  kappa_val <- tryCatch(kappa(X), error = function(e) NA_real_)



  list(

    vif_vector  = vif_vec,

    max_adj_vif = if (any(is.finite(vif_vec))) max(vif_vec, na.rm = TRUE) else NA_real_,

    kappa       = kappa_val,

    flags       = list(

      vif_ok   = if (any(is.finite(vif_vec))) max(vif_vec, na.rm = TRUE) <= vif_threshold else NA,

      kappa_ok = if (is.finite(kappa_val)) kappa_val < kappa_threshold else NA

    )

  )

}



# ---- main (exported) ----

#' One-page diagnostic scorecard for lm models (with influential rows)

#'

#' @param model Fitted \code{lm} object.

#' @param vif_threshold Flag if adjusted VIF > this (default 10).

#' @param kappa_threshold Flag if condition number >= this (default 30).

#' @param include_bp Logical; if TRUE, runs Breusch-Pagan test. Default FALSE.

#' @param include_reset Logical; if TRUE, runs Ramsey RESET test. Default FALSE.

#' @param lb_lags Integer vector of lags for Ljung-Box fallback (used inside autocorr).

#' @param export_csv Optional file path to save the scorecard table as CSV.

#' @param show_influence_rows Logical; if TRUE, prints row indices over the Cook's D cutoff.

#' @param max_rows_to_print Max influential row indices to print inline. Default 10.

#' @return Invisibly, a list with \code{table}, \code{influential_rows}, and cutoff.

#' @importFrom utils write.csv

#' @export

diag_scorecard <- function(model,

                           vif_threshold = 10,

                           kappa_threshold = 30,

                           include_bp = FALSE,

                           include_reset = FALSE,

                           lb_lags = 1:4,

                           export_csv = NULL,

                           show_influence_rows = TRUE,

                           max_rows_to_print = 10) {

  stopifnot(inherits(model, "lm"))



  fmt_p <- function(p) ifelse(is.na(p), "NA", sprintf("%.3f", p))

  flag  <- function(ok) if (is.na(ok)) "-" else if (ok) "PASS" else "ATTN"



  mf   <- stats::model.frame(model)

  y    <- names(mf)[1]

  res  <- stats::residuals(model)

  fit  <- stats::fitted(model)

  n    <- length(res)

  r2   <- summary(model)$r.squared

  aic  <- tryCatch(stats::AIC(model), error = function(e) NA_real_)

  fobj <- summary(model)$fstatistic

  ftxt <- if (!is.null(fobj)) {

    sprintf("F(%d,%d)=%.2f, p<%.001f", fobj[2], fobj[3], fobj[1],

            max(1e-3, stats::pf(fobj[1], fobj[2], fobj[3], lower.tail = FALSE)))

  } else "-"



  # Normality (skip for very large n)

  shap <- if (n <= 5000) tryCatch(stats::shapiro.test(res), error = function(e) NULL) else NULL

  shap_p <- if (is.null(shap)) NA_real_ else shap$p.value

  ok_norm <- is.na(shap_p) || shap_p > 0.05



  # Outliers / Influence

  rstud <- tryCatch(stats::rstudent(model), error = function(e) rep(NA_real_, n))

  n_out <- sum(abs(rstud) > 3, na.rm = TRUE)



  cook  <- tryCatch(stats::cooks.distance(model), error = function(e) rep(NA_real_, n))

  cook_thr <- 4 / n

  infl_idx <- which(is.finite(cook) & cook > cook_thr)

  n_inf <- length(infl_idx)

  ok_out <- n_out == 0

  ok_inf <- n_inf == 0



  infl_inline <- if (n_inf == 0) {

    "none"

  } else {

    if (n_inf <= max_rows_to_print) paste(infl_idx, collapse = ", ")

    else paste0(paste(infl_idx[1:max_rows_to_print], collapse = ", "), " ... +", n_inf - max_rows_to_print, " more")

  }

  # Multicollinearity (helper)

  mc <- diag_multicollinearity_tests(

    model,

    vif_threshold = vif_threshold,

    kappa_threshold = kappa_threshold

  )



  max_adj_vif <- if (!is.null(mc)) mc$max_adj_vif else NA_real_

  kappa_val   <- if (!is.null(mc)) mc$kappa       else NA_real_

  ok_vif      <- is.na(max_adj_vif) || max_adj_vif <= vif_threshold

  ok_kappa    <- is.na(kappa_val)   || kappa_val   <  kappa_threshold



  # Autocorrelation (helper: DW + BG/LB fallback)

  ac <- tryCatch(diag_autocorr_tests(model, lags = lb_lags), error = function(e) NULL)

  dw_stat <- if (!is.null(ac)) ac$statistic[ac$test == "Durbin-Watson"][1] else NA_real_

  ac_pvals <- if (!is.null(ac)) ac$p.value[ac$test != "Durbin-Watson"] else NA_real_

  ac_pmax  <- if (length(ac_pvals)) suppressWarnings(max(ac_pvals, na.rm = TRUE)) else NA_real_

  ok_ac    <- is.na(ac_pmax) || ac_pmax > 0.05



  # Optional: BP & RESET

  bp_p <- NA_real_; ok_bp <- NA

  if (isTRUE(include_bp) && requireNamespace("lmtest", quietly = TRUE)) {

    bp   <- tryCatch(lmtest::bptest(model), error = function(e) NULL)

    bp_p <- if (is.null(bp)) NA_real_ else bp$p.value

    ok_bp <- is.na(bp_p) || bp_p > 0.05

  }

  rs_p <- NA_real_; ok_rs <- NA

  if (isTRUE(include_reset) && requireNamespace("lmtest", quietly = TRUE)) {

    rs   <- tryCatch(lmtest::resettest(model, power = 2:3, type = "fitted"), error = function(e) NULL)

    rs_p <- if (is.null(rs)) NA_real_ else rs$p.value

    ok_rs <- is.na(rs_p) || rs_p > 0.05

  }



  # Scorecard table

  tab <- data.frame(

    Check   = c("Linearity (visual)",

                "Homoskedasticity (visual)",

                "Independence (visual)",

                "Normality",

                "Outliers",

                "Influence",

                "Autocorrelation",

                "Multicollinearity (VIF)",

                "Condition number",

                if (isTRUE(include_bp)) "Homoskedasticity (BP)" else NULL,

                if (isTRUE(include_reset)) "Functional form (RESET)" else NULL,

                "Overall fit"),

    Result  = c("Use Residuals vs Fitted plot",

                "Use Residuals vs Fitted plot",

                "Use Residuals vs Order plot",

                flag(ok_norm),

                flag(ok_out),

                flag(ok_inf),

                flag(ok_ac),

                flag(ok_vif),

                flag(ok_kappa),

                if (isTRUE(include_bp)) flag(ok_bp) else NULL,

                if (isTRUE(include_reset)) flag(ok_rs) else NULL,

                "-"),

    Details = c("Loess ~ flat; no systematic curve",

                "Variance looks constant (no funnel spread)",

                "No obvious runs/pattern; points oscillate around 0",

                if (is.na(shap_p)) "Shapiro skipped/NA" else paste0("Shapiro p=", fmt_p(shap_p)),

                paste0("|rstudent|>3: ", n_out),

                paste0("Cook's D >", sprintf('%.3f', cook_thr), " -> rows: ", infl_inline),

                paste0("DW~", ifelse(is.na(dw_stat), "NA", sprintf("%.2f", dw_stat)),

                       if (!is.na(ac_pmax)) paste0(", max BG/LB p=", fmt_p(ac_pmax)) else ""),

                if (is.na(max_adj_vif)) "-" else paste0("Max adj VIF=", sprintf("%.2f", max_adj_vif), " (<= ", vif_threshold, ")"),

                if (is.na(kappa_val)) "-" else paste0("kappa=", sprintf("%.2f", kappa_val), " (flag >= ", kappa_threshold, ")"),

                if (isTRUE(include_bp)) paste0("BP p=", fmt_p(bp_p)) else NULL,

                if (isTRUE(include_reset)) paste0("RESET p=", fmt_p(rs_p)) else NULL,

                paste0("R^2=", sprintf('%.3f', r2),

                       if (!is.na(aic)) paste0(", AIC=", sprintf("%.1f", aic)) else "",

                       if (!is.null(fobj)) paste0(", ", ftxt) else ""))

  )



  # Print

  cat("\n===============  diagnosR :: scorecard  ===============\n")

  cat("Model:  ", deparse(stats::formula(model)), "\n")

  cat("n = ", n, "\n\n", sep = "")

  print(tab, row.names = FALSE, right = FALSE)

  cat("=======================================================\n")



  # VIF table

  if (!is.null(mc) && !is.null(mc$vif_table) && nrow(mc$vif_table) > 0) {

    cat("\nGVIF table (car::vif)\n")

    out <- mc$vif_table[, c("variable","GVIF","Df","adj_VIF")]

    names(out) <- c("Variable","GVIF","Df","GVIF^(1/(2*Df))")

    print(out, row.names = FALSE)

  } else {

    cat("\nVIF not available.\n")

  }


  # Influential rows

  if (isTRUE(show_influence_rows) && n_inf > 0) {

    data_src <- tryCatch(model$model, error = function(e) NULL)

    cat("\nInfluential observations (Cook's D >", round(cook_thr, 3), ")\n", sep = "")

    if (!is.null(data_src)) {

      out <- data_src[infl_idx, , drop = FALSE]

      out$CooksD <- cook[infl_idx]

      out <- cbind(Row = infl_idx, out)

      print(out)

    } else {

      print(data.frame(Row = infl_idx, CooksD = cook[infl_idx]))

      cat("Original data not attached to model; showing indices only.\n")

    }

  }



  if (!is.null(export_csv)) utils::write.csv(tab, export_csv, row.names = FALSE)



  invisible(list(

    table = tab,

    influential_rows = infl_idx,

    influential_cutoff = cook_thr

  ))

}  # <-- single closing brace for diag_scorecard()




