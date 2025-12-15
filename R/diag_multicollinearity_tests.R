#' Multicollinearity diagnostics: VIF/GVIF and Condition Number

#'

#' Uses `car::vif` (handles factors via GVIF) and `kappa` on X'X.

#'

#' @param model Fitted `lm` object.

#' @param vif_threshold Flag threshold for adjusted-VIF (default 10).

#' @param kappa_threshold Flag threshold for condition number (default 30).

#' @param print Logical; if TRUE, prints the GVIF table and summary. Default TRUE.

#' @return list(vif_table, max_adj_vif, kappa, flags)

#' @importFrom car vif

#' @export

diag_multicollinearity_tests <- function(model,

                                         vif_threshold = 10,

                                         kappa_threshold = 30,

                                         print = TRUE) {

  stopifnot(inherits(model, "lm"))

  vif_table <- NULL



  # 1) Try car::vif matrix (GVIF, Df, GVIF^(1/(2*Df)))

v_mat <- try(suppressMessages(car::vif(model)), silent = TRUE)

  if (!inherits(v_mat, "try-error") && is.matrix(v_mat)) {

    vif_table <- data.frame(

      variable  = rownames(v_mat),

      GVIF      = if ("GVIF" %in% colnames(v_mat)) v_mat[, "GVIF"] else NA_real_,

      Df        = if ("Df"   %in% colnames(v_mat)) v_mat[, "Df"]   else 1L,

      adj_VIF   = if ("GVIF^(1/(2*Df))" %in% colnames(v_mat)) {

        v_mat[, "GVIF^(1/(2*Df))"]

      } else v_mat[, "GVIF"]^(1/(2 * v_mat[, "Df"])),

      row.names = NULL,

      check.names = FALSE

    )

  } else {

# 2) Fallback: adjusted-only vector
v_vec <- try(suppressMessages(car::vif(model, type = "predictor")), silent = TRUE)

if (!inherits(v_vec, "try-error") && is.numeric(v_vec)) {

  # car::vif() may return an unnamed numeric (or a numeric matrix)
  v_num <- as.numeric(v_vec)
  n <- length(v_num)

  var_names <- names(v_vec)

  # If names are missing, try to recover from model matrix; else create placeholders
  if (is.null(var_names) || length(var_names) != n) {
    mm <- try(stats::model.matrix(model), silent = TRUE)
    if (!inherits(mm, "try-error")) {
      cn <- colnames(mm)
      cn <- setdiff(cn, "(Intercept)")
      if (length(cn) == n) {
        var_names <- cn
      }
    }
  }
  if (is.null(var_names) || length(var_names) != n) {
    var_names <- paste0("V", seq_len(n))
  }

  vif_table <- data.frame(
    variable  = var_names,
    GVIF      = rep(NA_real_, n),
    Df        = rep(1L, n),
    adj_VIF   = v_num,
    row.names = NULL,
    check.names = FALSE
  )
}
}

  # 3) Last-resort: compute VIF from model matrix if still NULL

  if (is.null(vif_table)) {

    X <- try(stats::model.matrix(model), silent = TRUE)

    if (!inherits(X, "try-error")) {

      if ("(Intercept)" %in% colnames(X)) {

        X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]

      }

      keep <- apply(X, 2, function(z) is.finite(stats::var(z)) && stats::var(z) > 0)

      if (any(keep)) {

        X <- X[, keep, drop = FALSE]

        nm <- colnames(X)

        vifv <- rep(NA_real_, ncol(X)); names(vifv) <- nm

        if (ncol(X) == 1L) {

          vifv[] <- 1

        } else {

          for (j in seq_along(nm)) {

            yj <- X[, j]; Z <- X[, -j, drop = FALSE]

            fit <- try(stats::lm.fit(Z, yj), silent = TRUE)

            if (!inherits(fit, "try-error")) {

              rss <- sum(fit$residuals^2); tss <- sum((yj - mean(yj))^2)

              R2 <- if (tss > 0) 1 - rss/tss else NA_real_

              if (!is.na(R2) && R2 < 1) vifv[j] <- 1/(1 - R2)

            }

          }

        }

        vif_table <- data.frame(

          variable = names(vifv),

          GVIF     = NA_real_,

          Df       = 1L,

          adj_VIF  = as.numeric(vifv),

          row.names = NULL,

          check.names = FALSE

        )

      }

    }

  }



  # kappa

  Xk <- try(stats::model.matrix(model), silent = TRUE)

  if (!inherits(Xk, "try-error") && "(Intercept)" %in% colnames(Xk))

    Xk <- Xk[, setdiff(colnames(Xk), "(Intercept)"), drop = FALSE]

  kappa_val <- if (!inherits(Xk, "try-error") && ncol(Xk) > 0) as.numeric(kappa(Xk)) else NA_real_



  max_adj_vif <- if (!is.null(vif_table) && any(!is.na(vif_table$adj_VIF)))

    max(vif_table$adj_VIF, na.rm = TRUE) else NA_real_



  flags <- list(

    vif_ok   = is.na(max_adj_vif) || max_adj_vif <= vif_threshold,

    kappa_ok = is.na(kappa_val)   || kappa_val   <  kappa_threshold

  )



  # ---- Print if requested ----

  if (isTRUE(print) && !is.null(vif_table) && nrow(vif_table) > 0) {

    cat("\nGVIF table (car::vif)\n")

    out <- vif_table[, c("variable", "GVIF", "Df", "adj_VIF")]

    names(out) <- c("Variable", "GVIF", "Df", "GVIF^(1/(2*Df))")

    print(out, row.names = FALSE)

    cat(sprintf("\nMax adj VIF: %s\n",

                ifelse(is.na(max_adj_vif), "NA", sprintf("%.3f", max_adj_vif))))

    cat(sprintf("Kappa: %s\n",

                ifelse(is.na(kappa_val), "NA", sprintf("%.3f", kappa_val))))

  }



  invisible(list(

    vif_table = vif_table,

    max_adj_vif = max_adj_vif,

    kappa = kappa_val,

    flags = flags

  ))

}




