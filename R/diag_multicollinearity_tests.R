#' Multicollinearity diagnostics: VIF/GVIF and Condition Number

#'

#' Uses `car::vif` (handles factors via GVIF) and `kappa` on X'X.

#' @param model Fitted `lm` object.

#' @param vif_threshold Flag threshold for adjusted-VIF (default 10).

#' @param kappa_threshold Flag threshold for condition number (default 30).

#' @return A list with $vif_table and $summary (invisibly). Also prints a compact summary.

#' @importFrom stats var

#' @importFrom car vif

#' @export

# R/diag_multicollinearity_tests.R

diag_multicollinearity_tests <- function(model,

                                         vif_threshold = 10,

                                         kappa_threshold = 30) {


  # car::vif() with type='predictor' works with interactions

  v <- tryCatch(

    suppressMessages(vif(model, type = "predictor")),

    error = function(e) { message("VIF failed: ", e$message); NULL }

  )




  vif_table <- NULL

  if (!is.null(v)) {

    if (is.matrix(v)) {

      # GVIF table

      vif_table <- data.frame(

        variable = rownames(v),

        Df       = if ("Df" %in% colnames(v)) v[, "Df"] else 1L,

        GVIF     = if ("GVIF" %in% colnames(v)) v[, "GVIF"] else NA_real_,

        adj_VIF  = if ("GVIF^(1/(2*Df))" %in% colnames(v)) v[, "GVIF^(1/(2*Df))"] else as.numeric(v[,1]),

        row.names = NULL

      )

    } else if (is.numeric(v)) {

      # type='predictor' often returns a named numeric vector of adjusted VIFs

      vif_table <- data.frame(

        variable = names(v),

        Df       = 1L,

        GVIF     = NA_real_,

        adj_VIF  = as.numeric(v),

        row.names = NULL

      )

    }

  }



  # Condition number (exclude intercept)

  X <- tryCatch(stats::model.matrix(model), error = function(e) NULL)

  if (!is.null(X) && "(Intercept)" %in% colnames(X)) {

    X <- X[, setdiff(colnames(X), "(Intercept)"), drop = FALSE]

  }

  kappa_val <- if (!is.null(X) && ncol(X) > 0) as.numeric(kappa(X)) else NA_real_



  max_adj_vif <- if (!is.null(vif_table)) max(vif_table$adj_VIF, na.rm = TRUE) else NA_real_

  bad_vif_vars <- if (!is.null(vif_table)) paste(vif_table$variable[vif_table$adj_VIF > vif_threshold], collapse = ", ") else "None"



  cat("\n=== Multicollinearity diagnostics ===\n")

  if (!is.null(vif_table)) {

    print(vif_table, row.names = FALSE)

    cat(sprintf("\nMax adjusted VIF: %.3f  (threshold  %.1f)\n", max_adj_vif, vif_threshold))

    cat("Over threshold: ", if (bad_vif_vars == "") "None" else bad_vif_vars, "\n", sep = "")

  } else {

    cat("VIF unavailable (see message above).\n")

  }

  cat(sprintf("Condition number (kappa): %s  (flag if  %.0f)\n",

              ifelse(is.na(kappa_val), "NA", sprintf("%.3f", kappa_val)),

              kappa_threshold))

  cat("=====================================\n\n")



  invisible(list(vif_table = vif_table, max_adj_vif = max_adj_vif,

                 kappa = kappa_val,

                 bad_vars = if (bad_vif_vars %in% c("", "None")) character(0) else strsplit(bad_vif_vars, ",\\s*")[[1]]))

}


