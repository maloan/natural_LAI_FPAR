# Statistical helpers for trend analysis scripts

compute_ols_full <- function(d, y_col = "value", x_col = "year") {
  d <- d |> dplyr::filter(is.finite(.data[[y_col]]), is.finite(.data[[x_col]]))

  if (nrow(d) < 3) {
    return(tibble::tibble(
      n = nrow(d),
      slope = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      p_value = NA_real_,
      r2 = NA_real_,
      sig = NA_character_
    ))
  }

  m <- stats::lm(
    formula(paste(y_col, "~", x_col)),
    data = d
  )

  sm <- summary(m)
  p_val <- sm$coefficients[2, 4]

  compute_trend_ci(m, sm, p_val)
}

compute_trend_ci <- function(lm_fit, summary_obj = NULL, p_value = NULL) {
  if (is.null(summary_obj)) {
    summary_obj <- summary(lm_fit)
  }

  slope <- coef(lm_fit)[2]
  se <- summary_obj$coefficients[2, 2]
  tcrit <- stats::qt(0.975, df = lm_fit$df.residual)

  if (is.null(p_value)) {
    p_value <- summary_obj$coefficients[2, 4]
  }

  tibble::tibble(
    n = nrow(lm_fit$model),
    slope = slope,
    ci_lower = slope - tcrit * se,
    ci_upper = slope + tcrit * se,
    p_value = p_value,
    r2 = summary_obj$r.squared,
    sig = ifelse(p_value < 0.05, "*", "")
  )
}

# Predict OLS trend line with CIs for plotting
predict_trend_line <- function(d, y_col = "value", x_col = "year") {
  d <- d |> dplyr::filter(is.finite(.data[[y_col]]), is.finite(.data[[x_col]]))

  if (nrow(d) < 3) {
    return(tibble::tibble(
      !!x_col := d[[x_col]],
      fit = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    ))
  }

  m <- stats::lm(
    formula(paste(y_col, "~", x_col)),
    data = d
  )

  pr <- stats::predict(m, newdata = d, se.fit = TRUE)
  tcrit <- stats::qt(0.975, df = m$df.residual)

  tibble::tibble(
    !!x_col := d[[x_col]],
    fit = pr$fit,
    ci_lower = pr$fit - tcrit * pr$se.fit,
    ci_upper = pr$fit + tcrit * pr$se.fit
  )
}
