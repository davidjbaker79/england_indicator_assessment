#' Build and plot a Pareto front for species targets
#'
#' Convenience function for summarising optimisation results across a sequence
#' of species targets and constructing a Pareto-front plot relating the spatial
#' footprint of the solution to the proportion of species meeting the abundance
#' target.
#'
#' For each optimisation result, the function extracts the realised number of
#' species meeting the target (\code{N10}), the associated planning footprint,
#' and other state statistics. Where multiple optimisation runs produce the
#' same realised \code{N10}, only the solution with the smallest spatial
#' footprint is retained.
#'
#' A monotonic Pareto front is then constructed and fitted using a
#' shape-constrained additive model. The derivative of the fitted relationship
#' is used to classify sections of the front according to the magnitude of
#' marginal gains in species reaching the abundance target.
#'
#' @param pareto_front_on_K A list of optimisation results generated for
#'   progressively increasing species targets. Each element is expected to
#'   contain at least \code{K}, \code{best_n_cells}, and \code{best_state};
#'   \code{best_state} must contain \code{x} and \code{beta}. The object stored
#'   in the second element of each result must be compatible with
#'   \code{compute_state_fast()}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{[[1]]}}{A \code{ggplot} object showing the Pareto front between
#'     planning-domain footprint and the percentage of species reaching the
#'     abundance target.}
#'   \item{\code{[[2]]}}{A \code{data.frame} containing predictions from the
#'     fitted Pareto-front curve, its estimated derivative, and the associated
#'     marginal-gain classification.}
#' }
#'
#' @details
#' The actual spatial footprint of each solution is calculated from the selected
#' cells in \code{best_state$x} and the corresponding values in the external
#' object \code{area_km2}. This allows the plotted footprint to reflect the
#' actual area of partial boundary cells even where the optimisation itself
#' assigns boundary-adjusted costs.
#'
#' The function uses the realised value of \code{N10}, rather than the requested
#' value of \code{K}, to define the Pareto front because a requested target may
#' not always be feasible or found by the optimiser.
#'
#' After retaining the minimum-footprint solution for each realised value of
#' \code{N10}, the footprint is forced to be monotonically non-decreasing using
#' \code{cummax()}. This corrects small optimisation inefficiencies where a
#' solution found for a larger species target has a smaller footprint than a
#' previously identified solution for a smaller target.
#'
#' A monotone decreasing shape-constrained additive model is fitted to the
#' relationship between the proportional spatial footprint and the proportion
#' of species not meeting the target. Numerical derivatives of this fitted
#' relationship are calculated using the model linear-predictor matrix. Sections
#' of the curve are then classified as \code{"Strong gains"},
#' \code{"Intermediate gains"}, or \code{"Moderate gains"} according to
#' quantiles of the derivative.
#'
#' The function currently assumes that the objects \code{area_km2}, \code{nS},
#' and \code{total_area} are available in the calling environment. Here,
#' \code{nS} is the total number of species and \code{total_area} is the total
#' area of the planning domain in the same units as \code{area_km2}.
#'
#' @examples
#' \dontrun{
#' pareto_out <- build_pareto_front_plot(
#'   pareto_front_on_K = pareto_front_on_K_b1_e1
#' )
#'
#' pareto_out[[1]]
#' head(pareto_out[[2]])
#' }
#'
#' @importFrom dplyr arrange case_when filter group_by mutate slice_min ungroup
#' @importFrom ggplot2 aes coord_cartesian element_text geom_line geom_point
#'   ggplot scale_colour_manual scale_x_continuous scale_y_continuous theme
#'   theme_bw
#' @importFrom scam scam
#' @importFrom stats coef predict quantile
#'
#' @export
build_pareto_front_plot <- function(pareto_front_on_K) {
  # = Total England area extracted from each K_target i = 1
  pareto_statistics <- lapply(1:length(pareto_front_on_K), function(i) {
    print(i)

    tmp <- pareto_front_on_K[[i]]
    beta <- tmp$best_state$beta
    
    # Map to actual area because best assumes equal cost for partial cells
    best_area_actual_size <- sum(tmp$best_state$x * area_km2)
    
    # State statistics
    state_i <- compute_state_fast(tmp[[2]])

    # Output summary
    data.frame(
      K_target = tmp$K,
      beta = beta,
      N10 = state_i$N10,
      best_footprint_10km = tmp$best_n_cells,
      best_footprint_area_km = best_area_actual_size,
      gm_pct = state_i$gm_pct,
      gm_log = state_i$gm_log,
      U = state_i$U,
      D = round(state_i$D, 4) # handle floating point issue
    )
  })
  pareto_statistics <- do.call(rbind, pareto_statistics)
  pareto_statistics$N10pc <- pareto_statistics$N10 / nS
  # The K_target might not necessarily be possible (or found) so the N10 and
  # K_target values could differ. Therefore, The K_target is not used, and
  # instead the N10 value is the definitive number of species meeting the +10%
  # target.
  
  summary(pareto_statistics)

  # = Get unique N10 lowest area footprint
  pareto_by_N10 <- pareto_statistics %>% 
    dplyr::group_by(N10) %>% 
    dplyr::slice_min(best_footprint_area_km, n = 1, with_ties = FALSE) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(N10)

  # Sometimes the forward step finds a lower footprint, which by definition would
  # also satisfy the early N10, therefore we enforce a monotone to correct these
  # small inefficiencies. Obviously, running the optimisation longer might find a
  # more efficient solution to the earlier N10, but in trials it took a long time
  # and made no real difference to the overall pattern of the pareto front.
  pareto_by_N10$footprint_monotone <- cummax(pareto_by_N10$best_footprint_area_km)

  # ===
  # === Plot Pareto front - minimum footprint vs N10 (minimising D)
  # ===
  
  # Minimum solution along N10
  min_sol_N10 <- pareto_by_N10 %>%
    mutate(
      U = nS - N10,
      y = U / nS,
      x = footprint_monotone / total_area
    )

  # Prepare data
  smooth_df <- min_sol_N10 %>% 
    filter(!is.na(x), !is.na(y))

  # Fit monotone decreasing spline (adjust bs = "mpd" or "mpi" to match your front)
  scam_fit <- scam(y ~ s(x, bs = "mpd", k = 40), data = smooth_df)

  # Identify points where rates of change are greatest
  eps <- diff(range(smooth_df$x)) * 1e-5

  x_seq <- seq(min(smooth_df$x), max(smooth_df$x), length.out = 1000)
  deriv_df <- data.frame(x = x_seq)

  Xp <- predict(scam_fit, newdata = transform(deriv_df, x = x + eps), type = "lpmatrix")
  X0 <- predict(scam_fit, newdata = deriv_df, type = "lpmatrix")
  Xm <- predict(scam_fit, newdata = transform(deriv_df, x = x - eps), type = "lpmatrix")

  Xderiv <- (Xp - Xm) / (2 * eps)

  left_edge <- deriv_df$x - eps < min(smooth_df$x)
  right_edge <- deriv_df$x + eps > max(smooth_df$x)

  Xderiv[left_edge, ] <- (Xp[left_edge, ] - X0[left_edge, ]) / eps
  Xderiv[right_edge, ] <- (X0[right_edge, ] - Xm[right_edge, ]) / eps
  
  deriv_df <- deriv_df %>% 
    mutate(
      y = predict(scam_fit, newdata = deriv_df),
      derivative = as.vector(Xderiv %*% coef(scam_fit)),
      rate_class = cut(
        derivative,
        breaks = quantile(derivative, probs = c(0, 0.1, 0.5, 1), na.rm = TRUE),
        labels = c("Strong gains", "Intermediate gains", "Moderate gains"),
        include.lowest = TRUE
      ),
      highlighted = !is.na(rate_class),
      group = cumsum(
        c(
          0,
          diff(as.numeric(rate_class)) != 0 |
            diff(highlighted) != 0
        )
      )
    )

  deriv_df_plot <- deriv_df %>% 
    filter(highlighted)

  # Predict over a fine grid
  x_seq <- seq(min(smooth_df$x), max(smooth_df$x), length.out = 200)
  pred_df <- data.frame(x = x_seq) |>
    mutate(y = predict(scam_fit, newdata = data.frame(x = x_seq)))

  # == Create pareto front plot
  x <- with(min_sol_N10, x)
  y <- with(min_sol_N10, y)
  ok <- is.finite(x) & is.finite(y)
  
  # === For setting breaks
  nice_by <- function(x, n_breaks = 6) {
    rng <- diff(range(x, na.rm = TRUE))
    raw_by <- rng / n_breaks
    
    magnitude <- 10 ^ floor(log10(raw_by))
    residual <- raw_by / magnitude
    
    nice_residual <- dplyr::case_when(
      residual <= 1 ~ 1,
      residual <= 2 ~ 2,
      residual <= 5 ~ 5,
      TRUE ~ 10
    )
    
    nice_residual * magnitude
  }
  
  x_lim <- range(x[ok], na.rm = TRUE)
  y_lim <-  c(0,0.525)
  
  x_by <- nice_by(x_lim, n_breaks = 6)
  y_by <- nice_by(y_lim, n_breaks = 6)

  # Now plot
  pareto_plot <- ggplot(min_sol_N10) +
    geom_line(
      data = deriv_df_plot,
      aes(x = x, y = y, group = group, colour = rate_class),
      linewidth = 5,
      alpha = 0.5,
      lineend = "round"
    ) +
    scale_colour_manual(
      values = c(
        "Strong gains" = "firebrick",
        "Intermediate gains" = "goldenrod",
        "Moderate gains" = "dodgerblue3"
      ),
      name = "Marginal gains\n(rate of change)"
    ) +
    geom_point(aes(y = y, x = x),
      colour = "grey50",
      size = 1.2, alpha = 0.5
    ) +
    geom_line(aes(y = y, x = x), data = pred_df) +
    scale_y_continuous(
      name = "Species reaching +10% abundance target (%)",
      breaks = seq(
        floor(y_lim[1] / y_by) * y_by,
        ceiling(y_lim[2] / y_by) * y_by,
        by = y_by
      ),
      labels = function(x) 100 * (1 - x)
    ) +
    scale_x_continuous(
      name = "Footprint of planning domain (% land area)",
      breaks = seq(
        floor(x_lim[1] / x_by) * x_by,
        ceiling(x_lim[2] / x_by) * x_by,
        by = x_by
      ),
      labels = function(x) 100 * x
    ) +
    coord_cartesian(
      xlim = range(x[ok]),
      ylim = c(0,0.525)#range(y[ok])
    ) +
    theme_bw(base_size = 14) %+replace%
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.8, 0.8),
      legend.text = element_text(size = 12)
    )

  # Return
  list(pareto_plot, deriv_df)
  
}



