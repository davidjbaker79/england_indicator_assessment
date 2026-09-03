#' Add progressive species-target inclusion columns
#'
#' Adds binary columns indicating which species meet a specified abundance
#' target at progressively larger points along a Pareto front defined by the
#' number of species targeted.
#'
#' For each value in \code{proportions}, the function identifies the Pareto-front
#' solution whose target number of species is closest to that proportion of the
#' full species pool. Species whose multiplier \code{m_s} meets or exceeds
#' \code{threshold} in that solution are assigned a value of 1 in the
#' corresponding column; all other species retain a value of 0.
#'
#' @param asai_info_compiled A \code{data.table} or object coercible to a
#'   \code{data.table} containing species-level information. Must include a
#'   \code{BDI_NAME_CHECK} column matching names in \code{species_names}.
#' @param pareto_front_on_K A list of Pareto-front solutions. Each element used
#'   by the function must contain \code{$best_state$m_s}, a numeric vector of
#'   species-level abundance multipliers corresponding to \code{species_names}.
#' @param N10_list_index Numeric vector giving the target number of species
#'   associated with each element of \code{pareto_front_on_K}.
#' @param species_names Character vector of species names corresponding, in
#'   order, to values of \code{m_s} in each Pareto-front solution.
#' @param proportions Numeric vector of proportions of the total species pool
#'   for which progressive target columns should be generated. Defaults to
#'   \code{seq(0.3, 1, by = 0.1)}.
#' @param prefix Character string used as the prefix for generated column names.
#'   With the default \code{"N_"}, proportions of 0.3 and 1 generate columns
#'   \code{N_30} and \code{N_100}, respectively.
#' @param threshold Numeric scalar giving the minimum value of \code{m_s} used
#'   to classify a species as meeting its target. Defaults to \code{1.1}.
#'
#' @return Invisibly returns \code{asai_info_compiled} with additional binary
#'   columns indicating whether each species meets the target at each
#'   progressive Pareto-front solution.
#'
#' @details
#' For each requested proportion \eqn{p}, the function identifies the element of
#' \code{N10_list_index} closest to \eqn{pS}, where \eqn{S} is the total number
#' of species in \code{species_names}. Species satisfying
#' \code{m_s >= threshold} in the corresponding Pareto-front solution are then
#' marked as included.
#'
#' Because \code{data.table::setDT()} modifies objects by reference,
#' \code{asai_info_compiled} may also be modified in the calling environment.
#'
#' @examples
#' \dontrun{
#' asai_info_compiled <- add_progressive_targets(
#'   asai_info_compiled = asai_info_compiled,
#'   pareto_front_on_K = pareto_front_on_K,
#'   N10_list_index = N10_list_index,
#'   species_names = species_names
#' )
#' }
#'
#' @importFrom data.table setDT
#' @export
add_progressive_targets <- function(
    asai_info_compiled,
    pareto_front_on_K,
    N10_list_index,
    species_names,
    proportions = seq(0.3, 1, by = 0.1),
    prefix = "N_",
    threshold = 1.1
) {
  setDT(asai_info_compiled)
  
  n_species <- length(species_names)
  for (p in proportions) {
    col_nm <- paste0(prefix, sprintf("%02d", round(p * 100)))
    
    if (!col_nm %in% names(asai_info_compiled)) {
      asai_info_compiled[, (col_nm) := 0]
    }
    
    ind <- which.min(abs(N10_list_index - n_species * p))[1]
    spp_selected <- species_names[
      pareto_front_on_K[[ind]]$best_state$m_s >= threshold
    ]
    
    asai_info_compiled[
      BDI_NAME_CHECK %in% spp_selected,
      (col_nm) := 1
    ]
    
  }
  
  asai_info_compiled[asai_info_compiled$N_100 == 0,]
  spp_selected[agrep(spp_selected, pattern ="limosa")]
  
  invisible(asai_info_compiled)
}

#' Plot progressive prop species meeting targets by class
#' phy_select = asai_info_compiled_1174
plot_progressive_taxon_stack <- function(
    phy_select,
    taxon_col = "CLASS",
    target_cols = paste0("N_", seq(30, 100, by = 10)),
    as_proportion = TRUE,
    return_dt =FALSE
) {
  library(data.table)
  library(ggplot2)
  
  dt <- copy(as.data.table(phy_select))
  
  dt[
    ,
    first_included := fifelse(
      rowSums(.SD, na.rm = TRUE) == 0,
      "Not included",
      names(.SD)[max.col(.SD == 1, ties.method = "first")]
    ),
    .SDcols = target_cols
  ]
  
  dt[
    ,
    first_rank := match(first_included, target_cols)
  ]
  
  taxon_levels <- dt[
    ,
    .(mean_rank = mean(first_rank, na.rm = TRUE)),
    by = taxon_col
  ][
    order(mean_rank),
    get(taxon_col)
  ]
  
  # Number of species per taxon
  taxon_n <- dt[
    ,
    .N,
    by = taxon_col
  ]
  
  taxon_n[
    ,
    taxon_label := paste0(get(taxon_col), " (", N, ")")
  ]
  
  # Named lookup: original taxon name -> taxon label with sample size
  taxon_label_lookup <- setNames(
    taxon_n$taxon_label,
    taxon_n[[taxon_col]]
  )
  
  dt[
    ,
    first_included := factor(
      first_included,
      levels = c(target_cols, "Not included")
    )
  ]
  
  # Add labelled taxon column for plotting
  dt[
    ,
    taxon_label := taxon_label_lookup[as.character(get(taxon_col))]
  ]
  
  dt[
    ,
    taxon_label := factor(
      taxon_label,
      levels = rev(taxon_label_lookup[taxon_levels])
    )
  ]
  
  plot_dt <- dt[
    ,
    .N,
    by = c("taxon_label", "first_included")
  ]
  
  if (return_dt) {
    return(plot_dt)
  }
  
  if (as_proportion) {
    plot_dt[, value := N / sum(N), by = taxon_label]
    y_lab <- "Proportion of species"
  } else {
    plot_dt[, value := N]
    y_lab <- "Number of species"
  }
  
  # legend text
  levels(plot_dt$first_included) <-
    paste0(
      "N = ",
      sub("^N_", "", levels(plot_dt$first_included)),
      "%"
    )
  
  zissou_cols <- hcl.colors(
    n = dplyr::n_distinct(plot_dt$first_included),
    palette = "Zissou 1"
  ) |>
    colorspace::desaturate(amount = 0.3) |>
    colorspace::lighten(amount = 0.1)
  
  
  ggplot(
    plot_dt,
    aes(
      x = taxon_label,
      y = value,
      fill = first_included
    )
  ) +
    geom_col(
      position = position_stack(reverse = TRUE),
      width = 0.75
    ) +
    coord_flip() +
    scale_fill_manual(
      values = rev(zissou_cols),
      guide = guide_legend(
        reverse = FALSE,
        nrow = 2,
        byrow = TRUE,
        theme = theme(
          legend.title.position = "top",
          legend.title = element_text(
            size = 6,
            hjust = 0,
            lineheight = 0.9,
            margin = margin(b = 2, unit = "pt")
          )
        )
      )
    ) +
    scale_y_continuous(
      labels = if (as_proportion) {
        scales::label_number(
          accuracy = 1,
          scale = 100
        )
      } else {
        waiver()
      },
      expand = expansion(mult = c(0, 0.01))
    ) +
    labs(
      x = NULL,
      y = "Species (%)",
      fill = "First included by species target"
    ) +
    theme_bw(
      base_size = 7,
      base_family = "Arial"
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = 6,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 6,
        colour = "black"
      ),
      
      axis.ticks = element_line(linewidth = 0.18),
      axis.ticks.length = grid::unit(1.5, "mm"),
      
      panel.border = element_rect(
        colour = "black",
        linewidth = 0.18
      ),
      panel.grid.major.x = element_line(
        colour = "grey88",
        linewidth = 0.12
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "left",
      legend.box.just = "left",
      
      legend.text = element_text(size = 5.5),
      legend.key.width = grid::unit(3, "mm"),
      legend.key.height = grid::unit(3.5, "mm"),
      legend.spacing.x = grid::unit(0.5, "mm"),
      legend.margin = margin(
        t = 2,
        r = 0,
        b = 0,
        l = 0,
        unit = "pt"
      ),
      
      plot.margin = margin(
        t = 3,
        r = 8,
        b = 3,
        l = 3,
        unit = "pt"
      )
    )
  
}




