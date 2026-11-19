#' Solve a Pareto front incrementally across increasing species targets
#'
#' This function builds a sequence of footprint-minimising solutions for 
#' increasing values of \code{K}, where \code{K} is the number of species 
#' required to reach the abundance threshold of +10\%. The first target in 
#' \code{K_seq} is solved using \code{\link{solve_min_footprint_given_K}}, 
#' and each subsequent target is obtained by incrementally advancing the 
#' previous solution using \code{\link{advance_K_by_one}}.This approach exploits 
#' the largely nested structure of the optimisation problem such that solutions 
#' for \eqn{K+1} are expected to be close to those for \eqn{K}, often differing 
#' by only a small number of added and pruned cells.
#'
#' @param P_over_R Numeric matrix of proportional abundance gains, with rows
#'   corresponding to candidate cells and columns corresponding to species.
#' @param beta Numeric vector of baseline species abundances used to calculate
#'   multiplicative gains and determine whether species exceed the +10\%
#'   threshold.
#' @param GM_target Numeric scalar giving the minimum required geometric mean
#'   abundance gain across species.
#' @param K_seq Integer vector of target numbers of species required to meet the
#'   +10\% threshold. These should normally be supplied in increasing order.
#' @param cost Optional numeric vector of per-cell costs. If \code{NULL}, all
#'   cells are assumed to have unit cost.
#' @param randomise_within_cost_prune Logical; if \code{TRUE}, ties among cells
#'   with equal removal cost during pruning are broken randomly.
#' @param top_species Integer; number of candidate unsatisfied species closest to
#'   the +10\% threshold to consider when advancing from one \code{K} target to
#'   the next.
#' @param beam_width Integer; number of partial repair branches retained at each
#'   depth during the beam search used in \code{\link{advance_K_by_one}}.
#' @param max_depth Integer; maximum number of sequential cell additions allowed
#'   when searching for a repair that achieves the next \code{K} target.
#' @param top_expand Integer; number of candidate cells considered when
#'   expanding each partial repair branch.
#' @param verbose Logical; if \code{TRUE}, progress messages are printed during
#'   the initial solve and each subsequent advance in \code{K}.
#'
#' @return A named list with one element per value in \code{K_seq}. Each element
#'   is itself a list containing:
#'   \describe{
#'     \item{\code{K}}{The target number of species required to exceed the +10\%
#'     threshold.}
#'     \item{\code{best_state}}{The selected solution state for that target, or
#'     \code{NULL} if no feasible solution was found.}
#'     \item{\code{best_footprint}}{The footprint size, defined as the number of
#'     selected cells, or \code{NA_integer_} if no feasible solution was found.}
#'   }
#'
#' @details
#' The first target in \code{K_seq} is solved using the existing full solver.
#' Subsequent targets are then reached incrementally by repeatedly calling
#' \code{\link{advance_K_by_one}} until the current solution satisfies the next
#' requested \code{K}. After each successful advance, the resulting solution is
#' carried forward and used as the starting point for the next target.
#'
#' This strategy is intended to be more efficient than solving each \code{K}
#' independently when the feasible solutions are approximately nested across
#' increasing \code{K}.
#'
#' @seealso
#' \code{\link{solve_min_footprint_given_K}},
#' \code{\link{advance_K_by_one}},
#' \code{\link{prep_sparse_rows}}
#'
#' @examples
#' \dontrun{
#' K_seq <- seq(600, 700, by = 5)
#'
#' res <- solve_pareto_front_incremental(
#'   P_over_R = P_over_R,
#'   beta = beta,
#'   GM_target = 1.02,
#'   K_seq = K_seq,
#'   top_species = 10,
#'   beam_width = 10,
#'   max_depth = 4,
#'   top_expand = 12,
#'   verbose = TRUE
#' )
#'
#' sapply(res, `[[`, "best_footprint")
#' }
#'
#' @export
solve_pareto_front_incremental <- function(P_over_R,
                                           beta,
                                           GM_target,
                                           K_seq,
                                           cost = NULL,
                                           randomise_within_cost_prune = TRUE,
                                           top_species = 10,
                                           beam_width = 10,
                                           max_depth = 4,
                                           top_expand = 12,
                                           overshoot_lambda = 0,
                                           verbose = TRUE) {
  
  # == Set up
  
  # Make sure to use sparse representation of P_over_R
  sparse <- prep_sparse_rows(P_over_R)
  
  # Number of cells
  C <- nrow(P_over_R)
  
  # Set costs if not supplied
  if (is.null(cost)) cost <- rep(1, C)
  
  # Store results (one for each K)
  results <- vector("list", length(K_seq))
  names(results) <- paste0("K_", K_seq)
  
  # == First solve the unconstrained solution
  K0 <- K_seq[1]
  if (verbose) message(sprintf("Initial solve for K = %d", K0))
  
  # This does a full solve for minimum footprint with K0
  st <- solve_min_footprint_given_K(
    P_over_R = P_over_R,
    beta = beta,
    GM_target = GM_target,
    K = K0,
    cost = cost,
    incumbent = NULL,
    randomise_within_cost_prune = randomise_within_cost_prune
  )
  
  # Check whether the solution is valid
  K0_check <- !is.null(st) && ((st$gm_sumlog / length(st$m_s)) >= GM_target) && (st$N10 >= K0)
  
  # K0 results 
  results[[1]] <- list(
    K = K0,
    best_state = if (K0_check) st else NULL,
    best_area = if (K0_check) sum(cost[st$x == 1L]) else NA_real_,
    best_n_cells = if (K0_check) sum(st$x == 1L) else NA_integer_
  )
  
  # == Now, solve other K values initialised by K0 solution
  
  if (length(K_seq) >= 2) {
    
    
    for (i in seq.int(2, length(K_seq))) {
      
      Kt <- K_seq[i]
      if (verbose) message(sprintf("Advancing to K = %d", Kt))
      
      while (!is.null(st) && st$N10 < Kt) {
        
        st <- advance_K_by_one(
          state = st,
          target_K = st$N10 + 1L,
          GM_target = GM_target,
          sparse = sparse,
          cost = cost,
          top_species = top_species,
          beam_width = beam_width,
          max_depth = max_depth,
          top_expand = top_expand,
          overshoot_lambda = overshoot_lambda,
          randomise_within_cost_prune = randomise_within_cost_prune,
          fallback = "none",
          P_over_R = P_over_R
        )
        
      }
      
      # Check valid solution
      K_check <- !is.null(st) && ((st$gm_sumlog / length(st$m_s)) >= GM_target) &&(st$N10 >= Kt)
      
      # K results
      results[[i]] <- list(
        K = Kt,
        best_state = if (K_check) st else NULL,
        best_area = if (K_check) sum(cost[st$x == 1L]) else NA_real_,
        best_n_cells = if (K_check) sum(st$x == 1L) else NA_integer_
      )
    }
  }
  
  results
}


#' Calculate boundary-adjusted cell costs
#'
#' Helper function to account for partial cells along the boundary of a study
#' region. Without adjustment, clipped boundary cells can have substantially
#' lower costs than complete cells and may therefore disproportionately
#' contribute to an optimisation solution.
#'
#' The adjustment increases the cost of partial cells towards the cost of a
#' full cell. The strength of this adjustment is controlled by \code{eta}.
#' Cells whose area is sufficiently close to the full-cell area, as determined
#' by \code{partial_tol}, are treated as complete cells and are not adjusted.
#'
#' @param cell_area Numeric vector giving the area of each cell. Units are
#'   arbitrary but must be consistent with \code{full_cell_area}.
#' @param full_cell_area Optional numeric scalar giving the area of a complete
#'   cell. If \code{NULL}, the maximum value in \code{cell_area} is used.
#' @param eta Numeric scalar controlling the additional cost assigned to partial
#'   cells. \code{eta = 0} gives no adjustment (cost equals actual cell area),
#'   \code{eta = 1} assigns each partial cell the cost of a complete cell, and
#'   values between 0 and 1 provide an intermediate penalty. Default is
#'   \code{1}.
#' @param partial_tol Numeric scalar giving the minimum proportion of
#'   \code{full_cell_area} at which a cell is treated as complete and therefore
#'   left unadjusted. Default is \code{0.999}.
#'
#' @return A numeric vector of boundary-adjusted cell costs with the same
#'   length as \code{cell_area}.
#'
#' @details
#' For cells identified as partial, adjusted cost is calculated as
#'
#' \deqn{C_i = A_i + \eta(A_f - A_i)}
#'
#' where \eqn{A_i} is the observed cell area, \eqn{A_f} is the full-cell area,
#' and \eqn{\eta} controls the strength of the boundary-cell penalty.
#'
#' @examples
#' cell_area <- c(1, 1, 0.75, 0.25)
#'
#' # Treat all partial cells as having the cost of a full cell
#' make_boundary_adjusted_cost(cell_area, full_cell_area = 1, eta = 1)
#'
#' # Apply an intermediate penalty
#' make_boundary_adjusted_cost(cell_area, full_cell_area = 1, eta = 0.5)
#'
#' @export
make_boundary_adjusted_cost <- function(cell_area,
                                        full_cell_area = NULL,
                                        eta = 1,
                                        partial_tol = 0.999) {
  
  # Get full cell area
  if (is.null(full_cell_area)) {
    full_cell_area <- max(cell_area)
  }
  
  area_frac <- pmin(cell_area / full_cell_area, 1)
  
  # By changing eta the relative weight given to partial vs full cells changes:
  #  * eta = 0 gives current behaviour, cost = actual area
  #  * eta = 1 give all partial cells cost as much as full cells
  #  * 0 < eta < 1 gives a partial penalty, but not full equivalence
  cost_adj <- cell_area + eta * pmax(0, full_cell_area - cell_area)
  
  # Only apply to genuinely clipped cells
  cost_adj[area_frac >= partial_tol] <- cell_area[area_frac >= partial_tol]
  
  cost_adj
}


#' Convert weighted presence matrix to sparse row representation
#'
#' Converts a dense cell × species weighted presence matrix (\code{P_over_R})
#' into a sparse row-wise representation to enable fast incremental updates
#' during optimisation. For each cell, the function stores the indices of
#' species present and their corresponding weights.
#'
#' This representation is used to efficiently update species-level state
#' variables (e.g. treated range fraction, abundance multipliers, and
#' indicator values) when adding or removing individual cells during
#' optimisation, avoiding repeated full matrix operations.
#'
#' @param P_over_R A numeric matrix with dimensions \code{C × S}, where rows
#'   correspond to planning units (cells) and columns correspond to species.
#'   Entries represent the contribution of each cell to the treated fraction
#'   of a species' range (typically \code{P_cs / R_s}).
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{idx}{A list of length \code{C}; each element is an integer vector
#'       giving the column indices of species present in that cell.}
#'     \item{wts}{A list of length \code{C}; each element is a numeric vector
#'       of weights corresponding to \code{P_over_R[c, idx[[c]] ]}.}
#'   }
#'
#' @details
#' The sparse representation allows add/remove operations on cells to update
#' only the subset of species affected by that cell, reducing per-move
#' complexity from \eqn{O(S)} to \eqn{O(k)}, where \eqn{k} is the number of
#' species present in the cell.
#'
#' @seealso \code{\link{greedy_meet_gm_fast}},
#'   \code{\link{improve_under_lambda_fast}}
#'
#' @examples
#' \dontrun{
#' sparse <- prep_sparse_rows(P_over_R)
#' idx <- sparse$idx
#' wts <- sparse$wts
#' }
prep_sparse_rows <- function(P_over_R) {
  
  # Number of cells
  C <- nrow(P_over_R)

  # For storing idx and wts
  idx <- vector("list", C)
  wts <- vector("list", C)
  
  # Each cell extract..
  for (c in 1:C) {
    ii <- which(P_over_R[c, ] != 0)
    idx[[c]] <- ii
    wts[[c]] <- P_over_R[c, ii]
  }
  
  # Return a list
  list(idx = idx, wts = wts)
}


#' Solve for a minimum-footprint feasible solution given an N10 target (K)
#'
#' High-level wrapper that (1) builds a feasible solution meeting both the log-GM
#' target and the N10 target \code{K}, (2) prunes the selected set to reduce
#' footprint/cost while preserving feasibility, and (optionally) (3) performs a
#' swap-based local search followed by another pruning pass when \code{K} is near
#' the upper end of the feasible range.
#'
#' @param P_over_R Numeric matrix of dimension \code{C x S}. Rows are cells and
#'   columns are species, containing per-cell contributions used by the state
#'   update functions.
#' @param beta Numeric scalar passed through to \code{build_feasible_for_K()} /
#'   \code{init_state()}.
#' @param GM_target Numeric scalar. Minimum acceptable log-GM value.
#' @param K Integer scalar. Minimum acceptable \code{state$N10} value.
#' @param cost Optional numeric vector of length \code{C} giving per-cell costs.
#'   If \code{NULL}, defaults inside \code{build_feasible_for_K()} to all ones.
#' @param sparse Optional sparse row representation (e.g. from
#'   \code{prep_sparse_rows(P_over_R)}) containing \code{idx} and \code{wts}.
#'   If \code{NULL}, it is computed once and reused for all phases.
#' @param eps_abs,eps_rel Numeric scalars passed to \code{pick_eps_tie()} for
#'   stable tie-breaking during the feasibility-building greedy phase.
#' @param randomise_within_cost_prune Logical. Passed to
#'   \code{prune_feasible_for_K()} (and to the final pruning after swapping) to
#'   randomise removal order within equal-cost blocks.
#'
#' @return A feasible \code{state} object (as produced/updated by your state
#'   functions) or \code{NULL} if feasibility cannot be achieved (i.e., the greedy
#'   build phase exhausts candidates before reaching \code{GM_target} and \code{K}).
#'
#' @details
#' The function runs:
#' \enumerate{
#'   \item \code{\link{build_feasible_for_K}} to construct a feasible state.
#'   \item Feasibility check: returns \code{NULL} if \code{gm_sumlog / S < GM_target}
#'     or \code{N10 < K} (with \code{S = length(m_s)}).
#'   \item \code{\link{prune_feasible_for_K}} to remove redundant selected cells
#'     while maintaining feasibility.
#'   \item Optional refinement: if \code{K > 0.85 * S}, run
#'     \code{\link{swap_explore_then_prune}} (default \code{n_swaps = 10000}) to
#'     explore feasible swaps and then prune again.
#' }
#'
#' @seealso \code{\link{build_feasible_for_K}}, \code{\link{prune_feasible_for_K}},
#'   \code{\link{swap_explore_then_prune}}, \code{\link{prep_sparse_rows}}
#'
#' @export
solve_min_footprint_given_K <- function(P_over_R, beta, GM_target, K, cost=NULL, 
                                        eps_abs=1e-12, eps_rel=1e-6,
                                        randomise_within_cost_prune=FALSE,
                                        incumbent = NULL) {
  
  # Create sparse data
  sparse <- prep_sparse_rows(P_over_R)
  
  # Build the intial feasible solution by greedily adding cells based on marginal 
  # contribution to GM targets
  st <- build_feasible_for_K(P_over_R, beta, GM_target, K,
                             cost=cost, sparse=sparse,
                             eps_abs=eps_abs, eps_rel=eps_rel,
                             incumbent = incumbent)

  # If not feasible (i.e., ran out of candidates), return NULL
  if ((st$gm_sumlog / length(st$m_s)) < GM_target || st$N10 < K) return(NULL)
  
  # Prune initial feasible to K
  st <- prune_feasible_for_K(st, sparse, 
                             GM_target, K, 
                             randomise_within_cost=randomise_within_cost_prune)
  
  # Near boundary it can get harder to find solutions so do more swapping and pruning
  if (K > 0.85 * ncol(P_over_R)) {
    st <- swap_explore_then_prune(st, sparse, GM_target, K,
                                  randomise_within_cost=randomise_within_cost_prune)
  }
  
  # Return
  st
  
}


#' Build a feasible solution for a fixed N10 target (K)
#'
#' Constructs a (greedy) feasible selection vector by adding cells in two phases:
#' (A) first add cells until the geometric-mean objective (GM) reaches `GM_target`,
#' then (B) add cells until the number of species meeting the 10% threshold
#' (`N10`) reaches `K`, as cheaply as possible. Ties in greedy scores are broken
#' using absolute/relative tolerances.
#'
#' @param P_over_R Numeric matrix of dimension `C x S` giving per-cell contributions
#'   by species (typically a ratio such as production over requirement). Rows are
#'   cells and columns are species.
#' @param beta Numeric scalar. Parameter passed to `init_state()` (e.g., controls
#'   smoothing/aggregation in the GM calculation, depending on your implementation).
#' @param GM_target Numeric scalar. Target value for the GM constraint, compared to
#'   `gm_log_from_state(state)` (i.e., the algorithm ensures GM is at least this target,
#'   if possible).
#' @param K Integer scalar. Target for `state$N10` (number of species meeting the 10%
#'   threshold). Phase B adds cells until `state$N10 >= K` or no cells remain.
#' @param cost Optional numeric vector of length `C` giving per-cell costs used to
#'   score candidate additions. Defaults to all ones.
#' @param sparse Optional list as returned by `prep_sparse_rows(P_over_R)`, containing
#'   row-wise sparse structure:
#'   \itemize{
#'     \item `idx`: integer indices of non-zero species per cell
#'     \item `wts`: numeric weights for those species
#'   }
#'   Providing this avoids recomputing sparsity.
#' @param eps_abs,eps_rel Numeric scalars. Absolute and relative tolerances passed
#'   to `pick_eps_tie()` for stable tie-breaking among near-equal greedy scores.
#'
#' @return A `state` object as produced and updated by `init_state()` / `apply_add()`.
#'   At minimum this function assumes the returned state contains:
#'   \itemize{
#'     \item `x`: integer vector of length `C` (0/1) indicating selected cells
#'     \item `cost`: numeric vector of length `C`
#'     \item `N10`: integer scalar tracking the number of species meeting the 10% threshold
#'   }
#'   Additional fields depend on your state implementation.
#'
#' @details
#' The algorithm proceeds greedily:
#' \enumerate{
#'   \item \strong{Phase A (GM feasibility):} repeatedly add the unselected cell that
#'     maximizes \code{(d_gm / S) / cost} until \code{gm_log_from_state(state) >= GM_target}
#'     or no candidates remain.
#'   \item \strong{Phase B (N10 feasibility):} while \code{state$N10 < K}, prefer cells with
#'     positive \code{d_N10}. If none exist, add the cell that best reduces the total shortfall
#'     \code{D} (i.e., most negative \code{d_D}). In both cases a small GM term is used as a
#'     secondary tie-breaker.
#' }
#' If all cells are selected before reaching targets, the loop terminates and the
#' best-achieved state is returned.
#'
#' @seealso \code{\link{prep_sparse_rows}}, \code{\link{init_state}},
#'   \code{\link{delta_add}}, \code{\link{apply_add}}, \code{\link{gm_log_from_state}},
#'   \code{\link{pick_eps_tie}}
#'
#' @export
build_feasible_for_K <- function(P_over_R, beta, GM_target, K,
                                 cost=NULL, sparse=NULL,
                                 eps_abs=1e-12, eps_rel=1e-6,
                                 incumbent = NULL) {
  # Number of cells
  C <- nrow(P_over_R) 
  
  # Number of species
  S <- ncol(P_over_R)
  
  # Add costs (1 by default)
  if (is.null(cost)) cost <- rep(1, C)
  
  # Create sparse if not provided
  if (is.null(sparse)) sparse <- prep_sparse_rows(P_over_R)
  idx <- sparse$idx # Sparse cell species lists
  wts <- sparse$wts # Sparse weights
  
  # Compute initial state or use previous (none for K0)
  if (!is.null(incumbent)) {
    state <- incumbent
  } else {
    state <- init_state(C, S, beta, cost)
  }
  
  # Initial phase to meet GM constraint 
  # Greedy algorithm finds candidate cell with highest marginal gain to GM 
  while (gm_log_from_state(state) < GM_target) {
    
    # Get candidates and if none break (can't do anything)
    # Candidates are all cells for K = 0
    cand <- which(state$x == 0L)
    if (!length(cand)) break
    
    # Calculate potential gains from each candidate cell
    gains <- vapply(cand, function(c) {
      d <- delta_add(state, c, idx, wts)
     (d$d_gm / S) / state$cost[c]
    }, numeric(1))
    
    # Best cell to add (handling ties)
    j <- pick_eps_tie(gains, "max", eps_abs, eps_rel)
    
    # Add cell to state
    state <- apply_add(state, cand[j], idx, wts)
  }
  
  # Second phase aims to meet N10 constraint as cheaply as possible 
  # Add cells untill N10 >= K
  while (state$N10 < K) {
    
    # Candidate cells
    cand <- which(state$x == 0L)
    if (!length(cand)) break
    
    # Get the statistics for adding each candidate
    deltas <- lapply(cand, function(c) delta_add(state, c, idx, wts))
    
    # Extract the d_N10 statistics
    dN <- vapply(deltas, '[[', numeric(1), "d_N10")
    
    # If any dN great than 0
    if (any(dN > 0)) {
      
      # Find moves that increase N10
      score <- vapply(seq_along(cand), function(j) {
        d <- deltas[[j]]
        c <- cand[j]
        (d$d_N10 / state$cost[c]) + 1e-6 * (d$d_gm / S)
      }, numeric(1))
      
      # Best cell to add (handling ties)
      j <- pick_eps_tie(score, "max", eps_abs, eps_rel)
      
      # Add cell to state
      state <- apply_add(state, cand[j], idx, wts)
      
    } else {
      # No immediate N10 gains available so reduce total shortfall D
      # (more negative d_D is better)
      score <- vapply(cand, function(c) {
        d <- delta_add(state, c, idx, wts)
        (-d$d_D / state$cost[c]) + 1e-6 * (d$d_gm / S)
      }, numeric(1))
      
      # Best cell to add (handling ties)
      j <- pick_eps_tie(score, "max", eps_abs, eps_rel)
      
      # Add cell to state
      state <- apply_add(state, cand[j], idx, wts)
    }
  }
  
  # Return state
  state
  
}


#' Compute log-geometric-mean objective from a state object
#'
#' Returns the log geometric mean (GM) value implied by \code{state}, computed as
#' the mean of \code{log(m_s)} over species, where \code{gm_sumlog} is assumed to
#' store \eqn{\sum_s \log(m_s)}.
#'
#' @param state A state object containing at least:
#'   \itemize{
#'     \item \code{gm_sumlog}: numeric scalar, sum of \code{log(m_s)} over species
#'     \item \code{m_s}: numeric vector of per-species values whose log-GM is tracked
#'   }
#'
#' @return Numeric scalar giving \code{gm_sumlog / length(m_s)}.
#'
#' @examples
#' state <- list(gm_sumlog = sum(log(c(1, 2, 4))), m_s = c(1, 2, 4))
#' gm_log_from_state(state)
#'
#' @export
gm_log_from_state <- function(state) {
  
  state$gm_sumlog / length(state$m_s)
  
}


#' Compute incremental change in cached state from adding a cell
#'
#' Computes the incremental (delta) change in cached optimisation quantities
#' that would result from adding a planning unit (cell) to the current
#' solution. The function does not modify the state; it is intended for fast
#' candidate scoring during optimisation.
#'
#' The returned deltas are defined with respect to the cached quantities in
#' \code{state}, and are computed only for the subset of species present in
#' the candidate cell using the sparse row representation (\code{idx}, \code{wts}).
#'
#' @param state A cached optimisation state object (see \code{\link{init_state}}).
#' @param c Integer. Index of the cell to add (\code{state$x[c]} should be 0).
#' @param idx List of length \code{C}. For each cell \code{c}, \code{idx[[c]]}
#'   contains integer indices of species present in that cell (from
#'   \code{\link{prep_sparse_rows}}).
#' @param wts List of length \code{C}. For each cell \code{c}, \code{wts[[c]]}
#'   contains numeric weights corresponding to \code{P_over_R[c, idx[[c]]]}.
#'
#' @return A named list of deltas:
#'   \describe{
#'     \item{d_gm}{Numeric scalar. Change in \code{gm_sumlog}
#'       (i.e., \eqn{\Delta \sum_s \log(m_s)}).}
#'     \item{d_D}{Numeric scalar. Change in total shortfall \code{D}.}
#'     \item{d_N10}{Integer. Change in the number of species meeting the +10\%
#'       target (\code{N10}).}
#'     \item{d_cost}{Numeric scalar. Cost increment associated with adding the
#'       cell (typically \code{state$cost[c]}).}
#'   }
#'
#' @details
#' This function assumes the species response model \eqn{m_s = 1 + \beta p_s},
#' where adding a cell increases \eqn{p_s} by the row weights \eqn{P_{cs}/R_s}
#' for species present in that cell. Changes are computed only for those
#' species, reducing per-candidate cost from \eqn{O(S)} to \eqn{O(k)}, where
#' \eqn{k} is the number of species present in the cell.
#' 
#' @export
delta_add <- function(state, c, idx, wts) {
  
  # Get the species and weights for the candidate cell
  ii <- idx[[c]]
  ww <- wts[[c]]
  
  # If no species for that candidates return zeros
  if (length(ii) == 0) {
    return(list(d_gm=0, d_D=0, d_N10=0, d_cost=state$cost[c]))
  }
  
  # Post treatment abundance multiplier
  old_m <- state$m_s[ii] # old state
  new_p <- state$p_s[ii] + ww # addition of new cell contributions
  new_m <- 1 + state$beta * new_p # new abundance
  
  # How much does this addition improve the solution? 
  d_gm <- sum(log(new_m)) - sum(log(old_m))
  
  # Species target shortfalls
  old_short <- pmax(0, 1.10 - old_m)
  new_short <- pmax(0, 1.10 - new_m)
  d_D <- sum(new_short) - sum(old_short)
  
  # How much does this addition improve the N10?
  d_N10 <- sum(new_m >= 1.10) - sum(old_m >= 1.10)
  
  # Return statistics as a list
  list(d_gm=d_gm, d_D=d_D, d_N10=d_N10, d_cost=state$cost[c])
}


#' Compute incremental change in cached state from removing a cell
#'
#' Computes the incremental (delta) change in cached optimisation quantities
#' that would result from removing a planning unit (cell) from the current
#' solution. The function does not modify the state; it is intended for fast
#' feasibility checks and candidate scoring during optimisation.
#'
#' @param state A cached optimisation state object (see \code{\link{init_state}}).
#' @param c Integer. Index of the cell to remove (\code{state$x[c]} should be 1).
#' @param idx List of length \code{C}. Species indices per cell
#'   (from \code{\link{prep_sparse_rows}}).
#' @param wts List of length \code{C}. Weights per cell corresponding to
#'   \code{P_over_R[c, idx[[c]]]} (from \code{\link{prep_sparse_rows}}).
#'
#' @return A named list of deltas (see \code{\link{delta_add}}). For removals,
#'   \code{d_cost} is negative (cost reduction).
#'
#' @details
#' This function computes changes under the response model \eqn{m_s = 1 + \beta p_s},
#' where removing a cell decreases \eqn{p_s} by the row weights for species
#' present in that cell. Only affected species are updated.
#' 
#' @export
delta_remove <- function(state, c, idx, wts) {
  ii <- idx[[c]]
  ww <- wts[[c]]
  if (length(ii) == 0) return(list(d_gm=0, d_D=0, d_N10=0, d_cost=-state$cost[c]))
  
  old_m <- state$m_s[ii]
  new_p <- state$p_s[ii] - ww
  new_m <- 1 + state$beta * new_p
  
  d_gm <- sum(log(new_m)) - sum(log(old_m))
  
  old_short <- pmax(0, 1.10 - old_m)
  new_short <- pmax(0, 1.10 - new_m)
  d_D <- sum(new_short) - sum(old_short)
  
  d_N10 <- sum(new_m >= 1.10) - sum(old_m >= 1.10)
  
  list(d_gm=d_gm, d_D=d_D, d_N10=d_N10, d_cost=-state$cost[c])
}


#' Apply addition of a cell to the cached optimisation state
#'
#' Updates (mutates) a cached optimisation state object by adding a planning
#' unit (cell) to the current solution. Cached quantities (\code{p_s}, \code{m_s},
#' \code{gm_sumlog}, \code{D}, and \code{N10}) are updated incrementally using the
#' sparse row representation, avoiding full recomputation.
#'
#' @param state A cached optimisation state object (see \code{\link{init_state}}).
#' @param c Integer. Index of the cell to add.
#' @param idx List of species indices per cell (from \code{\link{prep_sparse_rows}}).
#' @param wts List of weights per cell (from \code{\link{prep_sparse_rows}}).
#'
#' @return The updated \code{state} object (invisibly in typical use).
#'
#' @details
#' This function assumes \code{state$x[c] == 0} at the time of addition.
#' No explicit check is performed for speed; callers should enforce this.
#' 
#' @export
apply_add <- function(state, c, idx, wts) {
  
  # Update state based on added cell
  ii <- idx[[c]]
  ww <- wts[[c]]
  if (length(ii) > 0) {
    
    # Species abundance increase
    old_m <- state$m_s[ii]
    state$p_s[ii] <- state$p_s[ii] + ww
    state$m_s[ii] <- 1 + state$beta * state$p_s[ii]
    new_m <- state$m_s[ii]
    
    # GM target
    state$gm_sumlog <- state$gm_sumlog - sum(log(old_m)) + sum(log(new_m))
    
    # Species target shortfall
    old_short <- pmax(0, 1.10 - old_m)
    new_short <- pmax(0, 1.10 - new_m)
    state$D <- state$D - sum(old_short) + sum(new_short)
    
    # N10
    state$N10 <- state$N10 - sum(old_m >= 1.10) + sum(new_m >= 1.10)
  }
  
  # Add the cell
  state$x[c] <- 1L
  
  # Return state
  state
}


#' Apply removal of a cell from the cached optimisation state
#'
#' Updates (mutates) a cached optimisation state object by removing a planning
#' unit (cell) from the current solution. Cached quantities (\code{p_s}, \code{m_s},
#' \code{gm_sumlog}, \code{D}, and \code{N10}) are updated incrementally using the
#' sparse row representation.
#'
#' @param state A cached optimisation state object (see \code{\link{init_state}}).
#' @param c Integer. Index of the cell to remove.
#' @param idx List of species indices per cell (from \code{\link{prep_sparse_rows}}).
#' @param wts List of weights per cell (from \code{\link{prep_sparse_rows}}).
#'
#' @return The updated \code{state} object.
#'
#' @details
#' This function assumes \code{state$x[c] == 1} at the time of removal.
#' No explicit check is performed for speed; callers should enforce this.
#'
#' @export
apply_remove <- function(state, c, idx, wts) {
  
  # Candidate species and weights
  ii <- idx[[c]]
  ww <- wts[[c]]
  
  # Update state given candidate to be removed
  if (length(ii) > 0) {
    
    # Abundance decrease
    old_m <- state$m_s[ii]
    state$p_s[ii] <- state$p_s[ii] - ww
    state$m_s[ii] <- 1 + state$beta * state$p_s[ii]
    new_m <- state$m_s[ii]
    
    # GM target
    state$gm_sumlog <- state$gm_sumlog - sum(log(old_m)) + sum(log(new_m))
    
    # Species targets
    old_short <- pmax(0, 1.10 - old_m)
    new_short <- pmax(0, 1.10 - new_m)
    state$D <- state$D - sum(old_short) + sum(new_short)
    
    # N10
    state$N10 <- state$N10 - sum(old_m >= 1.10) + sum(new_m >= 1.10)
  }
  
  # Remove cell
  state$x[c] <- 0L
  
  # Return
  state
}

#' Prune a feasible solution while preserving GM and N10 targets
#'
#' Attempts to reduce the selected set of cells in a feasible \code{state} by
#' greedily removing cells while maintaining feasibility with respect to:
#' \itemize{
#'   \item the log-GM constraint: \code{gm_log_from_state(state) >= GM_target}
#'   \item the N10 constraint: \code{state$N10 >= K}
#' }
#'
#' The pruning loop runs for up to \code{max_passes} passes, and within each pass
#' considers removals in an "expensive-first" order (highest cost removed first),
#' optionally randomising within equal-cost blocks (...but thought overly complex
#' for this analysis).
#'
#' @param state A state object containing at least \code{x}, \code{cost},
#'   \code{gm_sumlog}, \code{m_s}, and \code{N10}.
#' @param sparse A sparse row representation (e.g. from \code{prep_sparse_rows()})
#'   containing \code{idx} and \code{wts}.
#' @param GM_target Numeric scalar. Minimum acceptable log-GM value.
#' @param K Integer scalar. Minimum acceptable \code{state$N10} value.
#' @param randomise_within_cost Logical. If \code{TRUE}, randomises the removal
#'   order within blocks of identical cost, while preserving the overall
#'   expensive-first structure.
#' @param max_passes Integer. Maximum number of pruning passes. Stops early if
#'   a full pass removes nothing.
#'
#' @return Updated \code{state} after greedy pruning.
#'
#' @details
#' For each selected cell \code{c}, the function computes removal deltas via
#' \code{delta_remove(state, c, idx, wts)}. The removal is applied only if:
#' \itemize{
#'   \item \code{(state$gm_sumlog + d$d_gm) / S >= GM_target}
#'   \item \code{(state$N10 + d$d_N10) >= K}
#' }
#' 
#' @export
prune_feasible_for_K <- function(state, sparse, GM_target, K,
                                 randomise_within_cost=FALSE, 
                                 max_passes=100) {
  
  # Get sparse ids and weights
  idx <- sparse$idx
  wts <- sparse$wts
  
  # Species number
  S <- length(state$m_s)
  
  # Run across passes upto max pass = 1
  for (pass in 1:max_passes) {
    
    # Selected candidates
    cand_rem <- which(state$x == 1L)
    if (!length(cand_rem)) break
    
    # Rank baseline by expensive-first
    ord <- cand_rem[order(state$cost[cand_rem], decreasing = TRUE)]
    
    # Randomise within equal-cost blocks (preserves expensive-first structure)
    if (randomise_within_cost && length(ord) > 1) {
      cost_ord <- state$cost[ord]
      rl <- rle(cost_ord)
      starts <- cumsum(c(1, head(rl$lengths, -1)))
      blocks <- lapply(seq_along(rl$lengths), function(i) {
        ord[starts[i]:(starts[i] + rl$lengths[i] - 1)]
      })
      ord <- unlist(
        lapply(blocks, function(b) {
          b[sample.int(length(b))]
        }),
        use.names = FALSE
      )
      stopifnot(
        identical(sort(ord), sort(cand_rem))
      )
    }
    
    # Check candidates for removal 
    removed_any <- FALSE
    for (c in ord) {
      
      # Skip if already removed earlier in this pass
      if (state$x[c] == 0L) next
      
      # Find the best candidate for removal 
      d <- delta_remove(state, c, idx, wts)
      
      # Check that the GM target is still feasibility after removal
      if ((state$gm_sumlog + d$d_gm) / S < GM_target) next
      
      # Check that the N10 target is feasibility after removal
      if ((state$N10 + d$d_N10) < K) next
      
      # If feasible then remove
      state <- apply_remove(state, c, idx, wts)
      removed_any <- TRUE
    }
    
    # Break if none to remove
    if (!removed_any) break
  }
  
  # Return updated state
  state
  
}

#' Explore feasible swaps, then prune to tighten the solution footprint
#'
#' Performs a swap-based local search starting from a feasible \code{state}, then
#' calls \code{prune_feasible_for_K()} to attempt further cost/footprint reduction.
#'
#' Each iteration proposes swapping one currently selected cell out (\code{c_out})
#' with one currently unselected cell in (\code{c_in}), optionally constrained to
#' be a same-cost swap when costs vary. Proposed swaps must maintain feasibility:
#' \code{(gm_sumlog / S) >= GM_target} and \code{N10 >= K}. Among feasible swaps,
#' the algorithm accepts improvements in a squared-penalty objective, and can
#' optionally accept worse swaps under a simulated-annealing rule.
#'
#' @param state A feasible state object containing at least \code{x}, \code{cost},
#'   \code{gm_sumlog}, \code{m_s}, and \code{N10}.
#' @param sparse A sparse row representation (e.g. from \code{prep_sparse_rows()})
#'   containing \code{idx} and \code{wts}.
#' @param GM_target Numeric scalar. Minimum acceptable log-GM value.
#' @param K Integer scalar. Minimum acceptable \code{state$N10} value.
#' @param n_swaps Integer. Maximum number of swap proposals to evaluate.
#' @param cost Optional numeric vector of length \code{C}. If \code{NULL}, uses
#'   \code{state$cost}. When costs are not all equal, swaps are restricted to
#'   same-cost replacements (\code{cost[c_in] == cost[c_out]}).
#' @param w_gm Numeric scalar. Weight for the GM deviation term in the penalty.
#' @param w_n10 Numeric scalar. Weight for the N10 deviation term in the penalty.
#'   Typically larger than \code{w_gm} so N10 deviations are discouraged more strongly.
#' @param T0 Numeric scalar. Initial annealing temperature. If \code{0}, the search
#'   is deterministic (only accepts non-worsening swaps). If \code{> 0}, a worsening
#'   swap with penalty increase \eqn{\Delta} is accepted with probability
#'   \code{exp(-Delta / T)}, where \code{T = T0 / sqrt(i)} at iteration \code{i}.
#' @param randomise_within_cost Logical passed to \code{prune_feasible_for_K()}.
#'
#' @return A (typically improved) feasible \code{state} after swap exploration and pruning.
#'
#' @details
#' The penalty function used for comparing feasible states is:
#' \deqn{
#'   w_gm ((gm\_sumlog / S) - GM\_target)^2 + w\_n10 (N10 - K)^2
#' }
#' where \code{S = length(state$m_s)}.
#' 
#' @export
swap_explore_then_prune <- function(state, 
                                    sparse, 
                                    GM_target, 
                                    K,
                                    n_swaps_max = 10000, 
                                    cost = NULL,
                                    w_gm = 1, 
                                    w_n10 = 5,
                                    T0 = 0, # 0 means deterministic
                                    randomise_within_cost = TRUE,
                                    top_m_in = 10,
                                    min_swaps = 200,
                                    stall_accept = NULL,
                                    stall_improve = NULL,
                                    window = 200,
                                    min_accept_rate = 0.01,
                                    verbose = FALSE) {
  
  # Get sparse ids and weights
  idx <- sparse$idx 
  wts <- sparse$wts
  
  # Species number
  S <- length(state$m_s)
  
  # Add costs
  if (is.null(cost)) cost <- state$cost
  
  # Footprint
  footprint0 <- sum(cost[state$x == 1L])
  
  # Values to accept solution if it stalls
  if (is.null(stall_accept)) {
    stall_accept <- max(200L, 5L * footprint0)
  }
  if (is.null(stall_improve)) {
    stall_improve <- max(500L, 10L * footprint0)
  }
  
  # Squared penalty
  penalty <- function(st) {
    sgm  <- (st$gm_sumlog / S) - GM_target
    sn10 <- st$N10 - K
    w_gm * sgm^2 + w_n10 * sn10^2
  }
  p0 <- penalty(state)
  best_p <- p0
  
  # Trackers
  n_accept <- 0L
  n_improve <- 0L
  since_accept <- 0L
  since_improve <- 0L
  
  # Record of recent acceptances
  recent_accept <- rep(FALSE, window)
  recent_ptr <- 1L
  
  for (i in seq_len(n_swaps_max)) {
    in_cells <- which(state$x == 1L)
    out_cells <- which(state$x == 0L)
    if (!length(in_cells) || !length(out_cells)) break
    
    accepted <- FALSE
    
    # Choose c_out based on delta_remove
    d_out <- lapply(in_cells, function(c) delta_remove(state, c, idx, wts))
    dN_out <- vapply(d_out, `[[`, numeric(1), "d_N10")
    
    # Prefer removals that do not reduce N10
    safe <- which(dN_out == 0)
    cand_out <- if (length(safe)) in_cells[safe] else in_cells
    d_out_use <- if (length(safe)) d_out[safe] else d_out
    
    # Weight by low GM damage
    w_out <- vapply(d_out_use, function(d) exp(-abs(d$d_gm)), numeric(1))
    if (sum(w_out) > 0) {
      j_out <- sample.int(length(cand_out), 1, prob = w_out)
      c_out <- cand_out[j_out]
      
      # Remove temporarily
      st_rm <- apply_remove(state, c_out, idx, wts)
      
      # Choose c_in conditionally on what removal broke
      cand_in <- out_cells
      if (any(cost != cost[1])) {
        cand_in <- cand_in[cost[cand_in] == cost[c_out]]
      }
      
      if (length(cand_in)) {
        crit_species <- which(st_rm$m_s < 1.10)
        
        if (length(crit_species)) {
          hits_crit <- vapply(
            cand_in,
            function(c) any(idx[[c]] %in% crit_species),
            logical(1)
          )
          if (any(hits_crit)) cand_in <- cand_in[hits_crit]
        }
        
        if (length(cand_in)) {
          d_in <- lapply(cand_in, function(c) delta_add(st_rm, c, idx, wts))
          
          score_in <- vapply(seq_along(cand_in), function(j) {
            d <- d_in[[j]]
            1000 * d$d_N10 + 100 * (-d$d_D) + (d$d_gm / S)
          }, numeric(1))
          
          ord <- order(score_in, decreasing = TRUE)
          top <- ord[seq_len(min(top_m_in, length(ord)))]
          j_in <- sample.int(length(top), 1)
          c_in <- cand_in[top[j_in]]
          
          st2 <- apply_add(st_rm, c_in, idx, wts)
          
          # Final feasibility check
          if ((st2$gm_sumlog / S) >= GM_target && st2$N10 >= K) {
            p1 <- penalty(st2)
            
            if (p1 <= p0) {
              state <- st2
              p0 <- p1
              accepted <- TRUE
            } else if (T0 > 0) {
              T <- T0 / sqrt(i)
              if (runif(1) < exp(-(p1 - p0) / T)) { 
                state <- st2
                p0 <- p1
                accepted <- TRUE
              }
            }
          }
        }
      }
    }
    
    # Update counters
    if (accepted) {
      n_accept <- n_accept + 1L
      since_accept <- 0L
      
      if (p0 < best_p - 1e-12) {
        best_p <- p0
        n_improve <- n_improve + 1L
        since_improve <- 0L
      } else {
        since_improve <- since_improve + 1L
      }
    } else {
      since_accept <- since_accept + 1L
      since_improve <- since_improve + 1L
    }
    
    recent_accept[recent_ptr] <- accepted
    recent_ptr <- if (recent_ptr == window) 1L else recent_ptr + 1L
    
    # Stopping rule
    if (i >= min_swaps) {
      accept_rate <- mean(recent_accept)
      
      stop_now <- FALSE
      
      if (since_accept >= stall_accept) {
        stop_now <- TRUE
        if (verbose) {
          message(sprintf(
            "Stopping: no accepted swap for %d attempts (i=%d)",
            since_accept, i
          ))
        }
      }
      
      if (!stop_now && since_improve >= stall_improve) {
        stop_now <- TRUE
        if (verbose) {
          message(sprintf(
            "Stopping: no penalty improvement for %d attempts (i=%d)",
            since_improve, i
          ))
        }
      }
      
      if (!stop_now && since_accept >= window && accept_rate < min_accept_rate) {
        stop_now <- TRUE
        if (verbose) {
          message(sprintf(
            "Stopping: low recent acceptance rate %.3f at i=%d",
            accept_rate, i
          ))
        }
      }
      
      if (stop_now) break
    }
  }
  
  # Prune again
  prune_feasible_for_K(state, sparse, GM_target, K,
                       randomise_within_cost = randomise_within_cost)
  
}


#' Pick an index among near-optimal scores (epsilon tie-breaking)
#'
#' Identifies the optimum score (maximum or minimum), then treats any candidate
#' within an absolute/relative epsilon band of that optimum as a tie, and returns
#' one tied index sampled uniformly at random.
#'
#' @param scores Numeric vector of candidate scores (non-empty).
#' @param direction Character scalar; whether larger or smaller is better.
#'   One of \code{"max"} (default) or \code{"min"}.
#' @param eps_abs Numeric scalar. Absolute tolerance for determining ties.
#' @param eps_rel Numeric scalar. Relative tolerance (scaled by \code{abs(best)})
#'   for determining ties.
#'
#' @return Integer scalar index into \code{scores} giving the selected candidate.
#'
#' @details
#' Let \code{best} be \code{max(scores)} if \code{direction = "max"} (or
#' \code{min(scores)} if \code{direction = "min"}). Define
#' \code{eps = max(eps_abs, eps_rel * abs(best))}. A score is considered tied if:
#' \itemize{
#'   \item \code{scores >= best - eps} for \code{direction = "max"}
#'   \item \code{scores <= best + eps} for \code{direction = "min"}
#' }
#' One index is then drawn uniformly from the tied set using \code{sample.int()}.
#'
#' @export
pick_eps_tie <- function(scores, direction = c("max", "min"),
                         eps_abs = 1e-12, eps_rel = 1e-6) {
  
  # Do we want min or max
  direction <- match.arg(direction)
  
  if (!length(scores)) stop("scores must be non-empty")
  
  # Break any ties
  if (direction == "max") {
    best <- max(scores)
    eps <- max(eps_abs, eps_rel * abs(best))
    tie <- which(scores >= best - eps)
  } else {
    best <- min(scores)
    eps <- max(eps_abs, eps_rel * abs(best))
    tie <- which(scores <= best + eps)
  }
  
  # Sample among ties (uniform)
  tie[sample.int(length(tie), 1)]
  
}


#' Initialise optimisation state object
#'
#' Creates and returns an initial state object for the spatial optimisation
#' algorithm, representing a solution with no selected cells. The state
#' stores cached species-level and indicator-level quantities that allow
#' fast incremental updates when cells are added or removed during
#' optimisation.
#'
#' The state object is designed to be mutated in place by add/remove
#' operations and avoids repeated recomputation of global summaries such
#' as treated range fractions, abundance multipliers, and geometric-mean
#' indicator values.
#'
#' @param C Integer. Number of planning units (cells).
#' @param S Integer. Number of species.
#' @param beta Numeric scalar. Response scaling parameter linking the
#'   fraction of a species' occupied range treated to its proportional
#'   abundance change (\eqn{m_s = 1 + \beta p_s}).
#' @param cost Optional numeric vector of length \code{C} giving per-cell
#'   costs. If \code{NULL}, all cells are assigned equal cost.
#'
#' @return A list representing the initial optimisation state, with elements:
#'   \describe{
#'     \item{x}{Integer vector of length \code{C} indicating selected
#'       (\code{1}) and unselected (\code{0}) cells.}
#'     \item{cost}{Numeric vector of per-cell costs.}
#'     \item{p_s}{Numeric vector of length \code{S} giving the fraction of
#'       each species' occupied range currently treated.}
#'     \item{m_s}{Numeric vector of length \code{S} giving species-specific
#'       abundance multipliers.}
#'     \item{gm_sumlog}{Numeric scalar equal to
#'       \eqn{\sum_s \log(m_s)}, used to compute the geometric-mean indicator
#'       efficiently.}
#'     \item{D}{Numeric scalar giving the total shortfall from the +10\%
#'       species-level abundance target across all species.}
#'     \item{N10}{Integer. Number of species meeting or exceeding the +10\%
#'       abundance target.}
#'     \item{beta}{Numeric scalar, stored for convenience during updates.}
#'   }
#'
#' @details
#' The initial state corresponds to a solution with no conservation action
#' applied (\code{x = 0} for all cells). All species therefore have
#' \eqn{p_s = 0} and \eqn{m_s = 1}. The total shortfall \code{D} is initially
#' equal to \code{0.10 × S}, reflecting the +10\% target applied uniformly
#' across species.
#'
#' @export
init_state <- function(C, S, beta, cost = NULL) {
  
  # Apply costs
  if (is.null(cost)) cost <- rep(1, C)
  
  # Selected 
  x <- rep(0L, C)
  
  # Fraction of each species range treated
  p_s <- rep(0.0, S)
  
  # Species specific abundance multiplier
  m_s <- rep(1.0, S)  # 1 + beta*0
  
  # GM log
  gm_sumlog <- sum(log(m_s))  # = 0
  
  # Per species shortfall
  shortfall <- pmax(0, 1.10 - m_s)  # = 0.10 per species initially
  
  # Total shortfall
  D <- sum(shortfall)
  
  # N10 
  N10 <- sum(m_s >= 1.10)  # 0 initially
  
  # Return state
  list(
    x = x, cost = cost,
    p_s = p_s, m_s = m_s,
    gm_sumlog = gm_sumlog, D = D, N10 = N10,
    beta = beta
  )
  
}


#' Extract summary metrics from cached optimisation state
#'
#' Computes species-level and indicator-level summary metrics from a cached
#' optimisation state object without recomputing treated range fractions or
#' abundance multipliers. This function is intended for reporting and
#' diagnostics, not for use inside candidate-evaluation loops.
#'
#' @param state A state object created by \code{\link{init_state}} and updated
#'   by incremental add/remove operations (e.g. \code{\link{apply_add}},
#'   \code{\link{apply_remove}}).
#'
#' @return A list containing summary quantities derived from the cached state:
#'   \describe{
#'     \item{p_s}{Numeric vector giving the fraction of each species' occupied
#'       range currently treated.}
#'     \item{m_s}{Numeric vector giving species-specific abundance multipliers.}
#'     \item{gm_log}{Numeric scalar giving the log-scale geometric mean abundance
#'       indicator value.}
#'     \item{gm_pct}{Numeric scalar giving the proportional change in the
#'       geometric mean abundance indicator relative to baseline.}
#'     \item{D}{Numeric scalar giving the total shortfall from the +10\%
#'       species-level abundance target across all species.}
#'     \item{U}{Integer giving the number of species not meeting the +10\%
#'       abundance target.}
#'     \item{N10}{Integer giving the number of species meeting or exceeding the
#'       +10\% abundance target.}
#'   }
#'
#' @details
#' This function relies on cached values stored in the optimisation state
#' (including \code{gm_sumlog}, \code{D}, and \code{N10}) and therefore runs in
#' constant time with respect to the number of species. It should be used for
#' summarising completed solutions rather than during optimisation, where
#' incremental update functions provide substantially greater efficiency.
#'
#' @export
compute_state_fast <- function(state) {
  
  # Number of species
  S <- length(state$m_s)
  
  # GM on log scale
  gm_log <- state$gm_sumlog / S
  
  # GM percentage
  gm_pct <- exp(gm_log) - 1
  
  # Unmet species targets
  U <- S - state$N10
  
  # Return state
  list(
    p_s = state$p_s,
    m_s = state$m_s,
    gm_log = gm_log,
    gm_pct = gm_pct,
    D = state$D,
    U = U,
    N10 = state$N10
  )
  
}

#' Advance an existing solution to the next species target
#'
#' Starting from a current feasible solution state, attempts to construct a new
#' solution that satisfies a higher target number of species exceeding the +10\%
#' abundance threshold. The function does this by identifying candidate
#' unsatisfied species close to the threshold, applying targeted repair using a
#' beam search over possible cell additions, and pruning the resulting solution
#' where possible.
#'
#' This is intended for incremental optimisation where the solution for
#' \eqn{K+1} is expected to be similar to the solution for \eqn{K}, differing by
#' only a small number of added and removed cells.
#'
#' @param state A solution state object containing the current selected cells and
#'   cached quantities such as species-level gains, current \code{N10}, geometric
#'   mean summary, and any other fields required by the update helpers.
#' @param target_K Integer; target number of species required to exceed the
#'   abundance threshold.
#' @param GM_target Numeric scalar giving the minimum required geometric mean
#'   abundance gain across species.
#' @param sparse Sparse row representation of the candidate cell by species gain
#'   matrix, as returned by \code{\link{prep_sparse_rows}}.
#' @param cost Optional numeric vector of per-cell costs. If \code{NULL}, the
#'   function uses \code{state$cost}.
#' @param top_species Integer; number of unsatisfied species closest to the
#'   threshold to consider as candidate focal species for targeted repair.
#' @param threshold Numeric scalar giving the abundance multiplier threshold used
#'   to count species as satisfied. The default of \code{1.10} corresponds to
#'   species reaching +10\%.
#' @param beam_width Integer; number of partial repair branches retained at each
#'   layer of the beam search.
#' @param max_depth Integer; maximum number of sequential cell additions allowed
#'   in the beam search when attempting to repair the solution.
#' @param top_expand Integer; number of candidate cells considered when
#'   expanding each partial repair branch.
#' @param randomise_within_cost_prune Logical; if \code{TRUE}, ties among cells
#'   with equal removal cost during pruning are broken randomly.
#' @param fallback Character string specifying what to do if the targeted repair
#'   search fails. One of:
#'   \describe{
#'     \item{\code{"none"}}{Return \code{NULL}.}
#'     \item{\code{"full_solve"}}{Attempt a fallback call to
#'     \code{\link{solve_min_footprint_given_K}} using the current state as the
#'     incumbent.}
#'   }
#'
#' @return
#' Returns a repaired and pruned solution state satisfying \code{target_K}, if a
#' feasible solution is found. Returns \code{NULL} if no feasible solution is
#' found and \code{fallback = "none"}.
#'
#' @details
#' If the supplied \code{state} already satisfies both the geometric mean
#' constraint and \code{target_K}, the function simply prunes the solution and
#' returns it.
#'
#' Otherwise, the function:
#' \enumerate{
#'   \item identifies a set of candidate unsatisfied species closest to the
#'   threshold using \code{\link{candidate_species_for_next_K}},
#'   \item for each candidate species, attempts targeted repair using
#'   \code{\link{repair_state_for_species_beam}},
#'   \item compares the feasible repaired solutions by footprint size, and
#'   \item returns the smallest-footprint repaired state found.
#' }
#'
#' @export
advance_K_by_one <- function(state,
                             target_K,
                             GM_target,
                             sparse,
                             cost = NULL,
                             top_species = 5,
                             threshold = 1.10,
                             beam_width = 10,
                             max_depth = 4,
                             top_expand = 12,
                             overshoot_lambda = overshoot_lambda,
                             randomise_within_cost_prune = TRUE,
                             fallback = c("none", "full_solve"),
                             P_over_R = NULL) {
  
  # What to do if advancing by one can't find a solution
  fallback <- match.arg(fallback)
  
  # Number of species
  S <- length(state$m_s)
  
  # Set costs
  if (is.null(cost)) cost <- state$cost
  
  # If already feasible, just prune and return
  if ((state$gm_sumlog / S) >= GM_target && state$N10 >= target_K) {
    st <- prune_feasible_for_K(
      state, sparse,
      GM_target = GM_target,
      K = target_K,
      randomise_within_cost = randomise_within_cost_prune
    )
    return(st)
  }
  
  # Get candidate species for next K 
  # (these are the most likely fo their targets to be met)
  cand_species <- candidate_species_for_next_K(
    state = state,
    top_species = top_species,
    threshold = threshold
  )
  
  # If none then ...
  if (!length(cand_species)) return(NULL)
  
  # Tract state and footprint
  best_state <- NULL
  best_fp <- Inf
  
  # Just try to find solution targeting best candidate species
  for (s in cand_species) {
    
    #
    st_try <- repair_state_for_species_beam(
      state = state,
      focal_s = s,
      target_K = target_K,
      GM_target = GM_target,
      sparse = sparse,
      cost = cost,
      threshold = threshold,
      beam_width = beam_width,
      max_depth = max_depth,
      top_expand = top_expand,
      overshoot_lambda = overshoot_lambda,
      randomise_within_cost_prune = randomise_within_cost_prune
    )

    if (!is.null(st_try)) {
      fp <- sum(cost[st_try$x == 1L])
      if (fp < best_fp) {
        best_fp <- fp
        best_state <- st_try
      }
    }
  }
  
  if (!is.null(best_state)) return(best_state)
  
  # Fallback if local targeted repair fails
  if (fallback == "full_solve") {
    return(
      solve_min_footprint_given_K(
        P_over_R = P_over_R, # replace with actual object in outer wrapper
        beta = state$beta,
        GM_target = GM_target,
        K = target_K,
        cost = cost,
        randomise_within_cost_prune = randomise_within_cost_prune,
        incumbent = state
      )
    )
  }
  
  NULL
}

#' Identify candidate species to target when advancing to the next \eqn{K}
#'
#' Selects unsatisfied species that are closest to reaching the abundance
#' threshold, to guide targeted repair when moving from a solution for
#' \eqn{K} to a solution for \eqn{K+1}. Species are ranked by their shortfall
#' from the threshold, and the function returns the subset that are nearest to
#' being satisfied.
#'
#' This helper is intended for incremental optimisation, where the next
#' satisfied species is most likely to come from among those currently closest
#' to the threshold.
#'
#' @param state A solution state object containing species-level abundance
#'   multipliers in \code{state$m_s}.
#' @param top_species Integer; maximum number of candidate species to return.
#' @param threshold Numeric scalar giving the abundance multiplier threshold used
#'   to count species as satisfied. The default of \code{1.10} corresponds to
#'   species reaching +10\%.
#' @param slack_buffer Numeric scalar giving an optional tolerance around the
#'   smallest shortfall. Species whose shortfall is within
#'   \code{min(shortfall) + slack_buffer} are also considered as candidates
#'   before truncating to \code{top_species}.
#'
#' @return An integer vector of species indices corresponding to unsatisfied
#'   species closest to the threshold. Returns \code{integer(0)} if all species
#'   already satisfy the threshold.
#'
#' @details
#' The function first identifies species with \code{state$m_s < threshold}. For
#' these species, it computes the shortfall \code{threshold - state$m_s} and
#' ranks species from smallest to largest shortfall.
#'
#' Candidate species are then formed from:
#' \enumerate{
#'   \item the species with the smallest shortfalls, up to \code{top_species},
#'   and
#'   \item any additional species lying within \code{slack_buffer} of the best
#'   shortfall.
#' }
#'
#' The final output is truncated to at most \code{top_species} species.
#'
#' @seealso
#' \code{\link{advance_K_by_one}},
#' \code{\link{repair_state_for_species_beam}}
#'
#' @examples
#' \dontrun{
#' cand <- candidate_species_for_next_K(
#'   state = st_current,
#'   top_species = 5,
#'   threshold = 1.10,
#'   slack_buffer = 0.02
#' )
#'
#' cand
#' }
#'
#' @export
candidate_species_for_next_K <- function(state,
                                         top_species = 5,
                                         threshold = 1.10,
                                         slack_buffer = 0.02) {
  # Which species have unsatisfied targets
  unsat <- which(state$m_s < threshold)
  if (!length(unsat)) return(integer(0))
  
  # Compute shortfall
  shortfall <- threshold - state$m_s[unsat]
  
  # Primary candidates are those closest to threshold
  ord <- order(shortfall, decreasing = FALSE)
  cand <- unsat[ord]
  
  # Optional buffer to include all species within a small band of the best shortfall
  if (length(cand) > 1) {
    keep <- shortfall[ord] <= (min(shortfall) + slack_buffer)
    cand <- unique(c(cand[keep], cand[seq_len(min(top_species, length(cand)))]))
  } else {
    cand <- cand[1]
  }
  
  cand[seq_len(min(top_species, length(cand)))]
}


#' Score a candidate cell addition for a focal species
#'
#' Computes a heuristic score for adding a candidate cell to the current
#' solution state when targeting a particular unsatisfied species. The score
#' combines the cell's direct contribution to the focal species, its overlap
#' contribution to other near-threshold species, and its effect on the
#' geometric mean objective, scaled by cell cost.
#'
#' This helper is intended for use in targeted repair and beam-search expansion
#' when incrementally advancing from a solution for \eqn{K} to one for
#' \eqn{K+1}.
#'
#' @param state A solution state object containing the current selected cells
#'   and cached quantities required by \code{\link{delta_add}}.
#' @param c Integer index of the candidate cell to score.
#' @param focal_s Integer index of the focal species currently being targeted
#'   for repair.
#' @param sparse Sparse row representation of the candidate cell by species gain
#'   matrix, as returned by \code{\link{prep_sparse_rows}}.
#' @param near_species Optional integer vector of species indices representing
#'   additional near-threshold species that may benefit from the candidate cell.
#' @param cost Optional numeric vector of per-cell costs. If \code{NULL}, the
#'   function uses \code{state$cost}.
#' @param w_focal Numeric weight applied to the candidate cell's contribution to
#'   the focal species.
#' @param w_near Numeric weight applied to the candidate cell's contribution to
#'   the species listed in \code{near_species}.
#' @param w_gm Numeric weight applied to the candidate cell's change in the
#'   geometric mean summary, as returned by \code{\link{delta_add}}.
#'
#' @return A numeric scalar giving the heuristic score of adding the candidate
#'   cell. Larger values indicate more promising additions.
#'
#' @details
#' The score is calculated as:
#'
#' \deqn{
#' \frac{
#' w_{focal} \times \mathrm{focal\_gain} +
#' w_{near} \times \mathrm{near\_gain} +
#' w_{gm} \times \Delta gm
#' }{\mathrm{cost}_c}
#' }
#'
#' where:
#' \itemize{
#'   \item \code{focal_gain} is the summed contribution of cell \code{c} to the
#'   focal species,
#'   \item \code{near_gain} is the summed contribution of cell \code{c} to the
#'   species in \code{near_species}, and
#'   \item \code{d$d_gm} is the change in the geometric mean summary obtained
#'   from \code{\link{delta_add}}.
#' }
#'
#' This score is intended as a ranking heuristic rather than an exact objective,
#' and is used to prioritise which candidates to include in targeted repair
#' search.
#'
#' @seealso
#' \code{\link{delta_add}},
#' \code{\link{advance_K_by_one}},
#' \code{\link{repair_state_for_species_beam}}
#'
#' @examples
#' \dontrun{
#' sc <- score_add_cell_for_species(
#'   state = st_current,
#'   c = 42,
#'   focal_s = 7,
#'   sparse = sparse,
#'   near_species = c(3, 8, 11),
#'   cost = cost
#' )
#'
#' sc
#' }
#'
#' @export
score_add_cell_for_species <- function(state, c, focal_s, sparse,
                                       near_species = integer(0),
                                       cost = NULL,
                                       w_focal = 1.0,
                                       w_near = 0.25,
                                       w_gm = 1e-6) {
  idx <- sparse$idx
  wts <- sparse$wts
  if (is.null(cost)) cost <- state$cost
  
  ii <- idx[[c]]
  ww <- wts[[c]]
  
  focal_gain <- sum(ww[which(ii == focal_s)])
  
  near_gain <- 0
  if (length(near_species)) {
    near_gain <- sum(ww[ii %in% near_species])
  }
  
  d <- delta_add(state, c, idx, wts)
  
  (w_focal * focal_gain +
      w_near  * near_gain +
      w_gm    * d$d_gm) / cost[c]
}


#' Repair a solution for a focal species using beam search
#'
#' Attempts to construct a repaired solution state that satisfies a higher
#' species target by adding cells that benefit a specified focal species. The
#' search proceeds using a bounded beam search over partial sequences of cell
#' additions, retaining only the most promising partial repairs at each depth.
#'
#' This function is intended for incremental optimisation when advancing from a
#' solution for \eqn{K} to one for \eqn{K+1}, under the assumption that only a
#' small number of additional cells may be needed.
#'
#' @param state A solution state object containing the current selected cells and
#'   cached quantities such as \code{m_s}, \code{N10}, and \code{gm_sumlog}.
#' @param focal_s Integer index of the focal species being targeted for repair.
#' @param target_K Integer; target number of species required to exceed the
#'   abundance threshold.
#' @param GM_target Numeric scalar giving the minimum required geometric mean
#'   abundance gain across species.
#' @param sparse Sparse row representation of the candidate cell by species gain
#'   matrix, as returned by \code{\link{prep_sparse_rows}}.
#' @param cost Optional numeric vector of per-cell costs. If \code{NULL}, the
#'   function uses \code{state$cost}.
#' @param threshold Numeric scalar giving the abundance multiplier threshold used
#'   to count species as satisfied. The default of \code{1.10} corresponds to
#'   species reaching +10\%.
#' @param beam_width Integer; number of partial repair branches retained at each
#'   beam-search layer.
#' @param max_depth Integer; maximum number of sequential cell additions allowed
#'   in the search.
#' @param top_expand Integer; number of candidate cells considered when
#'   expanding each partial branch.
#' @param randomise_within_cost_prune Logical; if \code{TRUE}, ties among cells
#'   with equal removal cost during pruning are broken randomly.
#'
#' @return A repaired and pruned solution state satisfying \code{target_K}, if a
#'   feasible solution is found. Returns \code{NULL} if no feasible repair is
#'   found within the specified beam width and search depth.
#'
#' @details
#' The function first restricts attention to currently unselected cells that
#' contribute to the focal species. It then performs a beam search in which each
#' branch represents a partial sequence of added cells.
#'
#' At each depth:
#' \enumerate{
#'   \item each current branch is expanded by scoring candidate next-cell
#'   additions using \code{\link{score_add_cell_for_species}},
#'   \item only the top \code{top_expand} candidate cells are considered for each
#'   branch,
#'   \item each resulting child branch is either:
#'   \itemize{
#'     \item pruned and stored if it already satisfies both the geometric mean
#'     constraint and \code{target_K}, or
#'     \item retained as a partial repair with an internal branch score if it is
#'     not yet feasible,
#'   }
#'   \item the best \code{beam_width} partial repairs are retained for the next
#'   depth.
#' }
#'
#' Branch scoring is based on the remaining shortfall of the focal species, the
#' remaining shortfall of other near-threshold species, and a very small penalty
#' for the number of added cells. Lower branch scores are preferred.
#'
#' If one or more feasible repaired states are found at a given depth, the
#' smallest-footprint pruned state is returned immediately.
#'
#' This function does not attempt arbitrary global search; it is a targeted local
#' repair heuristic intended to exploit the approximate nesting of solutions
#' across increasing \code{K}.
#'
#' @seealso
#' \code{\link{score_add_cell_for_species}},
#' \code{\link{advance_K_by_one}},
#' \code{\link{prune_feasible_for_K}},
#' \code{\link{apply_add}}
#'
#' @examples
#' \dontrun{
#' st_repaired <- repair_state_for_species_beam(
#'   state = st_current,
#'   focal_s = 12,
#'   target_K = st_current$N10 + 1L,
#'   GM_target = 1.02,
#'   sparse = sparse,
#'   cost = cost,
#'   threshold = 1.10,
#'   beam_width = 10,
#'   max_depth = 4,
#'   top_expand = 12,
#'   randomise_within_cost_prune = TRUE
#' )
#'
#' if (!is.null(st_repaired)) {
#'   sum(cost[st_repaired$x == 1L])
#'   st_repaired$N10
#' }
#' }
#'
#' @export
repair_state_for_species_beam <- function(state,
                                          focal_s,
                                          target_K,
                                          GM_target,
                                          sparse,
                                          cost = NULL,
                                          threshold = 1.10,
                                          beam_width = 10,
                                          max_depth = 4,
                                          top_expand = 12,
                                          overshoot_lambda = 0,
                                          randomise_within_cost_prune = TRUE) {
  
  # Species and weights from sparse
  idx <- sparse$idx
  wts <- sparse$wts
  
  # Number of species
  S <- length(state$m_s)
  
  # Set costs
  if (is.null(cost)) cost <- state$cost
  
  # Candidates not selected
  unselected0 <- which(state$x == 0L)
  
  # And, of those, which include the focal species
  cand0 <- unselected0[vapply(unselected0, function(c) focal_s %in% idx[[c]], logical(1))]
  
  # If none then skip
  if (!length(cand0)) return(NULL)
  
  # Which species are close to threshold
  near_species <- which(state$m_s < threshold & state$m_s >= (threshold - 0.05))
  
  # Function to compute branch score
  branch_score <- function(st, added_cells) {
    shortfall_focal <- max(0, threshold - st$m_s[focal_s])
    near_shortfall <- sum(pmax(0, threshold - st$m_s[near_species]))
    # lower is better
    shortfall_focal + 0.2 * near_shortfall + 1e-6 * sum(cost[added_cells])
  }
  
  # Function to compute expand score
  expand_score <- function(st, c) {
    score_add_cell_for_species(
      state = st, c = c, focal_s = focal_s,
      sparse = sparse, near_species = near_species, cost = cost
    )
  }
  
  # Each beam element stores current state and added cells
  beam <- list(list(state = state, added = integer(0)))
  
  # Track states
  feasible_states <- list()
  feasible_score <- numeric(0)
  
  for (depth in 1:max_depth) {
  
    # Track children
    children <- list()
    
    for (b in beam) {
      
      # State and added
      st <- b$state
      added <- b$added
      
      # Unselected cells
      unselected <- which(st$x == 0L)
      
      # And of those have the species
      cand <- unselected[vapply(unselected, function(c)
        focal_s %in% idx[[c]], logical(1))]
      cand <- setdiff(cand, added)
      
      # Skip candidate
      if (!length(cand)) {
        next
      }
      
      # Compute scores for expansion
      sc <- vapply(cand, function(c)
        expand_score(st, c), numeric(1))
      ord <- order(sc, decreasing = TRUE)
      
      # Identify candidates
      cand <- cand[ord[seq_len(min(top_expand, length(cand)))]]
      
      c = cand[1]
      # Evalute candidates
      for (c in cand) {
        
        st2 <- apply_add(st, c, idx, wts)
        added2 <- c(added, c)
        
        # if feasible for target_K, prune and store
        if ((st2$gm_sumlog / S) >= GM_target && st2$N10 >= target_K) {
          stp <- prune_feasible_for_K(
            st2,
            sparse,
            GM_target = GM_target,
            K = target_K,
            randomise_within_cost = randomise_within_cost_prune
          )
          
          if ((stp$gm_sumlog / S) >= GM_target && stp$N10 >= target_K) {
            feasible_states[[length(feasible_states) + 1L]] <- stp
            feasible_score <- c(
              feasible_score,
              sum(cost[stp$x == 1L]) +
                overshoot_lambda * max(0, stp$N10 - target_K)
            )
          }
        } else {
          children[[length(children) + 1L]] <- list(
            state = st2,
            added = added2,
            score = branch_score(st2, added2)
          )
        }
      }
    }
    
    # If we found feasible repaired states, return the best one
    if (length(feasible_states)) {
      j <- which.min(feasible_score)
      return(feasible_states[[j]])
    }
    
    if (!length(children)) break
    
    child_scores <- vapply(children, `[[`, numeric(1), "score")
    keep <- order(child_scores, decreasing = FALSE)[seq_len(min(beam_width, length(children)))]
    beam <- children[keep]
  }
  
  NULL
}
