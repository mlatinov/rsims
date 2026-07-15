#' Create IDs for Nested Hierarchical Structures
#'
#' @description
#' Generates integer identifiers for arbitrary nested structures.
#'
#' This helper is designed for hierarchical simulations where lower-level
#' observations belong to higher-level groups.
#'
#' Examples of nested structures:
#'
#' * students nested within schools
#' * patients nested within wards within hospitals
#' * employees nested within departments within companies
#'
#' @param levels A named numeric vector or list describing the number of
#' observations at each level.
#'
#' The order should go from highest level to lowest level.
#'
#' Example:
#'
#' \code{
#' list(
#'   hospital = 20,
#'   ward = 4,
#'   patient = 30
#' )
#' }
#'
#' creates:
#'
#' \eqn{20 \times 4 \times 30 = 2400}
#' observations.
#'
#' @return
#' A data frame containing integer IDs for every level.
#'
#' @examples
#'
#' ids <- make_nested_ids(
#'   list(
#'     hospital = 2,
#'     ward = 3,
#'     patient = 5
#'   )
#' )
#'
#' head(ids)
#'
#' @keywords internal
make_nested_ids <- function(levels){

  if(is.null(names(levels))){
    stop("levels must be a named vector or list.")
  }

  levels <- as.list(levels)

  n <- prod(unlist(levels))

  ids <- data.frame(
    row_id = seq_len(n)
  )

  total_levels <- length(levels)

  for(i in seq_len(total_levels)){

    current_level <- names(levels)[i]

    repeat_times <- prod(unlist(levels[(i + 1):total_levels]))

    if(is.na(repeat_times)){
      repeat_times <- 1
    }

    values <- rep(
      seq_len(levels[[i]]),
      each = repeat_times
    )

    values <- rep(
      values,
      length.out = n
    )

    ids[[paste0(current_level, "_id")]] <- values
  }

  ids$row_id <- NULL

  ids
}



#' Create IDs for Crossed Designs
#'
#' @description
#' Generates identifiers for crossed experimental designs where every
#' combination of groups appears.
#'
#' Examples:
#'
#' * restaurants crossed with cuisines
#' * subjects crossed with measurement occasions
#' * machines crossed with operators
#'
#' @param factors A named numeric vector or list describing the number of
#' levels for each crossed factor.
#'
#' @return
#' A data frame containing all combinations of factor IDs.
#'
#' @examples
#'
#' ids <- make_crossed_ids(
#'   list(
#'     restaurant = 5,
#'     cuisine = 10
#'   )
#' )
#'
#' head(ids)
#'
#' @keywords internal
make_crossed_ids <- function(factors){

  if(is.null(names(factors))){
    stop("factors must be a named vector or list.")
  }

  factors <- as.list(factors)

  grid <- expand.grid(
    lapply(
      factors,
      seq_len
    )
  )

  names(grid) <- paste0(
    names(factors),
    "_id"
  )

  rownames(grid) <- NULL

  grid
}