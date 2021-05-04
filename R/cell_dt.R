#' @title Cell data table without cytoplasm mapping
#'
#' @description A \code{data.table} containing the information for plotting the cell structure.
#'
#' @format A \code{data.table} with 6 columns, which are:
#' \describe{
#' \item{x}{the x coordinate of the polygon}
#' \item{y}{the y coordinate of the polygon}
#' \item{pol}{index of the polygon}
#' \item{organ}{a string that identifies the subcellular structure}
#' \item{color}{the value to fill the polygon in an example}
#' \item{comb}{a string that uniquely identify the subpolygon}
#' }
"cell_dt"
