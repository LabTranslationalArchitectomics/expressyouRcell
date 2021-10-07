#' @title Neuron data table without cytoplasm mapping
#'
#' @description A \code{data.table} containing the information for plotting the neuronal structure.
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
"neuron_dt"

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

#' @title Legend for neuron plot data table
#'
#' @description A \code{data.table} containing the information for plotting the legend for organelles.
#'
#' @format A \code{data.table} with 6 columns, which are:
#' \describe{
#' \item{x}{the x coordinate of the polygon}
#' \item{y}{the y coordinate of the polygon}
#' \item{pol}{index of the polygon}
#' \item{subcell_struct}{a string that identifies the subcellular structure}
#' \item{comb}{a string that uniquely identify the subpolygon}
#' }
"legend_dt"

#' @title Example of input time point list
#'
#' @description An example list of two \code{data.table}s, one for each time point.
#' Each data.table must have at least a column of gene names named precisely “gene_symbol”.
#'
#' @format A list of \code{data.table}s with 6 columns, which are:
#' \describe{
#' \item{gene_symbol}{the x coordinate of the polygon}
#' \item{logFC}{the y coordinate of the polygon}
#' \item{sma}{index of the polygon}
#' \item{wt}{a string that identifies the subcellular structure}
#' \item{pval}{a string that uniquely identify the subpolygon}
#' \item{class}{a string that uniquely identify the subpolygon}
#' }
"example_list"
