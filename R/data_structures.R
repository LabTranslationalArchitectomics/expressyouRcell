#' @title Neuron data table
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

#' @title Cell data table
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

#' @title Fibroblast data table
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
"fibroblast_dt"

#' @title Microglia data table
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
"microglia_dt"

#' @title Lymphocyte data table
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
"lymphocyte_dt"

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

#' @title Example of gene_localization table created from a mus musculus genome annotation (GENCODE v.22)
#'
#' @description An example gene_localization table.
#'
#' @format A \code{data.table}s with 3 columns, which are:
#' \describe{
#' \item{ID}{the gene ontology term identifier}
#' \item{subcell_struct}{the subcellular structure name}
#' \item{gene_symbol}{the gene name}
#' }
"gene_loc_table_mm22"
