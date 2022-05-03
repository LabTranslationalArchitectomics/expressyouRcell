#' Create static cellular pictograms
#'
#' @description This function creates static cellular pictograms with the chosen coloring method for assigning colors to
#'   subcellular localizations.
#' @param timepoint_list A list of \code{data.tables}. Each must have at least a column named "gene_symbol". Optionally,
#'   the list of \code{data.tables} can also contain a numerical column with i) its expression level, in terms of read
#'   counts, count per million of reads (CPM) or reads per kilobase of gene per million (RPKM); ii) fold-changes and
#'   p-values from upstream differential analyses.
#' @param pictogram A character string with the name of the pictogram to be used.
#' The corresponding \code{data.table} with the polygon coordinates is loaded automaticaly. Default is generic cell pictogram.
#'   (\code{pictogram="cell"}, or \code{pictogram="neuron"}).
#' @param gene_loc_table A \code{data.table} with information for mapping genes to subcellular localizations. The
#'   localization of the genes can be either provided by the user or created through the \code{map_gene_localization}
#'   function.
#' @param coloring_mode Either "enrichment", "mean" or "median". Default is "enrichment". If "enrichment" is specified,
#'   gene set enrichment analysis restricted to the sub-tree ontology of cellular components is performed on the list of
#'   gene symbols entered input by the user. The color shade of each subcellular compartment is defined by the p-values
#'   from the Fisher's test, used to assess the statistical significance of the enrichment.
#'
#'   If "mean" or "median" are specified, a name of a numerical columns must also be provided as \code{col_name}
#'   parameter. This method computes the mean (or median) of values specified in the \code{col_name} column. The color
#'   of each cellular compartment is defined by the mean (or median) of the gene-specific values associated to each
#'   localization.
#' @param col_name A character string with the name of the column on which the user wants to base the color of cellular
#'   localizations when "median" or "mean" are the chosen coloring method. Default is null.
#' @param colors A list with named vectors of hex color codes for palette generation. The number of named vectors must
#'   match the vectors given as input with the \code{grouping_vars} parameters. Each vector should contain two color
#'   codes (for the starting and stopping point of the generated palettes). By default, the package automatically select
#'   the optimal palette based on the \code{coloring_mode} parameter
#' @param group_by A character string with the name of the categorical variable to group genes into separate subsets.
#'   The categorical variable name must be equal to a column name in the \code{data.tables} of the \code{timepoint_list}
#'   parameter.
#' @param grouping_vars A list with a named vector of characters for subselecting genes belonging only to the specified
#'   values of the variable \code{group_by}. The name of the vector must match the \code{group_by} parameter. The vector
#'   of characters must be present in the \code{data.tables} columns named as \code{group_by}.
#' @param ranges An optional \code{data.table} with values for specified intervals and associated colors. By default,
#'   the package automatically creates the most appropriate binning for the input data.
#' @param thr A numeric value specifying the cut-off value to be applied on the \code{col_name} column.
#' @param pval_col A character with the name of the column containing the statistical significance values.
#' @param pval_thr A numeric value with the cutoff value to be applied on the \code{pval_col} column.
#' @param legend A boolean value for choosing to plot the legend or not.
#' @return A list containing four data structures. The first is the localization_values \code{data.table}, with six
#'   columns, which reports for each subcellular component: i) its name, ii) the numeric value computed during the
#'   colour assignment step, iiI) a numeric code for grouping the cellular localizations by colour, iv) its associated
#'   colour shade, the v) the identifier of #' each dataset, and vi) the grouping variable. The second data structure is
#'   the ranges \code{data.table}, summarising the information on the ranges (e.g. start, end, colour, and labels) used
#'   to categorise each subcellular localization. The third data structure is the plot list, containing the graphical
#'   objects of class \code{ggplot} with the resulting cellular pictograms. The fourth data structure is the list of
#'   final_dt \code{data.table}, with the datasets used to plot the resulting cellular pictograms.
#'
#'
#' @import data.table
#'
#' @export
#'
color_cell <- function(timepoint_list,
                       pictogram="cell",
                       gene_loc_table,
                       coloring_mode='enrichment',
                       col_name=NULL,
                       colors=NULL,
                       group_by=NULL,
                       grouping_vars=NULL,
                       ranges=NULL,
                       thr=NULL,
                       pval_col=NULL,
                       pval_thr=NULL,
                       legend=FALSE){

    if (all(pictogram == "cell")){
        plot_data <- cell_dt
        w <- c(3,1)
    } else {
        if (all(pictogram == "neuron")) {
            plot_data <- neuron_dt
            w <- c(4,1)
        } else {
            if (all(pictogram == "fibroblast")) {
                plot_data <- fibroblast_dt
                w <- c(3,1)
            } else {
                if (all(pictogram == "microglia")) {
                    plot_data <- microglia_dt
                    w <- c(2,1)
                } else {
                    stop("No available pictogram with this name")
                }
            }
        }
    }


    if (!inherits(timepoint_list, "list")){
        tmp_timepoint_list <- list()
        tmp_timepoint_list[["1"]] <- timepoint_list
        timepoint_list <- tmp_timepoint_list
        #stop("Input data must be organized as a list of data.table(s)")
    }

    if (!all(unlist(lapply(timepoint_list, function(x) "gene_symbol" %in% colnames(x))))){
        stop("Input data must include a column of gene names named precisely \"gene_symbol\".")
    }

    if (coloring_mode == 'mean' || coloring_mode == 'median'){
        if (is.null(col_name)){
            stop("ERROR: col_name parameter is missing")
            #cat("ERROR: col_name parameter is missing")
        } else {

            if (!all(unlist(lapply(timepoint_list, function(x) col_name %in% colnames(x))))){
                stop("Check that every list contains a column named as specified in col_name param")
            } else {
                if (!is.null(group_by)){
                    # we are grouping by class of DEGs
                    if (!all(unlist(lapply(timepoint_list, function(x) group_by %in% colnames(x))))) {

                        if (is.null(thr) || is.null(pval_col) || is.null(pval_thr)){
                            stop("You want to group your data by the \"group_by\" parameter, but some data.table does not contain a column with the name you specified and I don't know how to create this column, since some parameters are missing (check thr, pval_col, pval_thr)")
                        } else {
                            timepoint_list <- lapply(timepoint_list, function(x) x[, eval(group_by) := "="
                                                                 ][get(col_name) <= -thr & get(pval_col) <= pval_thr, eval(group_by) := '-'
                                                                   ][get(col_name) >= thr & get(pval_col) <= pval_thr, eval(group_by) := '+'])

                        }


                    } else {
                        # random colors are chosen for each category in group_by column
                        if (is.null(colors)){
                            all_grouping_vars <- as.character(unique(unlist(lapply(timepoint_list, function(x) unique(x[, get(group_by)])))))

                            colors=list()
                            color_codes <- sample_colors(length(all_grouping_vars))
                            for (v in seq_len(length(all_grouping_vars))){
                                colors[[all_grouping_vars[v]]] <- c("white", color_codes[[v]])
                            }
                        }

                        if (!is.null(grouping_vars)){
                            colors <- colors[names(colors) %in% unlist(grouping_vars)]
                        }
                    }
                } else {
                    # not grouping by any categorical variable
                    if (is.null(colors)){
                        colors = c("white", "#296d98")
                    }
                }

                locdt_l <- plot_l <- finaldt_l <- list()
                for (tp in names(timepoint_list)){

                    cat(paste0("Creating cell pictogram for stage: ", tp, "\n"))

                    if (!is.null(group_by)){
                        # create a separate plot for each category
                        if (is.null(grouping_vars)) {
                            # all the categories are plotted
                            grouping_vars[[group_by]] <- as.character(unique(unlist(lapply(timepoint_list, function(x)
                                unique(x[, get(group_by)])))))
                        }

                        grouped_out <- list()

                        ranges <-  discrete_symmetric_ranges(timepoint_list = timepoint_list,
                                                             plot_data = plot_data,
                                                             gene_loc_table = gene_loc_table,
                                                             col_name = col_name,
                                                             grouping_vars = grouping_vars,
                                                             colors = colors,
                                                             coloring_mode = coloring_mode)
                        for (v in unlist(grouping_vars)){
                            genes <- timepoint_list[[tp]][get(group_by) == v]

                            colored_out <- assign_color_by_value(genes = genes,
                                                                 plot_data = plot_data,
                                                                 pictogram = pictogram,
                                                                 gene_loc_table = gene_loc_table,
                                                                 col_name = col_name,
                                                                 categorical_classes = ranges[get(group_by) == v],
                                                                 coloring_mode = coloring_mode)

                            colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp
                                                                                                         ][, eval(group_by) := v]

                            locdt_l[[paste0(tp,v)]] <- colored_out[["localization_values"]]
                            plot_l[[paste0(tp, v)]] <- colored_out[["plot"]]
                            finaldt_l[[paste0(tp, v)]] <- colored_out[["final_dt"]]


                        }
                    } else {
                        # create a plot regardless the classification
                        if (!is.null(grouping_vars)) {
                            # only specified categories are plotted, and genes are averaged regardless any classification

                            tmp <- data.table::data.table(gene_symbol=character(),
                                                          x=character(),
                                                          y=numeric())

                            setnames(tmp, old=c("x", "y"), new = c(names(grouping_vars), col_name))

                            for (v in unlist(grouping_vars)){
                                genes <- timepoint_list[[tp]][get(names(grouping_vars)) == v,
                                                              c("gene_symbol", names(grouping_vars), col_name),
                                                              with = FALSE]
                                tmp <- funion(tmp, genes)
                            }

                            genes <- tmp

                        } else {
                            # genes belonging to all categories are plotted
                            genes <- timepoint_list[[tp]]
                        }

                        ranges <- discrete_symmetric_ranges(timepoint_list = timepoint_list,
                                                            plot_data,
                                                            gene_loc_table,
                                                            col_name,
                                                            grouping_vars,
                                                            colors,
                                                            coloring_mode,
                                                            together=TRUE)

                        colored_out <- assign_color_by_value(genes = genes,
                                                             plot_data = plot_data,
                                                             pictogram = pictogram,
                                                             gene_loc_table = gene_loc_table,
                                                             col_name = col_name,
                                                             categorical_classes = ranges,
                                                             coloring_mode = coloring_mode,
                                                             together = TRUE)


                        colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp]
                        locdt_l[[tp]] <- colored_out[["localization_values"]]

                        plot_l[[paste0(tp)]] <- colored_out[["plot"]]
                        finaldt_l[[paste0(tp)]] <- colored_out[["final_dt"]]

                    }
                }

                output <- list()

                output[["localization_values"]] <- rbindlist(locdt_l)
                output[["ranges"]] <- ranges
                output[["legend_label"]] <- colored_out[["lab_title"]]

                if (legend == TRUE){
                    output[["plot"]]<- lapply(plot_l,
                                              function(x) ggarrange(x,
                                                                    plot_legend_organelles(dt_legend = legend_dt),
                                                                    widths = w))
                } else {
                    output[["plot"]]<- plot_l
                }

                output[["final_dt"]]<- finaldt_l
                output[["cell_type"]]<- pictogram
                output[["timepoint_list"]] <- timepoint_list

                return(output)
            }
        }
    } else {
        if (coloring_mode == "enrichment"){
            if (!is.null(group_by)){
                # we are grouping by class of DEGs
                if (!all(unlist(lapply(timepoint_list, function(x) group_by %in% colnames(x))))) {
                    if (is.null(thr) || is.null(pval_col) || is.null(pval_thr)){
                        stop("You want to group your data by the \"group_by\" parameter, but some data.table does not contain a column with the name you specified and I don't know how to create this column, since some parameters are missing (check thr, pval_col, pval_thr)")
                    } else {
                        timepoint_list <- lapply(timepoint_list, function(x) x[, eval(group_by) := "="
                        ][as.numeric(get(col_name)) <= -thr & as.numeric(get(pval_col)) <= pval_thr, eval(group_by) := '-'
                        ][as.numeric(get(col_name)) >= thr & as.numeric(get(pval_col)) <= pval_thr, eval(group_by) := '+'])

                    }
                } else {
                    if (is.null(colors)){
                        # random colors are chosen for each category in group_by column
                        all_grouping_vars <- as.character(unique(unlist(lapply(timepoint_list, function(x) unique(x[, get(group_by)])))))

                        colors=list()
                        color_codes <- sample_colors(length(all_grouping_vars))
                        for (v in seq_len(length(all_grouping_vars))){
                            colors[[all_grouping_vars[v]]] <- c("white", color_codes[[v]])
                        }
                    }

                    if (!is.null(grouping_vars)){
                        colors <- colors[names(colors) %in% unlist(grouping_vars)]
                    }
                }
            } else {
                # not grouping by any categorical variable
                if (is.null(colors)){
                    colors = c("white", "#296d98")
                }
            }

            locdt_l <- plot_l <- finaldt_l <- list()
            for (tp in names(timepoint_list)){

                cat(paste0("Creating cell pictogram for stage: ", tp, "\n"))

                if (!is.null(group_by)){
                    # create a separate plot for each category
                    if (is.null(grouping_vars)) {
                        # all the categories are plotted
                        grouping_vars[[group_by]] <- as.character(unique(unlist(lapply(timepoint_list, function(x) unique(x[, get(group_by)])))))
                    }

                    grouped_out <- list()

                    for (v in unlist(grouping_vars)){
                        genes <- timepoint_list[[tp]][get(group_by) == v]

                        colored_out <- assign_color_by_fdr(genes = genes,
                                                           plot_data = plot_data,
                                                           pictogram=pictogram,
                                                           gene_loc_table = gene_loc_table,
                                                           categorical_classes = NULL,
                                                           coloring_mode = coloring_mode)

                        colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp
                                                                                                     ][, eval(group_by) := v]
                        locdt_l[[paste0(tp,v)]] <- colored_out[["localization_values"]]
                        plot_l[[paste0(tp, v)]] <- colored_out[["plot"]]
                        finaldt_l[[paste0(tp, v)]] <- colored_out[["final_dt"]]

                    }



                } else {
                    # create a plot regardless the classification
                    if (!is.null(grouping_vars)) {
                        # only specified categories are plotted, and genes are averaged regardless any classification

                        tmp <- data.table::data.table(gene_symbol=character())

                        for (v in unlist(grouping_vars)){
                            genes <- timepoint_list[[tp]][get(names(grouping_vars)) == v,
                                                          c("gene_symbol"),
                                                          with = FALSE]
                            tmp <- funion(tmp, genes)
                        }


                        genes <- tmp

                    } else {
                        # genes belonging to all categories are plotted
                        genes <- timepoint_list[[tp]]
                    }

                    colored_out <- assign_color_by_fdr(genes = genes,
                                                       plot_data = plot_data,
                                                       pictogram=pictogram,
                                                       gene_loc_table = gene_loc_table,
                                                       categorical_classes = NULL,
                                                       coloring_mode = coloring_mode)
                    colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp]
                    locdt_l[[tp]] <- colored_out[["localization_values"]]

                    plot_l[[paste0(tp)]] <- colored_out[["plot"]]
                    finaldt_l[[paste0(tp)]] <- colored_out[["final_dt"]]
                }
            }
            output <- list()

            output[["localization_values"]] <- rbindlist(locdt_l)
            output[["ranges"]] <- colored_out$ranges
            output[["legend_label"]] <- colored_out[["lab_title"]]
            if (legend == TRUE){
                output[["plot"]]<- lapply(plot_l,
                                          function(x) ggarrange(x,
                                                                plot_legend_organelles(dt_legend = legend_dt),
                                                                widths = w))
            } else {
                output[["plot"]]<- plot_l
            }
            output[["final_dt"]]<- finaldt_l
            output[["cell_type"]]<- pictogram
            output[["timepoint_list"]] <- timepoint_list

            return(output)

            # assign_color_by_fdr(genes,
            #                     plot_data=plot_data,
            #                     gene_loc_table=gene_loc_table,
            #                     categorical_classes=NULL)

        } else {
            cat("Wrong coloring method")
        }

    }

}
