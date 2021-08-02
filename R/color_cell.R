#' Main function that creates static cellular pictogram with the desired
#' coloring method
#'
#' @description This function creates static cellular pictogram with the desired
#'   coloring method for assigning colors to subcellular localizations
#' @param timepoint_list A list of \code{data.table}s, one for each time point.
#'   Each one must have at least a column named "gene_symbol".
#' @param plot_data A \code{data.table} with the polygon coordinates to be
#'   plotted.
#' @param gene_loc_table A \code{data.table} with information for mapping genes
#'   to subcellular localizations.
#' @param coloring_mode Either "enrichment", "mean" or "median". Default is
#'   "enrichment". If "enrichment" is specified, enrichment analysis restricted
#'   to the sub-ontology of cellular components is performed on genes input by
#'   the user. Colors of each subcellular compartment are based on pvalues from
#'   the Fisher's test, used to assess the statistical significance of the
#'   enrichment.
#'
#'   If "mean" or "median" are specified, a name of a numerical columns must
#'   also be provided as \code{col_name} parameter. This method computes the
#'   mean (or median) of values specified in the \code{col_name} column.
#' @param col_name A character string with the name of the column on which the
#'   user wants to base the color of cellular localizations when "median" or
#'   "mean" are the chosen coloring method.
#' @param colors A character vector of hex color codes for palette generation.
#' @param group_by A character string with the name of the column with the
#'   categorical variable based on which the user can compute separated analysis
#'   on subset of genes grouped accordingly to the value of this variable.
#' @param grouping_vars A vector of characters for subselecting genes belonging
#'   only to the specified values of the varible \code{group_by}
#' @param ranges An optional \code{data.table} with values for specified
#'   intervals and associated colors.
#' @param thr a numeric value specifying the cutoff value to be applied on the
#'   col_name column
#' @param pval_col a character with the name of the column containing the
#'   statistical significance values
#' @param pval_thr a numeric value with the cutoff value to be applied on the
#'   pval_col column
#' @return
#'
#' @import data.table
#'
#' @export
#'
color_cell <- function(timepoint_list,
                       plot_data="cell",
                       gene_loc_table,
                       coloring_mode='enrichment',
                       col_name=NULL,
                       colors=NULL,
                       group_by=NULL,
                       grouping_vars=NULL,
                       ranges=NULL,
                       thr=NULL,
                       pval_col=NULL,
                       pval_thr=NULL){

    if (plot_data == "cell"){
        plot_data <- cell_dt
    } else {
        if (plot_data == "neuron") {
            plot_data <- neuron_dt_nocyto
        } else {
            stop("No available pictogram with this name")
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
                stop("check that every list contains a column named as specified in col_name param")
            } else {
                if (!is.null(group_by)){
                    # we are grouping by class of DEGs
                    if (!all(unlist(lapply(timepoint_list, function(x) group_by %in% colnames(x))))) {

                        if (is.null(thr) || is.null(pval_col) || is.null(pval_thr)){
                            stop("You want to group your data by the \"group_by\" parameter, but some data.table does not contain a column with the name you specified and I don't know how to create this column, since some parameters are missing (check thr, pval_col, pval_thr)")
                        } else {
                            timepoint_list <- lapply(timepoint_list, function(x) x[, eval(group_by) := "="
                                                                 ][get(col_name) <= thr & get(pval_col) <= pval_thr, eval(group_by) := '-'
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

                locdt_l <- plot_l <- list()
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
                                                                 gene_loc_table = gene_loc_table,
                                                                 col_name = col_name,
                                                                 categorical_classes = ranges[get(group_by) == v],
                                                                 coloring_mode = coloring_mode)

                            colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp
                                                                                                         ][, eval(group_by) := v]

                            locdt_l[[paste0(tp,v)]] <- colored_out[["localization_values"]]
                            plot_l[[paste0("plot_", tp, v)]] <- colored_out[["plot"]]


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

                        ranges <- discrete_symmetric_ranges(timepoint_list,
                                                            plot_data,
                                                            gene_loc_table,
                                                            col_name,
                                                            grouping_vars,
                                                            colors,
                                                            coloring_mode,
                                                            together=TRUE)

                        colored_out <- assign_color_by_value(genes = genes,
                                                              plot_data = plot_data,
                                                              gene_loc_table = gene_loc_table,
                                                              col_name = col_name,
                                                              categorical_classes = ranges,
                                                              coloring_mode = coloring_mode,
                                                              together = TRUE)


                        colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp]
                        locdt_l[[tp]] <- colored_out[["localization_values"]]

                        plot_l[[paste0("plot_", tp)]] <- colored_out[["plot"]]

                    }
                }

                output <- list()

                output[["localization_values"]] <- rbindlist(locdt_l)
                output[["ranges"]] <- ranges
                output[["plot"]]<- plot_l

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
                        ][get(col_name) <= thr & get(pval_col) <= pval_thr, eval(group_by) := '-'
                        ][get(col_name) >= thr & get(pval_col) <= pval_thr, eval(group_by) := '+'])

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

            locdt_l <- plot_l <- list()
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
                                                                gene_loc_table = gene_loc_table,
                                                                categorical_classes = NULL,
                                                                coloring_mode = coloring_mode)

                        colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp
                                                                                                     ][, eval(group_by) := v]
                        locdt_l[[paste0(tp,v)]] <- colored_out[["localization_values"]]
                        plot_l[[paste0("plot_", tp, v)]] <- colored_out[["plot"]]

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
                                                        gene_loc_table = gene_loc_table,
                                                        categorical_classes = NULL,
                                                        coloring_mode = coloring_mode)
                    colored_out[["localization_values"]] <- colored_out[["localization_values"]][, time_point := tp]
                    locdt_l[[tp]] <- colored_out[["localization_values"]]

                    plot_l[[paste0("plot_", tp)]] <- colored_out[["plot"]]
                }
            }
            output <- list()

            output[["localization_values"]] <- rbindlist(locdt_l)
            output[["ranges"]] <- ranges
            output[["plot"]]<- plot_l

            return(output)

            assign_color_by_fdr(genes,
                                plot_data=plot_data,
                                gene_loc_table=gene_loc_table,
                                categorical_classes=NULL)

        } else {
            cat("Wrong coloring method")
        }

    }

}
