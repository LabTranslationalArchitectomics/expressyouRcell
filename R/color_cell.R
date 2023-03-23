#' Create static cellular pictographs
#'
#' @description This function creates static cellular pictographs with the chosen coloring method for assigning colors to
#'   subcellular localizations.
#' @param timepoint_list A list of \code{data.tables}. Each must have at least a column named "gene_symbol". Optionally,
#'   the list of \code{data.tables} can also contain a numerical column with i) its expression level, in terms of read
#'   counts, count per million of reads (CPM) or reads per kilobase of gene per million (RPKM); ii) fold-changes and
#'   p-values from upstream differential analyses.
#' @param pictograph A character string, or a vector with multiple character strings, with the names of the pictographs
#'   to be used.
#'   Currently available pictographs are named "cell", "neuron", "microglia", "fibroblast", "macrophage" and "lymphocyte".
#'   The corresponding \code{data.table} with the polygon coordinates is loaded automatically.
#'   Default value is generic cell pictograph ("cell").
#'   Paraemeter values can be in the form of (\code{pictograph="cell"}, or \code{pictograph=c("cell","neuron")}).
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
#' @param enr_color_ranges An optional vector of color shades for defining the significance of the enrichment analysis.
#' This parameter works with the "enrichment" option of the \code{coloring_mode} parameter. Example of usage: c("palegreen", "palegreen4"), or c("white", "darkblue").
#' @param thr A numeric value specifying the cut-off value to be applied on the \code{col_name} column.
#' @param pval_col A character with the name of the column containing the statistical significance values.
#' @param pval_thr A numeric value with the cutoff value to be applied on the \code{pval_col} column.
#' @param legend A boolean value for choosing to plot the legend or not. Default is false.
#' @param scaling A  boolean value for choosing whether to scale values by row, or not. Default is false.
#' @param scale_libsize A  boolean value for choosing whether to scale values by library size (cpm), or not. Default is false.
#'
#' @return A list containing four data structures. The first is the localization_values \code{data.table}, with six
#'   columns, which reports for each subcellular component: i) its name, ii) the numeric value computed during the
#'   colour assignment step, iiI) a numeric code for grouping the cellular localizations by colour, iv) its associated
#'   colour shade, the v) the identifier of #' each dataset, and vi) the grouping variable. The second data structure is
#'   the ranges \code{data.table}, summarising the information on the ranges (e.g. start, end, colour, and labels) used
#'   to categorise each subcellular localization. The third data structure is the plot list, containing the graphical
#'   objects of class \code{ggplot} with the resulting cellular pictographs. The fourth data structure is the list of
#'   final_dt \code{data.table}, with the datasets used to plot the resulting cellular pictographs.
#' When multiple cell types are assigned to the pictograph parameter, the function return a single ggplot object for the arranged multiple cell type in the \code{merged_cells} item.
#'
#' @import data.table
#' @import patchwork
#' @import rlist
#'
#' @export
#'
color_cell <- function(timepoint_list,
                       pictograph="cell",
                       gene_loc_table,
                       coloring_mode='enrichment',
                       col_name=NULL,
                       colors=NULL,
                       group_by=NULL,
                       grouping_vars=NULL,
                       enr_color_ranges=NULL,
                       thr=NULL,
                       pval_col=NULL,
                       pval_thr=NULL,
                       legend=FALSE,
                       scaling=FALSE,
                       scale_libsize=FALSE){

    if (suppressWarnings(!all(lapply(timepoint_list, function(x) inherits(x, "data.table"))))){
        dt <- lapply(timepoint_list, function(x) inherits(x, "data.table"))
        not_dt <- names(dt[dt==FALSE])
        timepoint_list <- lapply(timepoint_list, function(x) data.table(x))
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

    cell_diagram <- widths_list <- list()
    for (p in pictograph){

        if (all(p == "cell")){
            plot_data <- cell_dt
            w <- c(3,1)
        } else {
            if (all(p == "neuron")) {
                plot_data <- neuron_dt
                w <- c(4,1)
            } else {
                if (all(p == "fibroblast")) {
                    plot_data <- fibroblast_dt
                    w <- c(3,1)
                } else {
                    if (all(p == "microglia")) {
                        plot_data <- microglia_dt
                        w <- c(2,1)
                    } else {
                        if (all(p == "macrophage")) {
                            plot_data <- macrophage_dt
                            w <- c(3,1)
                        } else {
                            if (all(p == "lymphocyte")) {
                                plot_data <- lymphocyte_dt
                                w <- c(3,1)
                            } else {
                                stop(paste0("No available pictograph defined as:", p))
                            }
                        }
                    }
                }
            }
        }

        widths_list[[p]] <- w

        if (coloring_mode == 'mean' || coloring_mode == 'median'){

            if (all(unlist(lapply(timepoint_list, function(x) x[, get(col_name)]-floor(x[, get(col_name)])==0)))) {
                warning("Integer values have been detected in the input data. For more reliable results, we strongly suggest you to scale your data before proceeding with the pictographic visulization.")
            }

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

                    if (scaling==TRUE){
                        tobescaled_allcols <- list.cbind(lapply(timepoint_list, function(x) x[, get(col_name)]))
                        tobescaled_allcols <- data.frame(tobescaled_allcols)
                        rownames(tobescaled_allcols) <- timepoint_list[[1]]$gene_symbol
                        scaled_allcols <- scale(t(tobescaled_allcols), center = TRUE, scale = TRUE)

                        scaled_tp_list <- list()
                        for (tp in colnames(scaled_allcols)){
                            scaled_tp_list[[tp]] <- data.table(scaled_allcols, keep.rownames = "gene_symbol"
                            )[, c("gene_symbol", tp), with=FALSE]

                        }
                        timepoint_list <- scaled_tp_list
                    }

                    if (scale_libsize==TRUE){
                        tobescaled_allcols <- list.cbind(lapply(timepoint_list, function(x) x[, get(col_name)]))
                        tobescaled_allcols <- data.frame(tobescaled_allcols)
                        rownames(tobescaled_allcols) <- timepoint_list[[1]]$gene_symbol
                        scaled_allcols <- tobescaled_allcols/colSums(tobescaled_allcols)

                        scaled_tp_list <- list()
                        for (tp in colnames(scaled_allcols)){
                            scaled_tp_list[[tp]] <- data.table(scaled_allcols, keep.rownames = "gene_symbol"
                            )[, c("gene_symbol", tp), with=FALSE]

                        }
                        timepoint_list <- scaled_tp_list
                    }





                    locdt_l <- plot_l <- finaldt_l <- list()
                    for (tp in names(timepoint_list)){

                        cat(paste0("Creating cell pictograph for stage: ", tp, "\n"))

                        if (!is.null(group_by)){
                            # create a separate plot for each category
                            if (is.null(grouping_vars)) {
                                # all the categories are plotted
                                grouping_vars[[group_by]] <- as.character(unique(unlist(lapply(timepoint_list, function(x)
                                    unique(x[, get(group_by)])))))
                            }

                            grouped_out <- list()


                            if (is.null(enr_color_ranges)){
                                ranges <-  discrete_symmetric_ranges(timepoint_list = timepoint_list,
                                                                     plot_data = plot_data,
                                                                     gene_loc_table = gene_loc_table,
                                                                     col_name = col_name,
                                                                     grouping_vars = grouping_vars,
                                                                     colors = colors,
                                                                     coloring_mode = coloring_mode)
                            } else {
                                    ranges <- enr_color_ranges
                            }


                            for (v in unlist(grouping_vars)){
                                genes <- timepoint_list[[tp]][get(group_by) == v]

                                colored_out <- assign_color_by_value(genes = genes,
                                                                     plot_data = plot_data,
                                                                     pictograph = p,
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

                            if (is.null(enr_color_ranges)){
                                ranges <- discrete_symmetric_ranges(timepoint_list = timepoint_list,
                                                                    plot_data,
                                                                    gene_loc_table,
                                                                    col_name,
                                                                    grouping_vars,
                                                                    colors,
                                                                    coloring_mode,
                                                                    together=TRUE)
                            } else {
                                ranges <- enr_color_ranges
                            }

                            colored_out <- assign_color_by_value(genes = genes,
                                                                 plot_data = plot_data,
                                                                 pictograph = p,
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
                    output[["cell_type"]]<- p
                    output[["timepoint_list"]] <- timepoint_list

                    cell_diagram[[p]] <- output

                    #return(output)
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

                    cat(paste0("Creating cell pictograph for stage: ", tp, "\n"))

                    if (!is.null(group_by)){
                        # create a separate plot for each category
                        if (is.null(grouping_vars)) {
                            # all the categories are plotted
                            grouping_vars[[group_by]] <- as.character(unique(unlist(lapply(timepoint_list, function(x) unique(x[, get(group_by)])))))
                        }

                        grouped_out <- list()

                        for (v in unlist(grouping_vars)){
                            genes <- timepoint_list[[tp]][get(group_by) == v]

                            if (is.null(enr_color_ranges)){
                                default_ranges <- data.table(start = c(0, 1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2),
                                                             end = c(1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2, 1),
                                                             values = seq(1:5),
                                                             colors = c("#40486e", "#296982", "#3484a3", "#adcdda", "grey50")
                                                             )[, lab := paste("<", .SD[, end]), by=values]

                                enr_ranges <- default_ranges
                            } else {
                                enr_ranges <- enr_color_ranges
                                start <- enr_ranges[1]
                                end <- enr_ranges[2]

                                if (start=="white" | start=="#FFFFFF"){
                                    col <- colorRampPalette(c(start, end))(5)
                                    col <- col[-1]
                                } else {
                                    col <- colorRampPalette(c(start, end))(4)
                                }


                                ranges <- data.table(start = c(0, 1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2),
                                                     end = c(1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2, 1),
                                                     values = seq(1:5),
                                                     colors = c(rev(col), "grey50")
                                                     )[, lab := paste("<", .SD[, end]), by=values]

                                enr_ranges <- ranges
                            }

                            colored_out <- assign_color_by_fdr(genes = genes,
                                                               plot_data = plot_data,
                                                               pictograph=,
                                                               gene_loc_table = gene_loc_table,
                                                               categorical_classes = enr_ranges,
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

                        if (is.null(enr_color_ranges)){
                            default_ranges <- data.table(start = c(0, 1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2),
                                                         end = c(1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2, 1),
                                                         values = seq(1:5),
                                                         colors = c("#40486e", "#296982", "#3484a3", "#adcdda", "grey50")
                                                         )[, lab := paste("<", .SD[, end]), by=values]
                            enr_ranges <- default_ranges
                        } else {
                            start <- enr_color_ranges[1]
                            end <- enr_color_ranges[2]

                            if (start=="white" | start=="#FFFFFF"){
                                col <- colorRampPalette(c(start, end))(5)
                                col <- col[-1]
                            } else {
                                col <- colorRampPalette(c(start, end))(4)
                            }


                            ranges <- data.table(start = c(0, 1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2),
                                                 end = c(1*10^-10, 1*10^-5, 1*10^-4, 5*10^-2, 1),
                                                 values = seq(1:5),
                                                 colors = c(rev(col), "grey50")
                                                 )[, lab := paste("<", .SD[, end]), by=values]

                            enr_ranges <- ranges
                        }

                        colored_out <- assign_color_by_fdr(genes = genes,
                                                           plot_data = plot_data,
                                                           pictograph=,
                                                           gene_loc_table = gene_loc_table,
                                                           categorical_classes = enr_ranges,
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
                output[["cell_type"]] <- p
                output[["timepoint_list"]] <- timepoint_list

                cell_diagram[[p]] <- output

                #return(output)

                # assign_color_by_fdr(genes,
                #                     plot_data=plot_data,
                #                     gene_loc_table=gene_loc_table,
                #                     categorical_classes=NULL)

            } else {
                cat("Wrong coloring method")
            }

        }
    }

    if (length(pictograph)>1){
        merged_list <- list()
        timepoints <- names(cell_diagram[[1]]$plot)
        for (tp in timepoints){
            merged <- patchwork::wrap_plots(lapply(cell_diagram, function(x) x[["plot"]][[tp]]+coord_fixed()))
            merged_list[[paste0(tp, "_merged")]] <- merged
        }
        cell_diagram[["merged_cells"]] <- merged_list
        return(cell_diagram)
    } else {
        return(cell_diagram[[1]])
    }
}
