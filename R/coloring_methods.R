compute_enrichment <- function(genes, plot_data, gene_loc_table, universe_set, categorical_classes=categorical_classes){
    localization_values <- data.table()
    # for each subcellular component computes FDR
    for (t in unique(plot_data$subcell_struct)){
        list_interest <- intersect(gene_loc_table[subcell_struct == t]$gene_symbol, universe_set)
        my_list <- genes

        LY <- length(intersect(my_list, list_interest))
        RY <- length(setdiff(list_interest, my_list))
        LN <- length(setdiff(my_list, list_interest))
        RN <- length(universe_set) - (LY + RY + LN)

        test <- fisher.test(matrix(c(LY, RY, LN, RN), 2, 2), alternative="greater")

        row <- data.table(subcell_struct = t, enr = test$estimate, pval = test$p.value)
        localization_values <- rbind(localization_values, row)
    }

    localization_values <- localization_values[order(pval, decreasing = FALSE)]
    localization_values <- localization_values[, `:=` ("value"=categorical_classes[pval > categorical_classes$start &
                                                                                       pval <= categorical_classes$end, values],
                                                       "color_grad"=categorical_classes[pval > categorical_classes$start &
                                                                                           pval <= categorical_classes$end, colors]), by=subcell_struct]

    final_dt <- merge.data.table(plot_data,
                                 localization_values[, c("subcell_struct", "color_grad", "value")],
                                 by.x = "subcell_struct",
                                 by.y = "subcell_struct")

    final_dt[, comb := factor(comb, levels = intersect(plot_data$comb, final_dt$comb))
             ][, value := factor(value, levels = rev(unique(localization_values$value)))]

    final_dt <- final_dt[order(comb)]

    return(list("final_dt"=final_dt, "localization_values"=localization_values))
}

#' Plot the neuron figure coloring subcellular localizations according to FDR
#' from the Fisher’s test, used to assess the statistical significance of the
#' enrichment.
#'
#' @description This function creates the plot of the neuron figure with colored
#'   subcellular localizations according to FDR values.
#' @param genes A character vector of gene names.
#' @param plot_data A \code{data.table} with the polygon coordinates to be
#'   plotted.
#' @param pictograph A character string with the name of the pictograph to be used.
#' @param gene_loc_table A \code{data.table} with information for mapping genes
#'   to subcellular localizations.
#' @param coloring_mode If "enrichment", computes the False Discovery Rate
#'   grouping genes by subcellular localization.
#' @param categorical_classes An optional \code{data.table} with values for
#'   specified intervals and associated colors. Default classification is
#'   divided in 8 intervals.
#'
#' @return \itemize{A list with two \code{data.table}s:
#'   \item{"final_dt"}{contains all the necessary information for plotting
#'   subcellular localizations with color values depending on the FDR;}
#'   \item{"localization_dt"}{contains a table with the enrichment scores and
#'   pvalues obtained for each subcellular localization} }
#'
#' @examples
#' #compute_enrichment(genes, plot_data, gene_loc_table)
#'
#' @import data.table
#' @export
assign_color_by_fdr <- function(genes, plot_data, pictograph, gene_loc_table, coloring_mode, categorical_classes=NULL){

    universe_set <- unique(gene_loc_table$gene_symbol)
    genes <- intersect(genes$gene_symbol, universe_set)

    res <- compute_enrichment(genes=genes,
                              plot_data,
                              universe_set = universe_set,
                              gene_loc_table=gene_loc_table,
                              categorical_classes)
    final_dt <- res$final_dt
    localization_values <- res$localization_values

    colors_vec <- categorical_classes$colors
    colors.spe <- colors_vec[sort(as.numeric(unique(as.vector(final_dt$value))), decreasing = TRUE)]

    lab <- as.expression(sapply(categorical_classes$lab, function(x) x))
    lab.spe <- lab[sort(as.numeric(unique(as.vector(final_dt$value))), decreasing = TRUE)]

    ecmx <- max(final_dt[subcell_struct == "extracellular_region"]$x)-(max(final_dt[subcell_struct == "extracellular_region"]$x)-min(final_dt[subcell_struct == "extracellular_region"]$x))/3

    ecmy <- (min(final_dt[subcell_struct == "extracellular_region"]$y)+max(final_dt[subcell_struct == "extracellular_region"]$y))/3

    bs=25
    p <- ggplot(final_dt, aes(x, y, color=color_grad, fill=value)) +
        scale_fill_manual(values = colors.spe, name="FDR", labels=lab.spe) +
        scale_color_manual(values = rep("black", length(unique(final_dt$subcell_struct)))) +
        scale_size_manual(values = rep(0.005, length(final_dt[, first(color_grad), by=subcell_struct]$V1))) +
        geom_polygon(aes(subgroup=comb)) +
        scale_y_reverse() +
        guides(color = FALSE) +
        theme_void()  +
        theme(legend.title = element_text(size=bs*0.9),
              legend.text = element_text(size=bs*0.9))

    p <- p + annotate("text", x=ecmx, y=ecmy, label="ECM", size=0.2*bs)


    return(list("plot"=p,
                "lab_title"="FDR",
                "localization_values"=localization_values,
                "final_dt"=final_dt,
                "ranges" = categorical_classes))
}

#' Assign colors to subcellular structures according to fold changes of genes
#'
#' @description This function colors subcellular structures based on mean (or
#'   median) of fold change of genes.
#' @param genes A \code{data.table} of gene names with associated log fold
#'   change values. Columns must be named "gene_symbol" and "logFC".
#' @param plot_data A \code{data.table} with the polygon coordinates to be
#'   plotted.
#' @param pictograph A character string with the name of the pictograph to be used.
#' @param gene_loc_table A \code{data.table} with information for mapping genes
#'   to subcellular localizations.
#' @param col_name  A character string with the name of the column on which the
#'   user wants to base the color of cellular localizations when "median" or
#'   "mean" are the chosen coloring method.
#' @param categorical_classes An optional \code{data.table} with values for
#'   specified intervals and associated colors. Default classification is
#'   divided in 8 intervals.
#' @param coloring_mode Either "mean" or "median". It computes the mean (or
#'   median) of log fold change values grouping genes by subcellular
#'   localization. Default is "mean".
#' @param together A boolean specifying whether or not genes are considered all
#'   together without any classification. Default value is false, so genes will
#'   be grouped by specified categories.
#' @return A list containing a ggplot2 object and two \code{data.table}s.
#'   "final_dt" contains data associated to the ggplot2 object, while
#'   "localization_values" contains the resulting mean (or median) fold changes
#'   for each subcellular localization.
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
assign_color_by_value <- function(genes, plot_data, pictograph, gene_loc_table, col_name, categorical_classes, coloring_mode="mean", together=FALSE){

    gene_loc_table <- gene_loc_table[subcell_struct %in% unique(plot_data$subcell_struct)]

    genes_sel <- merge.data.table(genes,
                                  gene_loc_table[, c("gene_symbol", "subcell_struct")],
                                  by = "gene_symbol")

    if (nrow(genes_sel)==0){
        localization_values <- data.table(subcell_struct = "NA")
    }

    if (coloring_mode == "mean"){
        localization_values <- genes_sel[, .(mean(get(col_name), na.rm = TRUE)), by=subcell_struct]

        if (all(localization_values$V1 > 0) || all(localization_values$V1 < 0)){
            localization_values <- localization_values[order(abs(V1), decreasing = TRUE)]
        } else {
            localization_values <- localization_values[order(V1, decreasing = FALSE)]
        }

        setnames(localization_values, old="V1", new=eval(coloring_mode))
    } else {
        if (coloring_mode == "median"){
            localization_values <- genes_sel[, .(median(get(col_name), na.rm = TRUE)), by=subcell_struct]

            if (all(localization_values$V1 > 0) || all(localization_values$V1 < 0)){
                localization_values <- localization_values[order(abs(V1), decreasing = TRUE)]
            } else {
                localization_values <- localization_values[order(V1, decreasing = FALSE)]
            }

            setnames(localization_values, old="V1", new=eval(coloring_mode))
        } else {
            cat("Unrecognized parameter, assigned default coloring_mode=\"mean\"")
        }
    }

    if (!together){

        f <- function(vec) {
            if(length(.categorical_classes <- which(abs(vec) > categorical_classes$start & abs(vec) <= categorical_classes$end))) .categorical_classes else NA
        }

        if (nrow(localization_values)!=0){
            localization_values$value <- categorical_classes$values[mapply(f, localization_values[, get(coloring_mode)])]
            localization_values$color_grad <- categorical_classes$colors[mapply(f, localization_values[, get(coloring_mode)])]
        } else {
            localization_values$value <- "500"
            localization_values$color_grad <- "grey90"
        }

        # localization_values <- localization_values[, `:=` ("value"=categorical_classes[(get(coloring_mode)) > abs(categorical_classes[[class]][, start]) &
        #                                                                                    (get(coloring_mode)) <= abs(categorical_classes[[class]][, end]), values],
        #                                                    "color_grad"=categorical_classes[(get(coloring_mode)) > abs(categorical_classes[[class]][, start]) &
        #                                                                                         (get(coloring_mode)) <= abs(categorical_classes[[class]][, end]), colors]), by=subcell_struct]

    } else {

        f <- function(vec) {
            if(length(.categorical_classes <- which(vec > categorical_classes$start & vec <= categorical_classes$end))) .categorical_classes else NA
        }

        if (nrow(localization_values)!=0){
            localization_values$value <- categorical_classes$values[mapply(f, localization_values[, get(coloring_mode)])]
            localization_values$color_grad <- categorical_classes$colors[mapply(f, localization_values[, get(coloring_mode)])]
        }
    }

    final_dt <- merge.data.table(plot_data,
                                 localization_values,
                                 by.x = "subcell_struct",
                                 by.y = "subcell_struct",
                                 all = TRUE)

    if (nrow(localization_values)!=0){
        final_dt[, comb := factor(comb, levels = intersect(plot_data$comb, final_dt$comb))
                ][, value := factor(value, levels = unique(localization_values$value))]
    } else {
        final_dt[, value  := 500
                 ][, value := factor(value)
                   ][, color_grad := "grey90"]
    }

    final_dt <- final_dt[order(comb)]

    if (all(localization_values[, get(coloring_mode)] > 0) || all(localization_values[, get(coloring_mode)] < 0)){
        mixed <- FALSE
        dec = TRUE
    } else {
        mixed <- TRUE
        dec = FALSE
    }

    lab_title <- paste0(col_name, " ", coloring_mode)

    colors_vec <- categorical_classes$colors
    colors.spe <- colors_vec[sort(as.numeric(unique(as.vector(final_dt$value))), decreasing = dec)]

    lab <- as.expression(sapply(categorical_classes$lab, function(x) x))
    lab.spe <- lab[sort(as.numeric(unique(as.vector(final_dt$value))), decreasing = dec)]

    na_val <- max(as.numeric(final_dt$value), na.rm = TRUE)
    final_dt <- final_dt[is.na(value), value := as.factor(na_val + 1)
                         ][is.na(color_grad), color_grad := "grey90"]

    if (length(levels(plot_data$subcell_struct)) > length(localization_values$subcell_struct)){
        colors.spe <- c(colors.spe,  "grey90")
    }

    nogreysquares <- copy(final_dt[color_grad != "grey90"])
    nogreysquares <- nogreysquares[, value := factor(value, levels = unique(localization_values$value))]

    ecmx <- max(final_dt[subcell_struct == "extracellular_region"]$x)-(max(final_dt[subcell_struct == "extracellular_region"]$x)-min(final_dt[subcell_struct == "extracellular_region"]$x))/3

    ecmy <- (min(final_dt[subcell_struct == "extracellular_region"]$y)+max(final_dt[subcell_struct == "extracellular_region"]$y))/3

    bs = 25
    p <- ggplot(final_dt, aes(x, y, color=value, fill=value)) +
        scale_fill_manual(values = colors.spe,
                          name = lab_title,
                          labels = lab.spe,
                          breaks = levels(nogreysquares$value)) +
        scale_color_manual(values = rep("black", length(unique(final_dt$subcell_struct))),
                           name = lab_title,
                           labels = lab.spe,
                           breaks = levels(nogreysquares$value)) +
        scale_size_manual(values = rep(0.005, length(final_dt[, first(color_grad), by=subcell_struct]$V1))) +
        geom_polygon(aes(subgroup=comb)) +
        scale_y_reverse() +
        #guides(color = FALSE) +
        theme_void() +
        theme(legend.title = element_text(size=bs*0.9),
              legend.text = element_text(size=bs*0.9))

    p

    if (pictograph == "neuron"){
        p <- p + annotate("text", x=ecmx, y=ecmy, label="ECM", size=0.2*bs)
    }

    return(list("plot"=p,
                "lab_title"=lab_title,
                "localization_values"=localization_values,
                "final_dt"=final_dt,
                "ranges" = categorical_classes))
}



