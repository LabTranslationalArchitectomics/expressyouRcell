compute_enrichment <- function(genes, plot_data, gene_loc_table, universe_set, categorical_classes=categorical_classes){
    localization_values <- data.table()
    # for each subcellular component computes FDR
    for (t in unique(plot_data$subcell_struct)){
        list_interest <- intersect(gene_loc_table[Description == t]$gene_symbol, universe_set)
        my_list <- genes

        LY <- length(intersect(my_list, list_interest))
        RY <- length(setdiff(list_interest, my_list))
        LN <- length(setdiff(my_list, list_interest))
        RN <- length(universe_set) - (LY + RY + LN)

        test <- fisher.test(matrix(c(LY, RY, LN, RN), 2, 2), alternative="greater")

        row <- data.table(Description = t, enr = test$estimate, pval = test$p.value)
        localization_values <- rbind(localization_values, row)
    }

    localization_values <- localization_values[order(pval, decreasing = FALSE)]
    localization_values <- localization_values[, `:=` ("value"=categorical_classes[pval > categorical_classes$start &
                                                                                       pval <= categorical_classes$end, values],
                                                       "enr_color"=categorical_classes[pval > categorical_classes$start &
                                                                                           pval <= categorical_classes$end, colors]), by=Description]

    final_dt <- merge.data.table(plot_data,
                                 localization_values[, c("Description", "enr_color", "value")],
                                 by.x = "subcell_struct",
                                 by.y = "Description")

    final_dt[, comb := factor(comb, levels = intersect(plot_data$comb, final_dt$comb))
    ][, value := factor(value, levels = rev(unique(localization_values$value)))]

    final_dt <- final_dt[order(comb)]

    return(list("final_dt"=final_dt, "localization_values"=localization_values))
}

#' Plot the neuron figure coloring subcellular localizations according to FDR
#'
#' @description This function creates the plot of the neuron figure with colored
#'   subcellular localizations according to FDR values.
#' @param gene_set A character vector of gene names.
#' @param plot_data A \code{data.table} with the polygon coordinates to be
#'   plotted.
#' @param gene_loc_table A \code{data.table} with information for mapping genes
#'   to subcellular localizations.
#' @param categorical_classes An optional \code{data.table} with values for
#'   specified intervals and associated colors. Default classification is
#'   divided in 8 intervals.
#'
#' @return \itemize{A list with two \code{data.table}s:
#' \item{"final_dt"}{contains all the necessary information for plotting
#' subcellular localizations with color values depending on the FDR;}
#' \item{"localization_dt"}{contains a table with the enrichment scores and
#' pvalues obtained for each subcellular localization} }
#'
#' @examples
#' #compute_enrichment(genes, plot_data, gene_loc_table)
#'
#' @import data.table
#' @export
assign_color_by_fdr <- function(gene_set, plot_data, gene_loc_table, categorical_classes=NULL){

    universe_set <- unique(gene_loc_table$gene_symbol)
    gene_set <- intersect(gene_set, universe_set)

    if (is.null(categorical_classes)){
        categorical_classes <- data.table(start = c(0, 1*10^-30, 1*10^-20, 1*10^-10, 1*10^-5, 1*10^-4, 1*10^-3, 5*10^-2),
                                          end = c(1*10^-30, 1*10^-20, 1*10^-10, 1*10^-5, 1*10^-4, 1*10^-3, 5*10^-2, 1),
                                          values = seq(1:8),
                                          colors = c("#31407c", "#4b5a95", "#737ead", "#296982", "#3484a3", "#5c9cb5", "#adcdda", "grey90")
        )[, lab := paste("<", .SD[, end]), by=values]

        categorical_classes <- categorical_classes[, lab := paste("<", .SD[, end]), by=values]
    }

    res <- compute_enrichment(genes=gene_set,
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

    bs=25
    p <- ggplot(final_dt, aes(x, y, color=enr_color, fill=value)) +
        scale_fill_manual(values = colors.spe, name="FDR", labels=lab.spe) +
        scale_color_manual(values = rep("black", length(unique(final_dt$subcell_struct)))) +
        scale_size_manual(values = rep(0.005, length(final_dt[, first(enr_color), by=subcell_struct]$V1))) +
        geom_polygon(aes(subgroup=comb)) +
        scale_y_reverse() +
        guides(color = FALSE) +
        theme_void()  +
        theme(legend.text=element_text(size=bs*0.7),
              legend.title=element_text(size=bs*0.7))

    return(list("plot"=p, "localization_values"=localization_values, "final_dt"=final_dt))
}

#' Assign colors to subcellular structures according to fold changes of genes
#'
#' @description This function colors subcellular structures based on mean (or median) of fold change of genes.
#' @param genes A \code{data.table} of gene names with associated log fold change values. Columns must be named "gene_symbol" and "logFC".
#' @param plot_data A \code{data.table} with the polygon coordinates to be plotted.
#' @param colors A character vector of hex color codes for palette generation.
#' @param gene_loc_table A \code{data.table} with information for mapping genes to subcellular localizations.
#' @param coloring_mode Either "mean" or "median". It computes the mean (or median) of log fold change values grouping genes
#' by subcellular localization. Default is "mean".
#' @return A list containing a ggplot2 object and two \code{data.table}s. "final_dt" contains data associated to the
#' ggplot2 object, while "localization_values" contains the resulting mean (or median) fold changes for each subcellular localization.
#' @examples
#' #assign_color_by_fc(genes=genes_down,plot_data = neuron_dt,gene_loc_table = gene_loc_table,colors=c("darkred", "white"),value="median")
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
assign_color_by_value <- function(genes, plot_data, gene_loc_table, col_name, categorical_classes, coloring_mode="mean", together=FALSE){

    gene_loc_table <- gene_loc_table[Description %in% unique(plot_data$subcell_struct)]

    genes_sel <- merge.data.table(genes,
                                  gene_loc_table[, c("gene_symbol", "Description")],
                                  by = "gene_symbol")

    if (coloring_mode == "mean"){
        localization_values <- genes_sel[, .(mean(get(col_name), na.rm = TRUE)), by=Description]
        localization_values <- localization_values[order(V1)]
        setnames(localization_values, old="V1", new=eval(coloring_mode))
    } else {
        if (coloring_mode == "median"){
            localization_values <- genes_sel[, .(median(get(col_name), na.rm = TRUE)), by=Description]
            localization_values <- localization_values[order(V1)]
            setnames(localization_values, old="V1", new=eval(coloring_mode))
        } else {
            cat("Unrecognized parameter, assigned default coloring_mode=\"mean\"")
        }
    }

    if (!together){

        f <- function(vec) {
            if(length(.categorical_classes <- which(abs(vec) > categorical_classes$start & abs(vec) <= categorical_classes$end))) .categorical_classes else NA
        }


        localization_values$value <- categorical_classes$values[mapply(f, localization_values[, get(coloring_mode)])]
        localization_values$color_grad <- categorical_classes$colors[mapply(f, localization_values[, get(coloring_mode)])]



        # localization_values <- localization_values[, `:=` ("value"=categorical_classes[(get(coloring_mode)) > abs(categorical_classes[[class]][, start]) &
        #                                                                                    (get(coloring_mode)) <= abs(categorical_classes[[class]][, end]), values],
        #                                                    "color_grad"=categorical_classes[(get(coloring_mode)) > abs(categorical_classes[[class]][, start]) &
        #                                                                                         (get(coloring_mode)) <= abs(categorical_classes[[class]][, end]), colors]), by=Description]

    } else {

        f <- function(vec) {
            if(length(.categorical_classes <- which(vec > categorical_classes$start & vec <= categorical_classes$end))) .categorical_classes else NA
        }

        localization_values$value <- categorical_classes$values[mapply(f, localization_values[, get(coloring_mode)])]
        localization_values$color_grad <- categorical_classes$colors[mapply(f, localization_values[, get(coloring_mode)])]
    }


    final_dt <- merge.data.table(plot_data,
                                 localization_values,
                                 by.x = "subcell_struct",
                                 by.y = "Description")

    final_dt[, comb := factor(comb, levels = intersect(plot_data$comb, final_dt$comb))
             ][, value := factor(value, levels = unique(localization_values$value))]

    final_dt <- final_dt[order(comb)]

    if (together){
        lab_title <- paste0(col_name, coloring_mode)
        dec = FALSE
    } else {
        lab_title <- paste0(col_name, " ", coloring_mode)
        dec = TRUE
    }

    colors_vec <- categorical_classes$colors
    colors.spe <- colors_vec[sort(as.numeric(unique(as.vector(final_dt$value))), decreasing = dec)]

    lab <- as.expression(sapply(categorical_classes$lab, function(x) x))
    lab.spe <- lab[sort(as.numeric(unique(as.vector(final_dt$value))), decreasing = dec)]

    final_dt[is.na(value), value := max(as.numeric(final_dt$value))+1
             ][is.na(color_grad), color_grad := "grey90"]

    if (length(levels(final_dt$value)) > length(colors.spe)){
        colors.spe <- c(colors.spe,  "grey90")
    }


    nogreysquares <- copy(final_dt[color_grad != "grey90"])
    nogreysquares <- nogreysquares[, value := factor(value, levels = unique(localization_values$value))]

    p <- ggplot(final_dt, aes(x, y, color=color_grad, fill=value)) +
        scale_fill_manual(values = colors.spe,
                          name = lab_title,
                          labels = lab.spe,
                          breaks = levels(nogreysquares$value)) +
        scale_color_manual(values = rep("black", length(unique(final_dt$subcell_struct)))) +
        scale_size_manual(values = rep(0.005, length(final_dt[, first(color_grad), by=subcell_struct]$V1))) +
        geom_polygon(aes(subgroup=comb)) +
        scale_y_reverse() +
        guides(color = FALSE) +
        theme_void()

    if (together){
        p <- p +
            theme(legend.background = element_rect(fill="lightgrey",
                                                   size=0.5,
                                                   linetype="solid"))
    }

    p

    return(list("localization_values"=localization_values,
                "plot"=p))
}



