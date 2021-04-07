#' Plot legend of organelles in the cell
#'
#' @description This function plots the legend
#' @param dt_legend A \code{data.table} with the polygon coordinates to be plotted in the legend figure
#' @return a ggplot object with the legend of labeled organelles
#' @examples
#' #plot_legend(dt_legend=legend)
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#'
#' @export
#'
plot_legend_organelles <- function(dt_legend){

  text_labels <- data.table(names=unique(dt_legend$subcell_struct),
                            x = max(dt_legend$x)+10,
                            y = (dt_legend[, max(.SD[, y]), by=subcell_struct]$V1 +
                                dt_legend[, min(.SD[, y]), by=subcell_struct]$V1)/2)

  text_labels <- text_labels[, names := stringr::str_replace(names, "\\_", "\n")]

  bs=25
  p <- ggplot(dt_legend, aes(x, y, color=comb)) +
    scale_color_manual(values = rep("black", length(dt_legend[, first(subcell_struct), by=comb]$V1))) +
    geom_polygon(fill=NA) +
    annotate("text",
              x=text_labels$x,
              y=text_labels$y,
              label=text_labels$names,
              hjust = 0,
              size=0.6*bs) +
    coord_cartesian(xlim = c(90,400)) +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none")
  return(p)
}

#' Plot cellular structure with default colors
#'
#' @description This function creates the plot of the specified cellular structure.
#' @param dt_complete A \code{data.table} with the polygon coordinates to be plotted.
#' @return a ggplot object with the cell figure with default colors.
#' @examples
#' #plot_cell(dt_complete)
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#'
#' @export
plot_cell <- function(dt_complete){
  p <- ggplot(dt_complete, aes(x, y, fill=comb, color=comb)) +
    scale_fill_manual(values = dt_complete[, first(color), by=comb]$V1) +
    scale_color_manual(values = rep("black", length(dt_complete[, first(color), by=comb]$V1))) +
    geom_polygon() +
    scale_y_reverse() +
    theme_void() +
    theme(legend.position = "none")

  return(p)
}

#' Create table for mapping genes to subcellular localization
#'
#' @description This function generates the \code{data.table} necessary to associate genes to subcellular localization.
#' @param gtf_path A character vector with the gtf file path.
#' @return A \code{data.table} with the gene ontology subcellular localization term for each gene.
#' @examples
#' #map_gene_localization(ensembl_gene_set)
#'
#' @import data.table
#' @import clusterProfiler
#' @import org.Mm.eg.db
#' @import DOSE
#'
#' @details
#' A gene annotatation file, in GTF format is required as input. On this complete
#' set of gene symbols, a gene ontology enrichment analysis is performed to associate
#' a gene with a term in the cellular component ontology. For this purpose, only the
#' sub-ontology of the cellular components is taken into consideration. This
#' step generates the gene-localization data.table, which maps each gene to the
#' locations in the cellular structures, either cellular compartments or
#' macromolecular complexes.
#'
#' @export
#'
map_gene_localization <- function(gtf_path, dataSource = NA, organism = NA){

  gtf <- rtracklayer::import(gtf_path)

  annotation_gene_names <- unique(gtf$gene_name)

  ensembl_entrez <- bitr(annotation_gene_names, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb=org.Mm.eg.db, drop=T)
  cc_complete <- as.data.table(as.data.frame(enrichGO(gene = ensembl_entrez[,2],
                               OrgDb=org.Mm.eg.db,
                               ont = "CC",
                               pvalueCutoff = 1,
                               pAdjustMethod = "none",
                               qvalueCutoff = 1,
                               minGSSize = 1,
                               maxGSSize = length(ensembl_entrez$ENTREZID),
                               readable = T)))

  cc_complete_dt <- cc_complete[, .(gene_symbol = unlist(tstrsplit(gene_symbol, "\\/", type.convert = TRUE))), by = c("ID", "Description")]
  cc_complete_dt <- cc_complete_dt[, Description := stringr::str_replace(Description, "\\ ", "\\_")]
  setnames(cc_complete_dt, old = "gene_symbol", new="gene_symbol")
  return(cc_complete_dt)
}

#' Create animation between stages
#'
#' @description This function creates animation between stages
#' @param sample_name A character vector with the name of the sample
#' @param colors A vector of colors to be assigned to the categorical classes
#' @param lab A vector of labels to be assigned to the categorical classes
#' @param output_folder A file path where to save output animated pictures
#'
#' @examples
#' #create_animation(sample_name, colors, lab, output_folder)
#'
#' @import data.table
#' @import ggpubr
#' @import magick
#'
#' @export
create_animation <- function(sample_name, colors, lab, plots_dt, n_frames, fps, output_folder=NULL){

  if (is.null(output_folder)){
    output_folder = getwd()
  }

  if (!dir.exists(file.path(output_folder, "gif"))){
    dir.create(path = file.path(output_folder, "gif"), recursive = T)
  }

  legend <- create_legend(color_vector = colors, lab_vector = lab)

  anim_stages <- list()
  for (stage in names(samples[[sample_name]])){
    cat(paste0("stage ", stage, "\n"))
    p <- ggarrange(plots_dt[[samples[[sample_name]][[stage]]]] +
                     guides(fill = FALSE) +
                     ggtitle(stage),
                   legend,
                   widths = c(3, 1))
    ggsave(p, filename = file.path(output_folder, "gif", paste0(stage, ".png")), device = "png", width = 9, height = 3)

    anim_stages[[stage]] <- file.path(output_folder, "gif", paste0(stage, ".png"))

    im <- image_scale(image_read(file.path(output_folder, "gif", paste0(stage, ".png"))))
    anim_stages[[stage]] <- im
    unlink(file.path(output_folder, "gif", paste0(stage, ".png")))
  }

  start.time <- Sys.time()
  animation <- image_resize(image_join(anim_stages), '1800x600!') %>%
    image_morph(frames = n_frames) %>%
    image_animate(optimize = TRUE, fps = fps)
  end.time <- Sys.time()

  time.taken <- end.time - start.time
  cat(paste0("animation time:", time.taken, "\n"))

  start.time <- Sys.time()
  image_write(animation,
              file.path(output_folder, "gif", paste0(sample_name, ".gif")),
              quality = 100)
  end.time <- Sys.time()

  time.taken <- end.time - start.time
  cat(paste0("saving time:", time.taken, "\n"))

  cat("gif saved in: ", file.path(output_folder, "gif", paste0(sample_name, ".gif")))
}

#' Create a legend for the animated pictures
#'
#' @description This function creates a legend for the animated gif files
#' @param color_vector A vector of hex color codes
#' @param lab_vector A vector of labels for the legend
#' @return a gtable with the legend to be added on the animated pictures
#'
#' @import ggplot2
#'
create_legend <- function(color_vector, lab_vector){

  gb <- data.table(x=seq(1:length(color_vector)),
             y=seq(1:length(color_vector)))

  # create legend by combing stages
  static_plot <- ggplot(data = gb, aes(x=x, y=y, fill=as.factor(x))) +
    scale_fill_manual(values = rev(color_vector), name="FDR", labels=rev(lab_vector)) +
    geom_bar(stat = "identity")
    #guides(color = FALSE) +
    #theme(legend.text = element_text(margin = margin(r = 15, l= 4, unit = "pt"))) # + theme_void()

  # extract legend
  tmp <- ggplot_gtable(ggplot_build(static_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#' Compute discrete and symmetric ranges into which classify FC values
#'
#' @description This function creates a table with ranges for mapping fold change values with associated colors and labels.
#' @param timepoint_list A list of \code{data.table}s, one for each time point. Each one must have at least two defined columns named respectively
#' "gene_symbol" and "logFC".
#' @param grouping_vars A character vector with classes of genes to be considered.
#' @param colors A character vector with two color codes for generating a color palettes with a shade for each interval.
#' @param value A character vector with two color codes for generating a color palettes composed of n_intervals colors.
#' @param together An boolean value specifying whether genes with different classification labels should be considered together (regardless their classification) or separately
#' @return A list of \code{data.table}s (one for each classes specified in the parameter classes) with categorical classes into which the fold change values grouped for each subcellular localization are mapped.
#' @examples
#' #assign_color_by_pval(genes=genes_down,plot_data = neuron_dt,gene_loc_table = gene_loc_table,colors=c("darkred", "white"))
#'
#' @import data.table
#'
discrete_symmetric_ranges <- function(timepoint_list,
                                      plot_data,
                                      gene_loc_table,
                                      col_name,
                                      grouping_vars,
                                      colors,
                                      coloring_mode,
                                      together=FALSE){

  if (together){
    groupedval <-  c()
    for (tp in names(timepoint_list)){

      if (!is.null(grouping_vars)){

        tmp <- data.table::data.table(gene_symbol=character(),
                                      gene_id=character(),
                                      logFC=numeric(),
                                      pvalue=numeric(),
                                      class=character(),
                                      Stime=character())

        for (c in grouping_vars){
          genes <- timepoint_list[[tp]][class == c]
          tmp <- funion(tmp, genes)
        }
      } else {
        genes <- timepoint_list[[tp]]
      }


      groupedval <- c(groupedval,
                     groupval_byloc(genes = genes,
                                   plot_data,
                                   gene_loc_table,
                                   col_name,
                                   coloring_mode))
    }

    sup <- max(abs(min(groupedval)), max(groupedval))

    fixed_ranges_dt <- data.table(start = head(c(0, 0.5, seq(1, ceiling(sup))), -1),
                                  end = c(0.5, seq(1, ceiling(sup))))
  } else {
    max_v <- min_v <- widths <- c()
    for (c in grouping_vars){
      groupedval <-  c()
      for (tp in names(timepoint_list)){

        genes <- timepoint_list[[tp]][class == c]

        groupedval <- c(groupedval,
                       groupval_byloc(genes = genes,
                                      plot_data,
                                      gene_loc_table,
                                      col_name,
                                      coloring_mode))
      }

      max_v[[c]] <- ceiling(max(abs(groupedval), na.rm = TRUE))
      min_v[[c]] <- floor(min(abs(groupedval), na.rm = TRUE))

      widths[[c]] <- abs(max_v[[c]]-min_v[[c]])

    }

    max_range_width <- which(widths == max(unlist(widths)))

    if (max(unlist(widths)) <= 8){
      inf <- unlist(min_v[max_range_width])
      sup <- unlist(max_v[max_range_width])
      fixed_ranges_dt <- data.table(start = head(seq(inf, sup), -1),
                                    end = seq(inf, sup)[-1])

    } else {
      w <- max(unlist(widths))
      binsize <- w / 8

      sizes <- c(1, 1.5, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000, 2000)
      diff <- sizes - binsize

      binsize <- sizes[which(diff == min(diff[diff >0]))]

      inf <- unlist(min_v[max_range_width])
      sup <- w + inf


      fixed_ranges_dt <- data.table(start = head(seq(inf, sup, by = binsize), -1),
                                    end = seq(inf, sup, by = binsize)[-1])

    }

  }

  dtlist <- list()
  for (c in grouping_vars){
    #cat(c)
    colfunc <- colorRampPalette(colors[[c]])

    fixed_ranges_dt <- copy(fixed_ranges_dt)
    fixed_ranges_dt <- fixed_ranges_dt[, values := seq(1, nrow(fixed_ranges_dt))
                                       ][, lab := paste("<", .SD[, round(end, 2)]), by=values]

    colors_vector <- colfunc(n = (nrow(fixed_ranges_dt)-1)*3+1)
    colors_vector <- colors_vector[seq(1, length(colors_vector), len=nrow(fixed_ranges_dt))]
    fixed_ranges_dt <- fixed_ranges_dt[, colors := colors_vector]

    dtlist[[c]] <- fixed_ranges_dt
  }

  if (together){
    dtlist[["-"]] <- dtlist[["-"]][, end := -end
                                   ][, start := -start
                                     ][, lab := paste("<", .SD[, round(start, 2)]), by=values
                                       ][order(start)]

    dt_together <- rbind(dtlist$`-`, dtlist$`+`)

    dt_together[, values := seq(1, nrow(dt_together))]

    dtlist[['together']] <- dt_together
  }

  return(dtlist)

}

#' Compute mean or median of values associated to genes for each cellular compartment
#'
#' @param genes A \code{data.table} of gene names with associated log fold change values. Columns must be named "gene_symbol" and "logFC".
#' @param plot_data A \code{data.table} with the polygon coordinates to be
#'   plotted.
#' @param gene_loc_table A \code{data.table} with information for mapping genes
#'   to subcellular localizations.
#' @param col_name A character string with the name of the column on which the user wants
#' to base the color of cellular localizations when "median" or "mean" are the chosen coloring method.
#' @param coloring_mode Either "mean" or "median". Default is "mean".
#'
#' A name of a numerical columns must also be provided as \code{col_name} parameter.
#' This method computes the mean (or median) of values specified in the \code{col_name} column.
#'
#'
#' @import data.table
#'
groupval_byloc <- function(genes, plot_data, gene_loc_table, col_name, coloring_mode="mean"){

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

  return(localization_values[, get(coloring_mode)])

}


#' Sample random colors
#'
#' @param n An integer specifying the number of random color to sample
#'
#' @import RColorBrewer
#'
sample_colors <- function(n){
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  col=sample(col_vector, n)
  return(col)
}

#' Main function that creates static cellular pictogram with the desired coloring method
#'
#' @description This function creates static cellular pictogram with the desired coloring method for assigning colors to subcellular localizations
#' @param timepoint_list A list of \code{data.table}s, one for each time point. Each one must have at least two defined columns named respectively
#' "gene_symbol" and "logFC".
#' @param plot_data A \code{data.table} with the polygon coordinates to be
#'   plotted.
#' @param gene_loc_table A \code{data.table} with information for mapping genes
#'   to subcellular localizations.
#' @param coloring_mode Either "enrichment", "mean" or "median". Default is "enrichment".
#' If "enrichment" is specified, enrichment analysis restricted to the sub-ontology of cellular components
#' is performed on genes input by the user. Colors of each subcellular compartment
#' are based on pvalues from the Fisher’s test, used to assess the statistical
#' significance of the enrichment.
#'
#' If "mean" or "median" are specified, a name of a numerical columns must also be provided as \code{col_name} parameter.
#' This method computes the mean (or median) of values specified in the \code{col_name} column.
#' @param col_name A character string with the name of the column on which the user wants
#' to base the color of cellular localizations when "median" or "mean" are the chosen coloring method.
#' @param colors A character vector of hex color codes for palette generation.
#' @param group_by A character string with the name of the column with the categorical variable based on which the user
#' can compute separated analysis on subset of genes grouped accordingly to the value of this variable.
#' @param grouping_vars A vector of characters for subselecting genes belonging only to the specified values of the
#' varible \code{group_by}
#' @param ranges An optional \code{data.table} with values for
#'   specified intervals and associated colors.
#' @return
#'
#' @import data.table
#'
#' @export
#'
color_cell <- function(timepoint_list,
                       plot_data,
                       gene_loc_table,
                       coloring_mode='enrichment',
                       col_name=NULL,
                       colors=NULL,
                       group_by=NULL,
                       grouping_vars=NULL,
                       ranges=NULL){

  if (coloring_mode == 'mean' || coloring_mode == 'median'){
    if (is.null(col_name)){
      cat("ERROR: col_name parameter is missing")
    } else {

      if (!all(unlist(lapply(timepoint_list, function(x) col_name %in% colnames(x))))){
        cat("check that every list contains a column named as specified in col_name param")
      } else {
        if (!is.null(group_by)){
        # we are grouping by class of DEGs
          if (!all(unlist(lapply(timepoint_list, function(x) group_by %in% colnames(x))))) {
            cat("check that every list contains a column named as specified in group_by param")
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
              colors <- colors[names(colors) %in% grouping_vars]
            }
          }
        } else {
          # not grouping by any categorical variable
          colors = c("white", "darkblue")

        }

        tp_out <- list()
        for (tp in names(timepoint_list)){

          cat(paste0("Creating cell pictogram for stage: ", tp, "\n"))

          if (!is.null(group_by)){
          # create a separate plot for each category
            if (is.null(grouping_vars)) {
              # all the categories are plotted
              grouping_vars <- as.character(unique(unlist(lapply(timepoint_list, function(x) unique(x[, get(group_by)])))))
            }

            grouped_out <- list()
            for (v in grouping_vars){
              genes <- timepoint_list[[tp]][get(group_by) == v]

              fixed_ranges_f <-  discrete_symmetric_ranges(timepoint_list,
                                                           plot_data,
                                                           gene_loc_table,
                                                           col_name,
                                                           grouping_vars,
                                                           colors,
                                                           coloring_mode)

              grouped_out[[v]] <- assign_color_by_value(genes = genes,
                                                        plot_data,
                                                        gene_loc_table,
                                                        col_name,
                                                        categorical_classes = fixed_ranges_f[[v]],
                                                        coloring_mode)


            }

            tp_out[[tp]] <- grouped_out

          } else {
            # create a plot regardless the classification
            if (!is.null(grouping_vars)) {
              # only specified categories are plotted, and genes are averaged regardless any classification
              for (v in grouping_vars){
                genes <- timepoint_list[[tp]][get(group_by) == v]
                tmp <- funion(tmp, genes)
              }
            } else {
              # genes belonging to all categories are plotted
              genes <- timepoint_list[[tp]]
            }

            fixed_ranges_f_together <- discrete_symmetric_ranges(timepoint_list,
                                                                 plot_data,
                                                                 gene_loc_table,
                                                                 col_name,
                                                                 grouping_vars,
                                                                 colors,
                                                                 coloring_mode,
                                                                 together=TRUE)


            tp_out[[tp]] <- assign_color_by_value(genes = genes,
                                                  plot_data,
                                                  gene_loc_table,
                                                  categorical_classes = fixed_ranges_f_together$together,
                                                  coloring_mode)
          }
        }
        return(tp_out)
      }
    }
  } else {
    if (coloring_mode == "enrichment"){

    } else {
      cat("Wrong coloring method")
    }

  }

}
