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
  bs=40
  static_plot <- ggplot(data = gb, aes(x=x, y=y, fill=as.factor(x))) +
    scale_fill_manual(values = rev(color_vector), name="FDR", labels=rev(lab_vector)) +
    geom_bar(stat = "identity")
  #guides(color = FALSE) +
  theme(legend.title = element_text(size=bs*0.9),
        legend.text = element_text(size=bs*0.9))
  #theme(legend.text = element_text(margin = margin(r = 15, l= 4, unit = "pt"))) # + theme_void()

  # extract legend
  tmp <- ggplot_gtable(ggplot_build(static_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

start <- Sys.time()


#' Main function that creates animated cellular pictogram with the desired
#' coloring method
#'
#' @description This function creates static animated pictogram with the desired
#'   coloring method for assigning colors to subcellular localizations
#' @param data List of data structures from \link[expressyouRcell]{color_cell}
#' @param seconds A numerical value specifying the duration of each transition.
#' @param fps  A numerical value specifying the number of frame in a second.
#' @param dir A character specifying the directory path where the datasets for
#'   the animation are stored.
#' @param height gif height in pixel.
#' @param width gif width in pixel.

#' @return
#'
#' @import data.table
#' @import gifski
#' @import ggpubr
#' @import ggplot2
#'
#' @export
#'
animate <- function(data, timepoints, seconds, fps, dir, names, height = 250, width = 700, filename="animation.gif", format = "gif"){

  if (!dir.exists(dir)){
    dir.create(dir, recursive = TRUE)
  }

  frame_path <- file.path(dir, "frames")
  finaldt_path <- file.path(dir, "final_dt")

  do.call(file.remove, list(list.files(frame_path, full.names = TRUE)))

  tot_frame <- seconds*fps

  stages <- data[["final_dt"]][timepoints]

  if (!all(timepoints %in% names(data[["final_dt"]]))){
    stop("Input data do not contain these datasets. Check the timepoints parameter.")
  }

  transition=1

  ranges <- data$ranges

  l <- create_legend(color_vector = ranges$colors, lab_vector = ranges$lab)
  pb = txtProgressBar(min = 0, max = tot_frame*(length(stages)-1), style = 3)

  # for each transition
  k=0
  for (i in rev(rev(seq_len(length(stages)))[-1])){
    #cat(paste0("Creating frames for transition", i))
    t1 = stages[[i]]
    t2 = stages[[i+1]]

    t1[color_grad == "grey90", color_grad := "#E5E5E5"]
    t2[color_grad == "grey90", color_grad := "#E5E5E5"]

    together <- cbind(t1, color2 = t2[, color_grad])[, combine := toupper(paste0(color_grad, color2))]

    color_shades <- together[, .N, by = c("color_grad", "color2")
                             ][, c("color_grad", "color2")]

    dt=data.table()

    for (r in seq_len(nrow(color_shades))){
      colll = colorRampPalette(color_shades[r])(tot_frame)
      dt <- rbind(dt, t(colll))
    }

    dt[, combine := paste0(V1, get(paste0("V", tot_frame)))]

    together <- merge.data.table(together, dt, all.x = TRUE, by = "combine")

    # create and save a plot for each frame
    for (j in seq_len(tot_frame)){
      setTxtProgressBar(pb,j+k)

      try <- setDT(unique(together[, c(paste0("V", j), "value"), with=FALSE]))

      temp <- together[, c("subcell_struct", "x","y","pol","color","comb",paste0("V", j),"value"), with=FALSE]
      setnames(temp, old=paste0("V", j), new="color_grad")

      bs=25

      xmin <- min(temp$x)
      xmax <- max(temp$x)
      ymin <- min(temp$y)

      tot_l <- xmax - xmin
      trans_l <- tot_l/(length(stages)-1)

      s <- seq_len(length(stages))
      xmin_trans <- xmin + ((s-1)*trans_l)

      labels <- data.table(name=names, xmin_trans, y = ymin - 50)

      frame_l <- trans_l/tot_frame

      xminf <- xmin
      xmaxf <- xminf+(trans_l*(i-1))+(frame_l*(j))

      ecmx <- max(temp[subcell_struct == "extracellular_region"]$x)-(max(temp[subcell_struct == "extracellular_region"]$x)-min(temp[subcell_struct == "extracellular_region"]$x))/3

      ecmy <- (min(temp[subcell_struct == "extracellular_region"]$y)+max(temp[subcell_struct == "extracellular_region"]$y))/3

      plot <- ggplot(data=temp,
                  aes(x, y)) +
        #scale_fill_manual(values = colors.spe, name="FDR", labels=lab.spe) +
        scale_fill_manual(values = temp$color_grad) +
        scale_color_manual(values = rep("black", length(unique(temp$subcell_struct)))) +
        scale_size_manual(values = rep(0.005, length(temp[, first(color_grad), by=subcell_struct]$V1))) +
        geom_polygon(aes(subgroup=comb, color=temp$color_grad), fill=temp$color_grad)  +
        scale_y_reverse() +
        guides(color = FALSE, fill=FALSE) +
        theme_void()  +
        geom_text(data=labels, aes(x=xmin_trans, y=y, label=name), size=0.3*bs) +
        annotate("text", x=ecmx, y=ecmy, label="ECM", size=bs*0.2) +
        geom_segment(aes(x=xminf, y=ymin-100, xend=xmaxf, yend=ymin-100), size = 0.2*bs, lineend = "butt", color="darkgrey") +
        geom_point(data = labels, aes(xmin_trans, y-50), size = 0.2*bs, shape=19, color="grey20") +
        theme(legend.title = element_text(size=bs*0.9),
              legend.text = element_text(size=bs*0.9),
              plot.title = element_text(size=bs*0.9))

      pl <- ggpubr::ggarrange(plot + guides(fill = FALSE),
                      l,
                      widths = c(3, 1))


      if (!dir.exists(frame_path)){
        dir.create(frame_path, recursive = TRUE)
      }

      if (transition<10){
        tr_n <- paste0(0, transition)
      } else {
        tr_n <- transition
      }

      if (j<10){
        j_n <- paste0(0, j)
      } else {
        j_n <- j
      }

      ggsave(pl, filename = file.path(frame_path, paste0(tr_n, "_", j_n, ".png")), width = 10, height = 4)
    }
    k=k+50

    transition=transition+1
  }

  cat("\n")
  lf <- list.files(frame_path, full.names = T)
  png_files <- stringr::str_sort(lf, numeric = TRUE)

  if (format == "gif"){
  gifski::gifski(png_files,
                 gif_file = file.path(dir, filename),
                 width = width,
                 height = height,
                 delay = seconds/tot_frame)
  } else {
    if (format == "video"){
      av::av_encode_video(input = png_files, output = file.path(dir, filename))
    } else {
      stop("Output format unrecognized")
    }
  }


  #do.call(file.remove, list(list.files(frame_path, full.names = TRUE)))
  #unlink(frame_path, recursive = TRUE)
  end <- Sys.time()
  cat(paste0("Total time: ", eval(end-start), "\n"))
  cat(paste0("saved in ", file.path(dir, filename)))
}



# seconds <- 2
# fps <- 25

# dt_toplot = t2
#
# colors_vec <- categorical_classes$colors
# colors.spe <- colors_vec[sort(as.numeric(unique(as.vector(dt_toplot$value))), decreasing = TRUE)]
#
# lab <- as.expression(sapply(categorical_classes$lab, function(x) x))
# lab.spe <- lab[sort(as.numeric(unique(as.vector(dt_toplot$value))), decreasing = TRUE)]
#
# bs= 25
# p <- ggplot(dt_toplot, aes(x, y, color=color_grad, fill=value)) +
#   scale_fill_manual(values = colors.spe, name="FDR", labels=lab.spe) +
#   scale_color_manual(values = rep("black", length(unique(dt_toplot$subcell_struct)))) +
#   scale_size_manual(values = rep(0.005, length(dt_toplot[, first(color_grad), by=subcell_struct]$V1))) +
#   geom_polygon(aes(subgroup=comb)) +
#   scale_y_reverse() +
#   guides(color = FALSE) +
#   theme_void()  +
#   theme(legend.title = element_text(size=bs*0.9),
#         legend.text = element_text(size=bs*0.9))
# p
