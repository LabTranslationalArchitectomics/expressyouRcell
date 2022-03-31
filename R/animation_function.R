#' Create a legend for the animated pictures
#'
#' @description This function creates a legend for the animated gif files.
#' @param color_vector A vector of hex color codes.
#' @param lab_vector A vector of labels for the legend.
#' @param title A vector for the legend title.
#' @return a gtable with the legend to be added on the animated pictures
#'
#' @import ggplot2
#'
create_legend <- function(color_vector, lab_vector, title){

  gb <- data.table(x=seq(1:length(color_vector)),
                   y=seq(1:length(color_vector)))

  # create legend by combing stages
  bs=30
  static_plot <- ggplot(data = gb, aes(x=x, y=y, fill=as.factor(x))) +
    scale_fill_manual(values = rev(color_vector), name=title, labels=rev(lab_vector)) +
    geom_bar(stat = "identity") +
  #guides(color = FALSE) +
  theme(legend.title = element_text(size=bs*0.5),
        legend.text = element_text(size=bs*0.5))
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
#' @description This function creates static animated pictogram with the desired coloring method for assigning colors to
#'   subcellular localizations
#' @param data List of data structures from \link[expressyouRcell]{color_cell}
#' @param timepoints A character vector with the names of the timepoints data structure that will be part of the final
#'   animation.
#' @param seconds A numerical value specifying the duration of each transition.
#' @param fps  A numerical value specifying the number of frame in a second.
#' @param input_dir A character specifying the directory path where the datasets for the animation are stored.
#' @param names A character vector specifying the names of the labels for creating the timeline on the top of the output
#'   gif or movie. The vector should match the timepoints label specified in the timepoints parameter above.
#' @param height output height in pixel.
#' @param width output width in pixel.
#' @param filename A character specifying the name of the output file, default is "animation"
#' @param format either "gif" or "video", which respecitivley save the animation as an animated gif picture or as a
#'   short movie in mp4 format. Default is "gif".
#'
#' @import data.table
#' @import gifski
#' @import av
#' @import ggpubr
#' @import ggplot2
#' @import gridExtra
#'
#' @export
#'
animate <- function(data, timepoints, seconds, fps, input_dir, names, height = 25, width = 30, filename="animation", format){

  if (!dir.exists(input_dir)){
    dir.create(input_dir, recursive = TRUE)
  }

  frame_path <- file.path(input_dir, "frames")
  finaldt_path <- file.path(input_dir, "final_dt")

  do.call(file.remove, list(list.files(frame_path, full.names = TRUE)))

  tot_frame <- seconds*fps

  stages <- data[["final_dt"]][timepoints]

  if (!all(timepoints %in% names(data[["final_dt"]]))){
    stop("Input data do not contain these datasets. Check the timepoints parameter.")
  }

  transition=1

  ranges <- data$ranges

  l <- create_legend(color_vector = ranges$colors,
                     lab_vector = ranges$lab,
                     title=data$legend_label)
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

      bs=30

      xmin <- min(temp$x)+20
      xmax <- max(temp$x)-20
      ymin <- min(temp$y)

      tot_l <- xmax - xmin
      trans_l <- tot_l/(length(stages)-1)

      s <- seq_len(length(stages))
      xmin_trans <- xmin + ((s-1)*trans_l)

      labels <- data.table(name=names, xmin_trans, y = ymin - 50)

      frame_l <- trans_l/tot_frame

      xminf <- xmin
      xmaxf <- xminf+(trans_l*(i-1))+(frame_l*(j))

      ecmx <- min(temp[subcell_struct == "extracellular_region"]$x) + (max(temp[subcell_struct == "extracellular_region"]$x)-
                                                                              min(temp[subcell_struct == "extracellular_region"]$x))/1.5

      ecmy <- (min(temp[subcell_struct == "extracellular_region"]$y)+
                 max(temp[subcell_struct == "extracellular_region"]$y))/3

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
        geom_text(data=labels, aes(x=xmin_trans, y=y, label=name), size=0.2*bs) +
        annotate("text", x=ecmx, y=ecmy, label="ECM", size=bs*0.2, angle=90) +
        geom_segment(aes(x=xminf, y=ymin-100, xend=xmaxf, yend=ymin-100), size = 0.2*bs, lineend = "butt", color="darkgrey") +
        geom_point(data = labels, aes(xmin_trans, y-50), size = 0.2*bs, shape=19, color="grey20") +
        theme(legend.title = element_text(size=bs*0.9),
              legend.text = element_text(size=bs*0.9),
              plot.title = element_text(size=bs*0.9))

      if (all(data$cell_type == "cell")){
        w <- c(3,1)
      } else {
        if (all(data$cell_type == "neuron")) {
          w <- c(4,1)
        } else {
          if (all(data$cell_type == "fibroblast")) {
            w <- c(3,1)
          } else {
            if (all(data$cell_type == "microglia")) {
              w <- c(2,1)
            } else {
              stop("No available pictogram with this name")
            }
          }
        }
      }

      pl <- grid.arrange(plot + guides(fill = FALSE), l, ncol=2, widths=w)

      # pl <- ggpubr::ggarrange(plot + guides(fill = FALSE),
      #                 l,
      #                 widths = w)


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

      ggsave(pl, filename = file.path(frame_path, paste0(tr_n, "_", j_n, ".png")),
             width = width,
             height = height)
    }
    k=k+(fps*seconds)

    transition=transition+1
  }

  cat("\n")
  lf <- list.files(frame_path, full.names = T)
  png_files <- stringr::str_sort(lf, numeric = TRUE)

  if (format == "gif"){
  gifski::gifski(png_files,
                 gif_file = file.path(input_dir, paste0(filename, ".gif")),
                 width = width,
                 height = height,
                 delay = seconds/tot_frame)
   } else {
    if (format == "video"){
      av::av_encode_video(input = png_files, output = file.path(input_dir, paste0(filename, ".mp4")))
    } else {
      stop("Output format unrecognized")
    }
  }


  #do.call(file.remove, list(list.files(frame_path, full.names = TRUE)))
  unlink(frame_path, recursive = TRUE)
  end <- Sys.time()
  cat(paste0("Total time: ", eval(end-start), "\n"))
  cat(paste0("saved in ", file.path(input_dir, paste0(filename, "\n"))))
}
