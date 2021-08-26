#library(apng)
library(gifski)
library(data.table)
library(ggplot2)
library(ggpubr)

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

start <- Sys.time()

seconds <- 2
fps <- 25
tot_frame <- seconds*fps

setwd(dir = "C:\\Users/martina/Desktop/trial.png/")
dir <- "C:\\Users/martina/Desktop/trial.png/package_ok/wound_enr/"
files <- list.files(dir, pattern = "-.RData")
stages <- stringr::str_sort(files, numeric = TRUE)

transition=1

for (i in rev(rev(seq_len(length(stages)))[-1])){
  t1 = stages[[i]]
  t2 = stages[[i+1]]
  
 
  load(file.path(dir, t1))
  t1 = tmp
  load(file.path(dir, t2))
  t2 = tmp

  t1[enr_color == "grey90", enr_color := "#E5E5E5"]
  t2[enr_color == "grey90", enr_color := "#E5E5E5"]
  
  together <- cbind(t1, color2 = t2[, enr_color])[, combine := toupper(paste0(enr_color, color2))]
  
  color_shades <- together[, .N, by = c("enr_color", "color2")
                           ][, c("enr_color", "color2")]
  
  dt=data.table()
  
  for (i in seq_len(nrow(color_shades))){
    colll = colorRampPalette(color_shades[i])(tot_frame)
    dt <- rbind(dt, t(colll))
  }
  
  dt[, combine :=paste0(V1, get(paste0("V", tot_frame)))]
  
  together <- merge.data.table(together, dt, all.x = TRUE, by = "combine")
  
  for (i in seq_len(tot_frame)){
    try <- setDT(unique(together[, c(paste0("V", i), "value"), with=FALSE]))
    
    temp <- together[, c("subcell_struct", "x","y","pol","color","comb",paste0("V", i),"value"), with=FALSE]
    setnames(temp, old=paste0("V", i), new="enr_color")
    
    load(file = "categorical_classes.RData")
    
    l <- create_legend(color_vector = categorical_classes$colors, lab_vector = categorical_classes$lab)
    
    categorical_classes[colors == "grey90", colors := "#E5E5E5"]
    categorical_classes$values <- as.factor(categorical_classes$values)
    
    categorical_classes2 <- merge.data.table(categorical_classes, try, by.x = "values", by.y = "value", all = TRUE)
    
    categorical_classes2[is.na(get(paste0("V", i))), paste0("V", i) := colors]
    categorical_classes2$colors <- NULL
    setnames(categorical_classes2, old = paste0("V", i), new = "colors")
    categorical_classes2$colors <- toupper(categorical_classes2$colors)
    categorical_classes <- categorical_classes2
    
    colors_vec <- toupper(categorical_classes$colors)
    colors.spe <- colors_vec[sort(as.numeric(unique(as.vector(temp$value))), decreasing = TRUE)]
    
    lab <- as.expression(sapply(categorical_classes$lab, function(x) x))
    lab.spe <- lab[sort(as.numeric(unique(as.vector(temp$value))), decreasing = TRUE)]
    
    
    
    bs=25
    p <- ggplot(temp, 
                aes(x, y, color=enr_color)) +
      #scale_fill_manual(values = colors.spe, name="FDR", labels=lab.spe) +
      scale_fill_manual(values = temp$enr_color) +
      scale_color_manual(values = rep("black", length(unique(temp$subcell_struct)))) +
      scale_size_manual(values = rep(0.005, length(temp[, first(enr_color), by=subcell_struct]$V1))) +
      geom_polygon(aes(subgroup=comb), fill=temp$enr_color) +
      scale_y_reverse() +
      guides(color = FALSE, fill=FALSE) +
      theme_void()  +
      theme(legend.title = element_text(size=bs*0.9),
            legend.text = element_text(size=bs*0.9))
    
    pl <- ggarrange(p +
              guides(fill = FALSE),
              #ggtitle(stage),
              l,
              widths = c(3, 1))
    
    ggsave(pl, filename = file.path(getwd(), "frames", paste0(transition, "_", i, ".png")), width = 10, height = 4)
  }
  
  
  transition=transition+1
}

lf <- list.files(file.path(getwd(), "frames"), full.names = T)
png_files <- stringr::str_sort(lf, numeric = TRUE)
gifski(png_files, gif_file = "animation2.gif", width = 800, height = 350, delay = 2/tot_frame)

end <- Sys.time()

end-start







# dt_toplot = t2
# 
# colors_vec <- categorical_classes$colors
# colors.spe <- colors_vec[sort(as.numeric(unique(as.vector(dt_toplot$value))), decreasing = TRUE)]
# 
# lab <- as.expression(sapply(categorical_classes$lab, function(x) x))
# lab.spe <- lab[sort(as.numeric(unique(as.vector(dt_toplot$value))), decreasing = TRUE)]
# 
# bs= 25
# p <- ggplot(dt_toplot, aes(x, y, color=enr_color, fill=value)) +
#   scale_fill_manual(values = colors.spe, name="FDR", labels=lab.spe) +
#   scale_color_manual(values = rep("black", length(unique(dt_toplot$subcell_struct)))) +
#   scale_size_manual(values = rep(0.005, length(dt_toplot[, first(enr_color), by=subcell_struct]$V1))) +
#   geom_polygon(aes(subgroup=comb)) +
#   scale_y_reverse() +
#   guides(color = FALSE) +
#   theme_void()  +
#   theme(legend.title = element_text(size=bs*0.9),
#         legend.text = element_text(size=bs*0.9))
# p
