#' Read svg files
#'
#' @description  Read coordinates from svg file and assign them predefined colors
#' @param f_svg A character string giving the path to svg file where polygon coordinates are saved
#' @param subcell_struct A character specifying the name of a subcellular structure
#' @param color A character specifying the color for a subcellular structure
#' @return A \code{data.table} with x and y coordinates, number of polygons and colors
#'
#' @import data.table
#' @import grImport2
#' @import rsvg
#'
#' @examples
#' #read_coordinates(f_svg = "file-path/actin.svg", subcell_struct="actin", color = colors_shapes[["actin"]])
read_coordinates <- function(file_svg) {
  file_svg <- rsvg::rsvg_svg(file_svg)
  rtc <- rawToChar(file_svg)
  pic <- grImport2::readPicture(rtc)
  pic@content[[1]]@content[[1]]@d@segments[[3]]

  if(exists("dt_pol")){rm(dt_pol)}
  for(j in 1:length(pic@content[[1]]@content)){
    for(i in 2:length(pic@content[[1]]@content[[j]]@d@segments)){
      dt_temp <- data.table(x = pic@content[[1]]@content[[j]]@d@segments[[i]]@x,
                            y = pic@content[[1]]@content[[j]]@d@segments[[i]]@y)[, pol := as.character(j)]
      if(!exists("dt_pol")){
        dt_pol <- dt_temp
      } else {
        dt_pol <- rbind(dt_pol, dt_temp)
      }
    }
  }

  dt_pol[, pol := factor(pol, levels = unique(pol))]
  dt_pol[order(pol)]

  return(dt_pol)
}

#' Create data.table for plotting the legend of organelles
#'
#' @description This function creates animation between stages
#' @param svg_folder A character vector with the path of the folder containing the svg files for the legend
#'
#' @examples
#' #create_legend_dt(svg_folder = "legend")
#'
#' @import data.table
#'
#' @return the \code{data.table} with the polygon coordinates to be plotted in the legend
#'
#' @export
create_legend_dt <- function(svg_folder){

  subcell_structures <- tstrsplit(list.files(file.path(getwd(), svg_folder)), '\\.')[[1]]

  coords_dt <- data.table()
  for (subcell_struct in subcell_structures){
    cat(paste0(subcell_struct, "\n"))
    dt_pol <- read_coordinates(file_svg = file.path(svg_folder, paste0(subcell_struct, '.svg')))
    dt_pol$subcell_struct <- subcell_struct
    coords_dt <- rbind(coords_dt, dt_pol)
  }

  coords_dt[, comb := factor(paste0(subcell_struct, pol), levels = unique(paste0(subcell_struct, pol)))]
  return(coords_dt)
}

#' Create neuron data.table to be plotted
#'
#' @description This function generates the main \code{data.table} containing all the coordinates of the polygons to be plotted
#' @param svg_folder A character string giving the path to folder where svg files are stored
#' @return the \code{data.table} with the polygon coordinates to be plotted
#' @examples
#' #create_neuron(svg_folder="single_shapes")
#'
#' @import data.table
#'
#'
#' @export
create_cell_dt <- function(svg_folder, colors_shapes, order_levels){

  coords_dt <- data.table()
  for (subcell_struct in names(colors_shapes)){
    cat(paste0(subcell_struct, "\n"))
    #f_svg <- rsvg::rsvg_svg(svg = file.path(svg_folder, paste0(subcell_struct, '.svg')))
    dt_pol <- read_coordinates(file_svg = file.path(svg_folder, paste0(subcell_struct, '.svg')))
    dt_pol$subcell_struct <- subcell_struct
    dt_pol$color <- colors_shapes[[subcell_struct]]
    coords_dt <- rbind(coords_dt, dt_pol)
  }

  coords_dt <- coords_dt[, subcell_struct := factor(subcell_struct,
                                               levels = order_levels)
                             ][, pol := factor(pol)
                               ][order(subcell_struct)
                                 ][, comb := factor(paste0(subcell_struct, pol), levels = unique(paste0(subcell_struct, pol)))]
  return(coords_dt)
}
