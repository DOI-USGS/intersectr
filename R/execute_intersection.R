#' Execute Intersection
#'
#' @description Traverses many time steps of data applying
#' intersection weights for each time step.
#'
#' @param nc_file character path or OPeNDAP URL to NetCDF source
#' @param variable_name character name of variable in NetCDF source to access
#' @param intersection_weights data.frame as output by
#' \code{\link{calculate_area_intersection_weights}}
#' @param cell_geometry sf data.frame with x_ind and y_ind collumns as output by
#' \code{\link{create_cell_geometry}}
#' @param x_var character variable with X axis data
#' @param y_var character variable with Y axis data
#' @param t_var character variable with T axis data
#' @param start_date not implemented
#' @param end_date not implemented
#' @importFrom dplyr right_join mutate summarise group_by ungroup
#' @importFrom RNetCDF open.nc dim.inq.nc var.inq.nc att.get.nc var.get.nc utcal.nc
#' @export
#' @examples
#'
#' library(dplyr)
#' library(sf)
#'
#' variable_name <- "precipitation_amount"
#' x_var <- "lon"
#' y_var <- "lat"
#' t_var <- "day"
#'
#' nc_file <- system.file("extdata/metdata.nc", package = "intersecter")
#' nc <- RNetCDF::open.nc(nc_file)
#'
#' x <- RNetCDF::var.get.nc(nc, x_var)
#' y <- RNetCDF::var.get.nc(nc, y_var)
#'
#' in_prj <- "+init=epsg:4326"
#' out_prj <- "+init=epsg:5070"
#'
#' geom <- read_sf(system.file("shape/nc.shp", package = "sf")) %>%
#'   st_transform(out_prj) %>%
#'   st_buffer(0)
#'
#' cell_geometry <- create_cell_geometry(x, y, in_prj, geom)
#'
#' data_source_cells <- st_sf(select(cell_geometry, grid_ids))
#' target_polygons <- st_sf(select(geom, CNTY_ID))
#'
#' st_agr(data_source_cells) <- "constant"
#' st_agr(target_polygons) <- "constant"
#'
#' area_weights <- calculate_area_intersection_weights(
#'   data_source_cells,
#'   target_polygons)
#'
#' intersected <- execute_intersection(nc_file, variable_name, area_weights,
#'                                     cell_geometry, x_var, y_var, t_var)
#'
#' x_inds <- seq(min(cell_geometry$x_ind), max(cell_geometry$x_ind), 1)
#' y_inds <- seq(min(cell_geometry$y_ind), max(cell_geometry$y_ind), 1)
#'
#' ids <- intersecter:::get_ids(length(x_inds), length(y_inds))
#'
#' grid_data <- RNetCDF::var.get.nc(nc, variable_name,
#'                                  start = c(min(x_inds), min(y_inds), 5),
#'                                  count = c(length(x_inds), length(y_inds), 1))
#'
#' grid_data <- data.frame(grid_data = matrix(t(grid_data), # note t() here!!!
#'                                            ncol = 1,
#'                                            byrow = TRUE),
#'                         grid_ids = matrix(ids, ncol = 1))
#'
#' grid_data$grid_data[grid_data$grid_data < 0] <- NA
#'
#' grid_data <- left_join(cell_geometry, grid_data, by = "grid_ids")
#'
#' geom_data <- select(geom, CNTY_ID) %>%
#'   left_join(data.frame(CNTY_ID = as.numeric(names(intersected)[2:ncol(intersected)]),
#'                        poly_data = as.numeric(intersected[5, 2:ncol(intersected)])), by = "CNTY_ID")
#'
#' geom_data <- st_transform(geom_data, st_crs(grid_data))
#'
#' breaks <- c(0, 0.5, 1, 2, 3, 5, 10)
#'
#' plot(grid_data$geometry)
#'
#' plot(grid_data["grid_data"], border = NA, breaks = breaks, add = TRUE)
#'
#' plot(geom_data["poly_data"], breaks = breaks, add = TRUE)
#'
#' plot(grid_data$geometry, lwd = 0.1, add = TRUE)

execute_intersection <- function(nc_file, variable_name,
                                 intersection_weights,
                                 cell_geometry,
                                 x_var, y_var, t_var,
                                 start_date = NULL, end_date = NULL) {

  names(intersection_weights) <- c("grid_ids", "poly_id", "w")

  nc <- open.nc(nc_file)
  nc_var_info <- var.inq.nc(nc, variable_name)

  if(nc_var_info$ndims != 3) stop("only 3d variables are supported")

  x_var_info <- var.inq.nc(nc, x_var)
  y_var_info <- var.inq.nc(nc, y_var)
  t_var_info <- var.inq.nc(nc, t_var)

  if(x_var_info$ndims > 1 | y_var_info$ndims > 1 | t_var_info$ndims > 1)
    stop("only 1d coordinate variables supported")

  t_dim_info <- dim.inq.nc(nc, t_var_info$dimids)

  dimid_order <- match(nc_var_info$dimids,
                       c(x_var_info$dimids,
                         y_var_info$dimids,
                         t_var_info$dimids))

  x_inds <- seq(min(cell_geometry$x_ind), max(cell_geometry$x_ind), 1)
  y_inds <- seq(min(cell_geometry$y_ind), max(cell_geometry$y_ind), 1)

  ids <- get_ids(length(x_inds), length(y_inds))

  # t_inds <- which(t_vals >= start_date & t_vals <= end_date)

  if(is.null(start_date) & is.null(end_date)) {
    out_nrows <- t_dim_info$length
    out_ncols <- length(unique(intersection_weights[,2][[1]]))

    out <- matrix(nrow = out_nrows, ncol = out_ncols)

    for(i in 1:out_nrows) {
      try_backoff({
        i_data <- var.get.nc(nc, variable_name,
                             start = c(min(x_inds), min(y_inds), i)[dimid_order],
                             count = c(length(x_inds), length(y_inds), 1)[dimid_order])

        i_data <- data.frame(d = matrix(t(i_data), # note t() here!!!
                                        ncol = 1,
                                        byrow = TRUE),
                             grid_ids = matrix(ids, ncol = 1))

        i_data <- right_join(i_data,
                             intersection_weights,
                             by = c("grid_ids")) %>%
          mutate(d = d * w) %>%
          group_by(poly_id) %>%
          summarise(d = sum(d, na.rm = T)/sum(w, na.rm = T)) %>%
          ungroup()

        out[i, ] <- i_data$d
      })
    }
  } else {
    stop("date filtering not supported yet")
  }
  time_units <- att.get.nc(nc, t_var_info$name, "units")
  time_vals <- var.get.nc(nc, t_var_info$name)
  time_steps <- utcal.nc(time_units, time_vals, "c")
  time_steps <- data.frame(time_stamp = time_steps)

  out <- data.frame(out)

  out_names <- as.character(unique(i_data$poly_id))

  out <- select(out, which(!is.na(out_names)))

  out_names <- out_names[which(!is.na(out_names))]

  names(out) <- out_names

  return(cbind(time_steps, out))
}
