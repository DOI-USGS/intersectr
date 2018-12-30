#' Create Cell Geometry
#'
#' @description Creates cell geometry from vectors of x and y positions.
#'
#' @param col_coords numeric center positions of x indices
#' @param row_coords numeric center positions of y indices
#' @param prj character proj4 string for x and y
#' @param geom sf data.frame with geometry that cell geometry should cover
#' @param buffer_dist a distance to buffer the cell geometry in units of geom projection
#'
#' @details Intersection is performed with cell centers then geometry is constructed.
#' A buffer may be required to fully cover geometry with cells.
#'
#' @importFrom sf st_sf st_as_sf st_as_sfc st_bbox st_transform st_buffer st_intersection
#' st_join st_contains "st_crs<-"
#' @importFrom dplyr filter
#' @importFrom stars st_as_stars st_dimensions
#' @export
#' @examples
#' library(RNetCDF)
#' library(sf)
#' variable_name <- "precipitation_amount"
#' nc <- open.nc(system.file("extdata/metdata.nc", package = "intersectr"))
#' nc_file <- file.inq.nc(nc)
#' nc_dim <- lapply(seq(0, nc_file$ndims-1), function(d) dim.inq.nc(nc, d))
#' nc_var <- lapply(seq(0, nc_file$nvars-1), function(v) var.inq.nc(nc, v))
#' var_names <- sapply(nc_var, function(v) v$name)
#' print(var_names)
#' var <- which(var_names == variable_name)
#' var_atts <- lapply(seq(0, nc_var[[var]]$natts - 1), function(a) att.inq.nc(nc, var - 1, a))
#' var_att_names <- sapply(var_atts, function(a) a$name)
#' coordinates <- att.get.nc(nc, var - 1, which(var_att_names == "coordinates") - 1)
#' print(coordinates)
#' # Could iomplement general handling for this... hand set for now.
#' # Either use coordinate variables or standard names for auxiliary coordinate variables.
#' x_var <- "lon"
#' y_var <- "lat"
#'
#' grid_mapping <- att.get.nc(nc, var - 1, "grid_mapping")
#'
#' if (!is.null(grid_mapping)) {
#'   prj <- ncdfgeom::get_prj(gm_atts)
#' } else {
#'   prj <- "+init=epsg:4326"
#'   warning("No grid mapping found. Assuming WGS84")
#' }
#'
#' x <- var.get.nc(nc, x_var)
#' y <- var.get.nc(nc, y_var)
#'
#' geom <- read_sf(system.file("shape/nc.shp", package = "sf"))
#'
#' cell_geometry <- create_cell_geometry(x, y, prj, geom, 0)
#'
#' plot(cell_geometry$geometry, lwd = 0.25)
#' plot(st_transform(geom$geometry, prj), add = TRUE)
#'
create_cell_geometry <- function(col_coords, row_coords, prj, geom = NULL, buffer_dist = 0) {

  # For 2d lat/lon
  # x <- matrix(rep(c(1:ncol(lon)), nrow(lon)),
  #             nrow = nrow(lon), ncol = ncol(lon),
  #             byrow = TRUE)
  #
  # y <- matrix(rep(c(1:nrow(lon)), ncol(lon)),
  #             nrow = nrow(lon), ncol = ncol(lon),
  #             byrow = FALSE)
  #
  # sf_points <- data.frame(x = matrix(x, ncol = 1),
  #                         y = matrix(y, ncol = 1),
  #                         lon = matrix(lon, ncol = 1),
  #                         lat = matrix(lat, ncol = 1))

  # cell size
  dif_dist_x <- diff(col_coords)
  dif_dist_y <- diff(row_coords)

  if (any(diff(dif_dist_x) > 1e-10) | any(diff(dif_dist_y) > 1e-10)) {
    stop("only regular rasters supported")
  } else {
    dif_dist_x <- mean(dif_dist_x)
    dif_dist_y <- mean(dif_dist_y)
  }

  sf_points <- construct_points(col_coords, row_coords, prj)

  if (!is.null(geom)) {
    # intersect in projection of geometry
    sf_points_filter <- st_intersection(
      st_transform(sf_points, st_crs(geom)),
      st_buffer(st_as_sfc(st_bbox(geom)), buffer_dist))

    # grab all the rows and cols needed.
    sf_points <- filter(sf_points,
                          col_ind >= min(sf_points_filter$col_ind) &
                          col_ind <= max(sf_points_filter$col_ind) &
                          row_ind >= min(sf_points_filter$row_ind) &
                          row_ind <= max(sf_points_filter$row_ind))
  }

  col_coords <- col_coords[min(sf_points$col_ind):max(sf_points$col_ind)]
  row_coords <- row_coords[min(sf_points$row_ind):max(sf_points$row_ind)]

  sf_polygons <- get_ids(length(col_coords), length(row_coords))
  dim(sf_polygons) <- c(x = length(col_coords), y = length(row_coords))

  sf_polygons <- st_as_stars(list(sf_polygons = sf_polygons),
                         dimensions = st_dimensions(x = as.numeric(col_coords),
                                                    y = as.numeric(row_coords),
                                                    .raster = c("x", "y"),
                                                    cell_midpoints = TRUE))

  st_crs(sf_polygons) <- st_crs(prj)

  sf_polygons <- st_as_sf(sf_polygons, as_points = FALSE)

  names(sf_polygons)[1] <- "grid_ids"

  sf_polygons <- st_join(sf_polygons, sf_points, join = st_contains)

  sf::st_agr(sf_polygons) <- "constant"

  return(sf_polygons)
}

construct_points <- function(x, y, prj) {
  x_vals <- matrix(x, nrow = length(y), ncol = length(x),
                   byrow = T)
  y_vals <- matrix(y, nrow = length(y), ncol = length(x),
                   byrow = F)
  col_ind <- matrix(rep(c(1:ncol(x_vals)), nrow(x_vals)),
                  nrow = nrow(x_vals), ncol = ncol(x_vals),
                  byrow = TRUE)
  row_ind <- matrix(rep(c(1:nrow(x_vals)), ncol(x_vals)),
                  nrow = nrow(x_vals), ncol = ncol(x_vals),
                  byrow = FALSE)

  sf_points <- st_as_sf(data.frame(x = matrix(x_vals, ncol = 1),
                                   y = matrix(y_vals, ncol = 1),
                                   col_ind = matrix(col_ind, ncol = 1),
                                   row_ind = matrix(row_ind, ncol = 1)),
                        coords = c("x", "y"),
                        crs = prj,
                        agr = "constant")
}

get_ids <- function(x_size, y_size) {
  matrix(as.numeric(seq(1, x_size * y_size)),
         nrow = x_size, ncol = y_size,
         byrow = T)
}
