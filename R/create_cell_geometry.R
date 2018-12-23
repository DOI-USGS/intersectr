#' Create Cell Geometry
#'
#' @description Creates cell geometry from vectors of x and y positions.
#'
#' @param x numeric center positions of x indices
#' @param y numeric center positions of y indices
#' @param prj character proj4 string for x and y
#'
#' @importFrom sf st_sf st_as_sf st_as_sfc st_bbox st_transform st_buffer st_geometry st_union
#' st_voronoi st_cast st_intersection st_join st_contains
#'
#' @export
#' @examples
#' library(RNetCDF)
#' variable_name <- "precipitation_amount"
#' nc <- open.nc(system.file("extdata/metdata.nc", package = "intersecter"))
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
#' if(!is.null(grid_mapping)) {
#'   in_prj <- ncdfgeom::get_prj(gm_atts)
#' } else {
#'   in_prj <- "+init=epsg:4326"
#'   warning("No grid mapping found. Assuming WGS84")
#' }
#'
#' x <- var.get.nc(nc, x_var)
#' y <- var.get.nc(nc, y_var)
#'
#' geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))
#'
#' out_prj <- "+init=epsg:5070"
#'
#' cell_geometry <- create_cell_geometry(x, y, in_prj, out_prj, geom, 1000)
#'
#' plot(cell_geometry$geometry, lwd = 0.25)
#' plot(st_transform(geom$geometry, out_prj), add = TRUE)
#'
create_cell_geometry <- function(x, y, in_prj, out_prj, geom = NULL, buffer_dist = 0) {

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

  x_vals <- matrix(x, nrow = length(y), ncol = length(x),
                   byrow = T)
  y_vals <- matrix(y, nrow = length(y), ncol = length(x),
                   byrow = F)
  # ids <- matrix(seq(1, length(x) * length(y)),
  #               nrow = length(y), ncol = length(x),
  #               byrow = T)
  x_ind <- matrix(rep(c(1:ncol(x_vals)), nrow(x_vals)),
                  nrow = nrow(x_vals), ncol = ncol(x_vals),
                  byrow = TRUE)

  y_ind <- matrix(rep(c(1:nrow(x_vals)), ncol(x_vals)),
                  nrow = nrow(x_vals), ncol = ncol(x_vals),
                  byrow = FALSE)

  sf_points <- st_as_sf(data.frame(x = matrix(x_vals, ncol = 1),
                                   y = matrix(y_vals, ncol = 1),
                                   # grid_id = matrix(ids, ncol = 1),
                                   x_ind = matrix(x_ind, ncol = 1),
                                   y_ind = matrix(y_ind, ncol = 1)),
                        coords = c("x", "y"),
                        crs = in_prj,
                        agr = "constant") %>%
    st_transform(out_prj)

  if(!is.null(geom)) {
    crop_box <- st_transform(geom, out_prj) %>%
      st_bbox() %>%
      st_as_sfc() %>%
      st_buffer(buffer_dist)
    sf_points <- st_intersection(sf_points, crop_box)
  } else {
    crop_box <- st_as_sfc(st_bbox(sf_points)) # should maybe buffer half a cell...
  }

  sf_polygons <- st_geometry(sf_points) %>%
    st_union() %>%
    st_voronoi() %>%
    st_cast() %>%
    st_buffer(0) %>%
    st_intersection(crop_box) %>%
    st_sf()

  sf_polygons <- st_join(sf_polygons, sf_points, join = st_contains) %>%
    mutate(grid_id = 1:nrow(sf_polygons)) # see other method for grid id above.

}
