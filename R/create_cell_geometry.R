#' Create Cell Geometry
#'
#' @description Creates cell geometry from vectors of x and y positions.
#'
#' @param x numeric center positions of x indices
#' @param y numeric center positions of y indices
#' @param prj character proj4 string for x and y
#' @param geom sf data.frame with geometry that cell geometry should cover
#' @param buffer_dist a distance to buffer the cell geometry in units of geom projection
#'
#' @details Intersection is performed with cell centers then geometry is constructed.
#' A buffer may be required to fully cover geometry with cells.
#'
#' @importFrom sf st_sf st_as_sf st_as_sfc st_bbox st_transform st_buffer st_geometry st_union
#' st_voronoi st_cast st_intersection st_join st_contains st_convex_hull st_distance
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
create_cell_geometry <- function(x, y, prj, geom = NULL, buffer_dist = 0) {

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
  dif_dist_x <- diff(x)
  dif_dist_y <- diff(y)

  if (any(diff(dif_dist_x) > 1e-10) | any(diff(dif_dist_y) > 1e-10)) {
    stop("only regular rasters supported")
  } else {
    dif_dist_x <- mean(dif_dist_x)
    dif_dist_y <- mean(dif_dist_y)
  }

  sf_points <- construct_points(x, y, prj)

  if (!is.null(geom)) {
    # intersect in projection of geometry
    sf_points_filter <- st_intersection(
      st_transform(sf_points, st_crs(geom)),
      st_buffer(st_as_sfc(st_bbox(geom)), buffer_dist))

    # grab all the rows and cols needed.
    sf_points <- filter(sf_points,
                          x_ind >= min(sf_points_filter$x_ind) &
                          x_ind <= max(sf_points_filter$x_ind) &
                          y_ind >= min(sf_points_filter$y_ind) &
                          y_ind <= max(sf_points_filter$y_ind))
  }

  x <- x[min(sf_points$x_ind):max(sf_points$x_ind)]
  y <- y[min(sf_points$y_ind):max(sf_points$y_ind)]

  sf_polygons <- get_ids(length(x), length(y))
  dim(sf_polygons) <- c(x = length(y), y = length(x))

  sf_polygons <- st_as_stars(list(sf_polygons = sf_polygons),
                         dimensions = st_dimensions(x = as.numeric(y),
                                                    y = as.numeric(x),
                                                    .raster = c("x", "y")))

  sf_polygons <- st_as_sf(sf_polygons, as_points = FALSE)

  sf::st_crs(sf_polygons) <- st_crs(prj)

  # sf_polygons <- st_as_sf(
  #   stars::st_as_stars(raster::raster(sf_polygons,
  #                                     xmn = min(x) - 0.5 * dif_dist_x,
  #                                     xmx = max(x) + 0.5 * dif_dist_x,
  #                                     ymn = min(y) + 0.5 * dif_dist_y,
  #                                     ymx = max(y) - 0.5 * dif_dist_y,
  #                                     crs = prj)),
  #   as_points = FALSE)

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
  x_ind <- matrix(rep(c(1:ncol(x_vals)), nrow(x_vals)),
                  nrow = nrow(x_vals), ncol = ncol(x_vals),
                  byrow = TRUE)
  y_ind <- matrix(rep(c(1:nrow(x_vals)), ncol(x_vals)),
                  nrow = nrow(x_vals), ncol = ncol(x_vals),
                  byrow = FALSE)

  sf_points <- st_as_sf(data.frame(x = matrix(x_vals, ncol = 1),
                                   y = matrix(y_vals, ncol = 1),
                                   x_ind = matrix(x_ind, ncol = 1),
                                   y_ind = matrix(y_ind, ncol = 1)),
                        coords = c("x", "y"),
                        crs = prj,
                        agr = "constant")
}

get_ids <- function(x_size, y_size) {
  matrix(as.numeric(seq(1, x_size * y_size)),
         nrow = y_size, ncol = x_size,
         byrow = T)
}
