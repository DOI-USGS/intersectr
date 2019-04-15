#' Create Cell Geometry
#'
#' @description Creates cell geometry from vectors of X and Y positions.
#'
#' @param X_coords numeric center positions of X axis indices
#' @param Y_coords numeric center positions of Y axis indices
#' @param prj character proj4 string for x and y
#' @param geom sf data.frame with geometry that cell geometry should cover
#' @param buffer_dist a distance to buffer the cell geometry in units of geom projection
#'
#' @details Intersection is performed with cell centers then geometry is constructed.
#' A buffer may be required to fully cover geometry with cells.
#'
#' @importFrom sf st_sf st_as_sf st_transform st_buffer st_join st_contains
#' "st_crs<-" st_union st_bbox
#' @importFrom dplyr filter
#' @importFrom stars st_as_stars st_dimensions
#' @export
#' @examples
#' nc <- RNetCDF::open.nc(system.file("extdata/metdata.nc", package = "intersectr"))
#' ncmeta::nc_vars(nc)
#' variable_name <- "precipitation_amount"
#' cv <- ncmeta::nc_coord_var(nc, variable_name)
#'
#' x <- RNetCDF::var.get.nc(nc, cv$X, unpack = TRUE)
#' y <- RNetCDF::var.get.nc(nc, cv$Y, unpack = TRUE)
#'
#' prj <- ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc))
#'
#' geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))
#' geom <- sf::st_transform(geom, 5070)
#'
#' cell_geometry <- create_cell_geometry(x, y, prj, geom, 0)
#'
#' plot(sf::st_geometry(cell_geometry), lwd = 0.25)
#' plot(sf::st_transform(sf::st_geometry(geom), prj), add = TRUE)
#'
create_cell_geometry <- function(X_coords, Y_coords, prj, geom = NULL, buffer_dist = 0) {

  X_ind <- Y_ind <- NULL

  if (!is.array(X_coords) || length(dim(X_coords)) == 1) {
    array_mode <- FALSE

    # cell size
    dif_dist_x <- diff(X_coords)
    dif_dist_y <- diff(Y_coords)

    if (any(diff(dif_dist_x) > 1e-10) | any(diff(dif_dist_y) > 1e-10)) {
      stop("only regular rasters or curvilinear grids supported")
    } else {
      dif_dist_x <- mean(dif_dist_x)
      dif_dist_y <- mean(dif_dist_y)
    }

  } else {
    array_mode <- TRUE
  }

  if(!is.null(geom)) {
    geom <- st_buffer(st_union(geom), buffer_dist)
    req_bbox <- st_bbox(st_transform(geom, prj))

    if(array_mode) {
      warning("Found curvilinear grid, caution, limited testing.")
      X_inds <- which(X_coords > req_bbox$xmin & X_coords < req_bbox$xmax, arr.ind = TRUE)
      Y_inds <- which(Y_coords > req_bbox$ymin & Y_coords < req_bbox$ymax, arr.ind = TRUE)

      matches <- dplyr::intersect(as.data.frame(X_inds), as.data.frame(Y_inds))

      X_inds <- seq(min(matches[, 1]), max(matches[, 1]), 1)
      Y_inds <- seq(min(matches[, 2]), max(matches[, 2]), 1)

      X_coords <- X_coords[X_inds, Y_inds]
      Y_coords <- Y_coords[X_inds, Y_inds]
    } else {

      # Grab stuff in bounding box.
      X_inds <- which(X_coords > req_bbox$xmin & X_coords < req_bbox$xmax)
      Y_inds <- which(Y_coords > req_bbox$ymin & Y_coords < req_bbox$ymax)

      if(length(X_inds) == 0 | length(Y_inds) == 0)
        stop("Data and geometry not found to intersect. Check projection?")

      X_coords <- X_coords[X_inds]
      Y_coords <- Y_coords[Y_inds]
    }

    sf_points <- construct_points(X_coords, Y_coords, X_inds, Y_inds, prj, array_mode = array_mode)

  } else {
    sf_points <- construct_points(X_coords, Y_coords, prj = prj, array_mode = array_mode)
  }

  if(array_mode) {
    sf_polygons <- get_ids(nrow(X_coords), ncol(X_coords))

    dim(sf_polygons) <- c(x = nrow(X_coords), y = ncol(X_coords))

    sf_polygons <- st_as_stars(list(sf_polygons = sf_polygons),
                               dimensions = st_dimensions(x = seq(1:nrow(X_coords)),
                                                          y = seq(1:ncol(X_coords)),
                                                          .raster = c("x", "y"),
                                                          cell_midpoints = TRUE))

    sf_polygons <- st_as_stars(sf_polygons, curvilinear = list(x = X_coords,
                                                               y = Y_coords))
  } else {

    sf_polygons <- get_ids(length(X_coords), length(Y_coords))
    dim(sf_polygons) <- c(x = length(X_coords), y = length(Y_coords))

    sf_polygons <- st_as_stars(list(sf_polygons = sf_polygons),
                               dimensions = st_dimensions(x = as.numeric(X_coords),
                                                          y = as.numeric(Y_coords),
                                                          .raster = c("x", "y"),
                                                          cell_midpoints = TRUE))
  }
  st_crs(sf_polygons) <- st_crs(prj)

  sf_polygons <- st_as_sf(sf_polygons, as_points = FALSE)

  names(sf_polygons)[1] <- "grid_ids"
  sf_polygons$grid_ids <- as.integer(sf_polygons$grid_ids)

  sf_polygons <- st_join(sf_polygons, sf_points, join = st_contains)

  sf::st_agr(sf_polygons) <- "constant"

  return(sf_polygons)
}

construct_points <- function(x, y, x_ind = NULL, y_ind = NULL, prj, array_mode) {

  if(array_mode) {
    if(is.null(x_ind)) {
      x_ind <- matrix(rep(c(1:ncol(x)), nrow(x)),
                      nrow = nrow(x), ncol = ncol(x),
                      byrow = TRUE)

      y_ind <- matrix(rep(c(1:nrow(x)), ncol(x)),
                      nrow = nrow(x), ncol = ncol(x),
                      byrow = FALSE)
    }
    x_vals <- x
    y_vals <- y
  } else {

    x_vals <- matrix(x, nrow = length(y), ncol = length(x),
                     byrow = T)
    y_vals <- matrix(y, nrow = length(y), ncol = length(x),
                     byrow = F)

    if(is.null(x_ind)) {
      x_ind <- c(1:ncol(x_vals))
      y_ind <- c(1:nrow(x_vals))
    }
  }

  if(length(x_ind) == nrow(x_vals)) {
    warning("Rows and columns flipped? Check output for valid indices.")
    X_ind <- matrix(rep(x_ind, ncol(x_vals)),
                    nrow = ncol(x_vals), ncol = nrow(x_vals),
                    byrow = TRUE)
    Y_ind <- matrix(rep(y_ind, nrow(x_vals)),
                    nrow = ncol(x_vals), ncol = nrow(x_vals),
                    byrow = FALSE)
  } else {
    X_ind <- matrix(rep(x_ind, nrow(x_vals)),
                    nrow = nrow(x_vals), ncol = ncol(x_vals),
                    byrow = TRUE)
    Y_ind <- matrix(rep(y_ind, ncol(x_vals)),
                    nrow = nrow(x_vals), ncol = ncol(x_vals),
                    byrow = FALSE)
  }

  sf_points <- st_as_sf(data.frame(x = matrix(x_vals, ncol = 1),
                                   y = matrix(y_vals, ncol = 1),
                                   X_ind = matrix(X_ind, ncol = 1),
                                   Y_ind = matrix(Y_ind, ncol = 1)),
                        coords = c("x", "y"),
                        crs = prj,
                        agr = "constant")
}

#' get data source array IDs
#'
#' Use IDs to map a NetCDF data source array onto cell geometry IDs
#'
#' @param X_dim_size size of dimension containing X coordinates (typically columns)
#' @param Y_dim_size size of dimension containing Y coordinates (typically rows)
#' @noRd
#' @examples
#' get_ids(5, 7)
#' matrix(get_ids(5,7), ncol = 1, byrow = FALSE)

get_ids <- function(X_dim_size, Y_dim_size) {
  matrix(as.numeric(seq(1, X_dim_size * Y_dim_size)),
         nrow = X_dim_size, ncol = Y_dim_size,
         byrow = FALSE)
}
