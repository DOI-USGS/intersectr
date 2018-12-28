#' Find NetCDF Variable by attribute
#' @description Given an attribute name and potentially a value,
#' searches for and returns variables with the desired attribute.
#'
#' @param x data.frame returned by ncmeta::nc_atts, open NetCDF object,
#' or character file path or url to be opened with RNetCDF::open.nc
#' @param attribute character the attribute name to search for variables with
#' @param value character defaults to any only return variables that have the
#' attribute with the given value
#' @param strict boolean if TRUE, only exact matches of value will be returned
#'
#' @importFrom RNetCDF open.nc close.nc
#' @importFrom ncmeta nc_atts
#' @export
#'
#' @examples
#' nc <- system.file("extdata/metdata.nc", package = "intersectr")
#'
#' find_var_by_att(nc, "coordinates")
#'
#' find_var_by_att(nc, "units")
#'
#' find_var_by_att(nc, "units", "degrees", strict = FALSE)
#'
#' find_var_by_att(nc, "units", "degrees", strict = TRUE)
#'
#' find_var_by_att(nc, "units", "degrees_east", strict = TRUE)
#'
find_var_by_att <- function(x, attribute, value = ".*", strict = TRUE) {

  open_nc <- FALSE
  if (is.character(x)) {
    x <- open.nc(x)
    open_nc <- TRUE
  }

  if (inherits(x, "NetCDF")) {
    atts <- nc_atts(x)
  } else if (inherits(x, "data.frame")) {
    atts <- x
  }

  if (strict) value <- paste0("^", value, "$")

  atts <- atts[atts$name == attribute, ]
  atts <- atts[grepl(value, atts$value), ]

  if (open_nc) close.nc(x)

  return(atts$variable)
}

#' Get Grid Mapping
#'
#' @description Get the grid mapping from a NetCDF file
#' @param x data.frame returned by ncmeta::nc_atts, open NetCDF object,
#' or character file path or url to be opened with RNetCDF::open.nc
#' @export
#' @examples
#' nc <- system.file("extdata/metdata.nc", package = "intersectr")
#' get_grid_mapping(nc)
#'
#' nc <- system.file("extdata/daymet.nc", package = "intersectr")
#' get_grid_mapping(nc)
#'
get_grid_mapping <- function(x) {

  open_nc <- FALSE
  if (is.character(x)) {
    x <- open.nc(x)
    open_nc <- TRUE
  }

  if (inherits(x, "NetCDF")) {
    atts <- nc_atts(x)
  } else if (inherits(x, "data.frame")) {
    atts <- x
  }

  if (open_nc) close.nc(x)

  gm_att <- "grid_mapping"
  grid_mapping_vars <- find_var_by_att(atts, gm_att)

  if (length(grid_mapping_vars) == 0) {
    warning(paste("No variables with a grid mapping found.\n",
                  "Defaulting to WGS84 Lon/Lat"))
    return(list(grid_mapping_name = "latitude_longitude",
                semi_major_axis = 6378137,
                inverse_flattening = 298.257223563,
                longitude_of_prime_meridian = 0))
  }

  grid_mapping_var <- unique(atts$variable[atts$name == "grid_mapping_name"])

  if (length(grid_mapping_var) > 1) {
    stop("Found more than one grid mapping variable. Only one is supported.")
  }

  grid_mapping <- atts[atts$variable == grid_mapping_var, ]

  return(stats::setNames(grid_mapping$value, grid_mapping$name))
}
