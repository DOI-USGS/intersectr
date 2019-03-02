#' Get Timerange
#' @description Helper function to get start and end time stamps from a NetCDF source.
#' @param nc_file character path or url to a NetCDF source
#' @param t_var character name of time variable
#' @export
#' @importFrom RNetCDF open.nc var.inq.nc utcal.nc close.nc att.get.nc var.get.nc
#'
#' @examples
#' nc_file <- system.file("extdata/metdata.nc", package = "intersectr")
#' get_timerange(nc_file, ncmeta::nc_coord_var(nc_file, "precipitation_amount")$T)

get_timerange <- function(nc_file, t_var) {
  nc <- open.nc(nc_file)
  T_var_info <- var.inq.nc(nc, t_var)
  time_steps <- utcal.nc(att.get.nc(nc, T_var_info$name, "units"),
                         var.get.nc(nc, T_var_info$name, unpack = TRUE),
                         "c")
  close.nc(nc)
  return(list(start = time_steps[1], end = time_steps[length(time_steps)]))
}

#' Write output to NetCDF
#' @description Helper function to write output to NetCDF
#' @param nc_file character path or url to a NetCDF source
#' @param int_data data.frame as output by execute_intersection
#' @param variable_name character name for data variable
#' @param variable_units character units for data variable
#' @param variable_long_name character optional long name for data variable
#' @param lats numeric vector giving representative latitude for each time series
#' @param lons numeric vector giving representative longitude for each time series
#' @export

write_ncdf <- function(nc_file, int_data,
                       variable_name, variable_units, variable_long_name = NULL,
                       lats = NULL, lons = NULL) {
  times <- int_data[, 1]
  int_data <- int_data[, 2:ncol(int_data), drop = FALSE]

  instance_names <- names(int_data)

  if(is.null(lats) | is.null(lons)) {
    warning("inserting fake latitude and longitude values")
    lats <- rep(0, ncol(int_data))
    lons <- rep(0, ncol(int_data))
  }

  if(is.null(variable_long_name)) variable_long_name <- ""

  out_file <- ncdfgeom::write_timeseries_dsg(nc_file = nc_file,
                                 instance_names = instance_names,
                                 lats = lats,
                                 lons = lons,
                                 times = times,
                                 data = int_data,
                                 data_unit = variable_units,
                                 data_metadata = list(name = variable_name,
                                                      long_name = variable_long_name))
  return(out_file)
}
