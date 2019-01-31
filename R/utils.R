#' Get Timerange
#' @description Helper function to get start and end time stamps from a NetCDF source.
#' @param nc_file character path or url to a NetCDF source
#' @param t_var character name of time variable
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
