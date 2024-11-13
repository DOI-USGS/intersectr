suppressMessages({
  library("sf")
  library("dplyr")
  library("RNetCDF")
  library("ncmeta")
  library("ncdfgeom")
  })

# Code to run the writer function in isolation
# file_handle <- "test.nc"
# var_meta <- list(name = "test", units = "mm", long_name = "test_long", ids = c("1", "2", "3"))
# test_data <- data.frame(times = out_table$time_stamp, "1" = c(1,2,3,4,5), "2" = c(6,7,8,9,10), "3" = c(11,12,13,14,15))
# names(test_data) <- c("times", "1", "2", "3")
#
# size <- c(nrow(test_data), (ncol(test_data) - 1))
#
# nc <- write_incremental(file_handle, 0, size, var_meta)
# for(i in 1:size[1]){
# write_incremental(nc, i, size, var_meta, unname(test_data[i, ]))
# }
# write_incremental(nc, -1, size, var_meta)
write_incremental <- function(file_handle,
                              step,
                              size,
                              var_meta,
                              timestep_data = NULL) {

  date_origin <- "days since 1900-01-01"

  if(step == 0) {
    nc <- RNetCDF::create.nc(file_handle, clobber = FALSE, large = TRUE, prefill = TRUE)
    RNetCDF::dim.def.nc(nc, "time", size[1], unlim = FALSE)
    RNetCDF::dim.def.nc(nc, "hru", size[2], unlim = FALSE)

    RNetCDF::var.def.nc(nc, "time", "NC_INT", "time")
    RNetCDF::att.put.nc(nc, "time", "calendar", "NC_CHAR", "standard")
    RNetCDF::att.put.nc(nc, "time", "units", "NC_CHAR", date_origin)

    RNetCDF::var.def.nc(nc, "hru", "NC_INT", "hru")
    RNetCDF::att.put.nc(nc, "hru", "long_name", "NC_CHAR", "Hydrologic Response Unit ID (HRU)")
    RNetCDF::var.put.nc(nc, "hru", as.numeric(var_meta$ids))

    RNetCDF::var.def.nc(nc, var_meta$name, "NC_FLOAT", dimensions = c("hru", "time"))
    RNetCDF::att.put.nc(nc, var_meta$name, "units", "NC_CHAR", var_meta$units)
    RNetCDF::att.put.nc(nc, var_meta$name, "long_name", "NC_CHAR", var_meta$long_name)
    RNetCDF::att.put.nc(nc, var_meta$name, "missing_value", "NC_FLOAT", -9999)

    RNetCDF::close.nc(nc)

    nc <- RNetCDF::open.nc(file_handle, write = TRUE, prefill = FALSE)

    return(nc)

  } else if(step == -1) {
    RNetCDF::close.nc(file_handle)
    return(invisible(NULL))
  } else {
    time_step <- RNetCDF::utinvcal.nc(date_origin, timestep_data[[1]])
    RNetCDF::var.put.nc(file_handle, "time", time_step, step, 1)

    RNetCDF::var.put.nc(file_handle, var_meta$name,
                        timestep_data[[2]],
                        start = c(NA, step), count = c(NA, 1))
    return(invisible(file_handle))
  }
}


test_na <- function(nc_file, variable_name, geom) {
  nc_coord_vars <- ncmeta::nc_coord_var(nc_file, variable_name)

  suppressWarnings(nc_prj <- ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc_file)))

  nc <- rnz::open_nz(nc_file)
  X_coords <- rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE)
  X_coords <- seq(from = X_coords[1],
                  to = X_coords[length(X_coords)],
                  along.with = X_coords)

  Y_coords <- rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE)
  Y_coords <- seq(from = Y_coords[1],
                  to = Y_coords[length(Y_coords)],
                  along.with = Y_coords)

  suppressWarnings(cell_geometry <- create_cell_geometry(X_coords = X_coords,
                                                         Y_coords = Y_coords,
                                                         prj = nc_prj,
                                                         geom = geom,
                                                         buffer_dist = 0.05))
  data_source_cells <- sf::st_sf(dplyr::select(cell_geometry, grid_ids))
  target_polygons <- sf::st_sf(dplyr::select(geom, CNTY_ID))

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons,
    allow_lonlat = TRUE)

  intersected <- execute_intersection(nc_file,
                                      variable_name,
                                      area_weights,
                                      cell_geometry,
                                      nc_coord_vars$X,
                                      nc_coord_vars$Y,
                                      nc_coord_vars$T,
                                      start_datetime = "1999-01-01 00:00:00",
                                      end_datetime = "1999-01-01 00:00:00")

  return(intersected)

}
