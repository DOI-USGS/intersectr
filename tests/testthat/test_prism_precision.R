context("precise execute intersection")

# There are three sample files with data for 1999.
# This function is used for all three. See calls at the bottom.

test_precise <- function(nc_file, variable_name) {
  nc_coord_vars <- ncmeta::nc_coord_var(nc_file, variable_name)

  geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))[5, ]

  suppressWarnings(nc_prj <- ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc_file)))

  nc <- RNetCDF::open.nc(nc_file)
  X_coords <- RNetCDF::var.get.nc(nc, nc_coord_vars$X, unpack = TRUE)
  X_coords <- seq(from = X_coords[1],
                  to = X_coords[length(X_coords)],
                  along.with = X_coords)

  Y_coords <- RNetCDF::var.get.nc(nc, nc_coord_vars$Y, unpack = TRUE)
  Y_coords <- seq(from = Y_coords[1],
                  to = Y_coords[length(Y_coords)],
                  along.with = Y_coords)

  suppressWarnings(cell_geometry <- create_cell_geometry(X_coords = X_coords,
                                                         Y_coords = Y_coords,
                                                         prj = nc_prj,
                                                         geom = geom,
                                                         buffer_dist = 0.05))

  geom_2 <- cell_geometry[20, ]
  names(geom_2)[1] <- "grid_ids2"

  data_source_cells <- sf::st_sf(dplyr::select(cell_geometry, grid_ids))
  target_polygons <- sf::st_sf(dplyr::select(geom_2, grid_ids2))
  sf::st_agr(data_source_cells) <- "constant"
  sf::st_agr(target_polygons) <- "constant"

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons)

  intersected <- execute_intersection(nc_file,
                                      variable_name,
                                      area_weights,
                                      cell_geometry,
                                      nc_coord_vars$X,
                                      nc_coord_vars$Y,
                                      nc_coord_vars$T,
                                      start_datetime = "1999-01-01 00:00:00",
                                      end_datetime = "1999-12-31 23:59:59")

  dates <- RNetCDF::var.get.nc(nc, "time", unpack = TRUE)
  date_units <- RNetCDF::att.get.nc(nc, "time", "units")
  dates <- RNetCDF::utcal.nc(date_units, dates, type = "c")
  date_ind <- which(as.character(dates) == "1999-01-01")

  nc_var_info <- RNetCDF::var.inq.nc(nc, variable_name)

  X_var_info <- RNetCDF::var.inq.nc(nc, nc_coord_vars$X)
  Y_var_info <- RNetCDF::var.inq.nc(nc, nc_coord_vars$Y)
  T_var_info <- RNetCDF::var.inq.nc(nc, nc_coord_vars$T)

  dimid_order <- match(nc_var_info$dimids,
                       c(X_var_info$dimids,
                         Y_var_info$dimids,
                         T_var_info$dimids))

  vals <- RNetCDF::var.get.nc(nc, variable_name,
                              start = c(geom_2$X_ind, geom_2$Y_ind, date_ind)[dimid_order],
                              count = c(1, 1, 12)[dimid_order], unpack = TRUE)

  test_that("values come back 1:1 for precise grid cell with time filter", {
    expect_equal(as.numeric(vals), as.numeric(intersected[[as.character(geom_2$grid_ids2)]]))
  })

  intersected <- execute_intersection(nc_file,
                                      variable_name,
                                      area_weights,
                                      cell_geometry,
                                      nc_coord_vars$X,
                                      nc_coord_vars$Y,
                                      nc_coord_vars$T)

  test_that("values come back 1:1 for precise grid cell with time filter", {
    expect_equal(as.numeric(vals), as.numeric(intersected[[as.character(geom_2$grid_ids2)]]))
  })
}

test_precise(system.file("extdata/bcsd_obs_1999_borked.nc", package = "intersectr"),
             "pr")

test_precise(system.file("extdata/bcsd_obs_1999.nc", package = "intersectr"),
             "pr")

test_precise(system.file("extdata/prism.nc", package = "intersectr"),
             "ppt")
