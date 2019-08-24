context("1d projected")

test_that("1d projected", {
  library(ncmeta)
  nc_file <- system.file("extdata/daymet.nc", package = "intersectr")

  nc_var <- nc_vars(nc_file)
  variable_name <- nc_var$name[1]

  suppressWarnings(nc_coord_vars <- nc_coord_var(nc_file, variable_name))

  x_var <- nc_coord_vars$X
  y_var <- nc_coord_vars$Y
  t_var <- nc_coord_vars$T

  nc <- RNetCDF::open.nc(nc_file)

  x <- RNetCDF::var.get.nc(nc, x_var, unpack = TRUE) * 1000
  y <- RNetCDF::var.get.nc(nc, y_var, unpack = TRUE) * 1000

  geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf")) %>%
    st_transform(5070) %>%
    dplyr::filter(CNTY_ID == 2156)

  in_prj <- ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc))

  cell_geometry <- suppressWarnings(
    create_cell_geometry(x, y, in_prj, geom, 1000))

  expect_true(nrow(cell_geometry) == 3626)
  expect_true(all(c("grid_ids", "X_ind", "Y_ind") %in% names(cell_geometry)))

  data_source_cells <- st_sf(select(cell_geometry, grid_ids))
  target_polygons <- st_sf(select(geom, CNTY_ID))

  sf::st_agr(data_source_cells) <- "constant"
  sf::st_agr(target_polygons) <- "constant"

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons)

  intersected <- execute_intersection(nc_file,
                                      variable_name,
                                      area_weights,
                                      cell_geometry,
                                      x_var, y_var, t_var)

  expect_true(all(names(intersected) %in% c("time_stamp", "2156")))
  expect_true(nrow(intersected) == 5)

  intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                      cell_geometry, x_var, y_var, t_var,
                                      start_datetime = "1999-09-14 00:00:00",
                                      end_datetime = "1999-09-16 00:00:00")

  expect_true(nrow(intersected) == 3)
})

test_that("1d projected", {
  # library(ncmeta)
  #
  # nc_file <- "https://thredds.daac.ornl.gov/thredds-daymet/dodsC/daymet-v3-agg/na.ncml"
  # variable_name <- "prcp"
  #
  # suppressWarnings(nc_coord_vars <- nc_coord_var(nc_file, variable_name))
  #
  # nc <- RNetCDF::open.nc(nc_file)
  # nc_prj <- ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc_file))
  # X_coords <- RNetCDF::var.get.nc(nc, nc_coord_vars$X[2], unpack = TRUE)
  # Y_coords <- RNetCDF::var.get.nc(nc, nc_coord_vars$Y[2], unpack = TRUE)

  test_data <- readRDS("data/daymet_cell_test.rds")

  geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf")) %>%
    st_transform(5070) %>%
    dplyr::filter(CNTY_ID == 2156)

  cell_geometry <-
    create_cell_geometry(X_coords = test_data$X_coords,
                         Y_coords = test_data$Y_coords,
                         prj = test_data$nc_prj,
                         geom = geom)

  expect_true(nrow(cell_geometry) == 3384)

})
