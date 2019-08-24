context("1d lat/lon")

test_that("1d lat/lon", {
  library(ncmeta)
  nc_file <- system.file("extdata/metdata.nc", package = "intersectr")

  nc_var <- nc_vars(nc_file)
  variable_name <- nc_var$name[1]

  nc_coord_vars <- nc_coord_var(nc_file, variable_name)

  x_var <- nc_coord_vars$X
  y_var <- nc_coord_vars$Y
  t_var <- nc_coord_vars$T

  nc <- RNetCDF::open.nc(nc_file)

  x <- RNetCDF::var.get.nc(nc, x_var, unpack = TRUE)
  y <- RNetCDF::var.get.nc(nc, y_var, unpack = TRUE)

  geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf")) %>%
    st_transform(5070)

  geom <- geom[5, ]

  in_prj <- "+init=epsg:4326"

  cell_geometry <- suppressWarnings(
    create_cell_geometry(x, y, in_prj, geom, 500))

  expect_true(nrow(cell_geometry) == 180)

  cell_geometry <- suppressWarnings(
    create_cell_geometry(x, y, in_prj, geom, 1000))

  expect_true(nrow(cell_geometry) == 210)

  expect_true(all(c("grid_ids", "X_ind", "Y_ind") %in% names(cell_geometry)))

  data_source_cells <- st_sf(select(cell_geometry, grid_ids))
  target_polygons <- st_sf(select(geom, CNTY_ID))

  sf::st_agr(data_source_cells) <- "constant"
  sf::st_agr(target_polygons) <- "constant"

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons)

  intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                      cell_geometry, x_var, y_var, t_var)

  expect_true(all(names(intersected) %in% c("time_stamp", "1832")))
  expect_true(nrow(intersected) == 5)

  intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                      cell_geometry, x_var, y_var, t_var,
                                      start_datetime = "2018-09-13 00:00:00",
                                      end_datetime = "2018-09-14 00:00:00")

  expect_true(nrow(intersected) == 2)

  tf <- tempfile()

  expect_warning(out_f <- write_ncdf(nc_file = tf, int_data = intersected,
                      variable_name = "test_varname", variable_units = "test_varunits"),
                 "inserting fake latitude and longitude values")

  expect_true(file.exists(tf))

  unlink(tf)

  test_file <- execute_intersection(nc_file, variable_name, area_weights,
                       cell_geometry, x_var, y_var, t_var,
                       writer_fun = write_incremental,
                       out_file = tf,
                       out_var_meta = list(name = "test",
                                           long_name = "long_test",
                                           units = "mm"))
  expect_true(file.exists(tf))

  unlink(tf)

  dap <- intersectr:::get_dap_url(min(cell_geometry$X_ind), max(cell_geometry$X_ind),
                                  min(cell_geometry$Y_ind), max(cell_geometry$Y_ind),
                                  nc_file, variable_name, x_var, y_var, t_var,
                                  start_datetime = "2018-09-13 00:00:00",
                                  end_datetime = "2018-09-14 00:00:00")

  expect_true(grepl("metdata.nc\\?lon\\[316:336\\],lat\\[73:82\\],day\\[2:3\\],precipitation_amount\\[2:3\\]\\[73:82\\]\\[316:336\\]", dap))

})
