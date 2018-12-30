context("1d projected")

test_that("1d projected", {
  variable_name <- "prcp"
  nc_file <- system.file("extdata/daymet.nc", package = "intersectr")
  nc <- RNetCDF::open.nc(nc_file)
  x_var <- "x"
  y_var <- "y"
  t_var <- "time"

  x <- RNetCDF::var.get.nc(nc, x_var) * 1000
  y <- RNetCDF::var.get.nc(nc, y_var) * 1000

  geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf")) %>%
    st_transform(5070) %>%
    dplyr::filter(CNTY_ID == 2156)

  grid_mapping <- get_grid_mapping(nc)

  in_prj <- ncdfgeom::get_prj(grid_mapping)

  cell_geometry <- suppressWarnings(
    create_cell_geometry(x, y, in_prj, geom, 1000))

  expect(nrow(cell_geometry) == 4212)
  expect(all(c("grid_ids", "col_ind", "row_ind") %in% names(cell_geometry)))

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

  expect(all(names(intersected) %in% c("time_stamp", "2156")))
  expect(nrow(intersected) == 5)

  intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                      cell_geometry, x_var, y_var, t_var,
                                      start_datetime = "1999-09-14 00:00:00",
                                      end_datetime = "1999-09-16 00:00:00")

  expect(nrow(intersected) == 2)

  expect_equal(intersected$`2156`, c(0, 10.3973), tolerance = 0.001)
})
