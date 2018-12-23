context("create_cell_geometry")

test_that("1d lat/lon", {
  variable_name <- "precipitation_amount"
  nc <- RNetCDF::open.nc(system.file("extdata/metdata.nc", package = "intersecter"))
  x_var <- "lon"
  y_var <- "lat"

  x <- RNetCDF::var.get.nc(nc, x_var)
  y <- RNetCDF::var.get.nc(nc, y_var)

  geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))

  in_prj <- "+init=epsg:4326"
  out_prj <- "+init=epsg:5070"

  cell_geometry <- create_cell_geometry(x, y, in_prj, out_prj, geom, 1000)

  expect(nrow(cell_geometry) == 15361)
  expect(all(c("grid_id", "x_ind", "y_ind") %in% names(cell_geometry)))

  })
