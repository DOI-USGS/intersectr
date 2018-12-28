context("utils")

test_that("get_var_by_att", {
  nc <- system.file("extdata/metdata.nc", package = "intersectr")

  expect(find_var_by_att(nc, "coordinates") == "precipitation_amount")

  expect(length(find_var_by_att(nc, "units")) == 4)

  expect(length(find_var_by_att(nc, "units", "degrees", strict = FALSE)) == 2)

  expect(length(find_var_by_att(nc, "units", "degrees", strict = TRUE)) == 0)

  expect(length(find_var_by_att(nc, "units", "degrees_east", strict = TRUE)) == 1)

  nc <- RNetCDF::open.nc(nc)

  expect(find_var_by_att(nc, "coordinates") == "precipitation_amount")

  expect(find_var_by_att(ncmeta::nc_atts(nc), "coordinates") == "precipitation_amount")
})

test_that("get_grid_mapping", {

  nc <- system.file("extdata/metdata.nc", package = "intersectr")
  expect_warning(gm <- get_grid_mapping(nc),
                 "No variables with a grid mapping found. Defaulting to WGS84 Lon/Lat")
  expect_equal(gm, list(grid_mapping_name = "latitude_longitude",
                        semi_major_axis = 6378137,
                        inverse_flattening = 298.257223563,
                        longitude_of_prime_meridian = 0))

  nc <- system.file("extdata/daymet.nc", package = "intersectr")
  gm <- get_grid_mapping(nc)

  expect(all(list(grid_mapping_name = "lambert_conformal_conic",
                  longitude_of_central_meridian = -100,
                  latitude_of_projection_origin = 42.5,
                  false_easting = 0,
                  false_northing = 0,
                  standard_parallel = c(25, 60),
                  semi_major_axis = 6378137,
                  inverse_flattening = 298.257223563,
                  longitude_of_prime_meridian = 0) %in% gm))

  expect_is(get_grid_mapping(ncmeta::nc_atts(nc)), "list")

  expect_is(get_grid_mapping(RNetCDF::open.nc(nc)), "list")
})
