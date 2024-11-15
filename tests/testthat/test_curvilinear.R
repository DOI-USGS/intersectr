context("curvilinear")

test_that("curvilinear", {
  nc_file <- system.file("nc/test_stageiv_xyt.nc", package = "stars")

  nc_var <- nc_vars(nc_file)
  variable_name <- nc_var$name[1]

  nc_coord_vars <- nc_coord_var(nc_file, variable_name)

  nc <- rnz::open_nz(nc_file)

  geom <- st_transform(sf::read_sf(system.file("shape/nc.shp", package = "sf")), 5070)[5, ]

  cell_geometry <- suppressWarnings(
    create_cell_geometry(rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE),
                         rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE),
                         "+init=epsg:4326"))

  cell_geometry2 <- suppressWarnings(
    create_cell_geometry(rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE),
                         rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE),
                         "+init=epsg:4326", st_transform(cell_geometry, 5070)))

  expect_true(all(cell_geometry$X_ind == cell_geometry2$X_ind))

  expect_equivalent(st_coordinates(cell_geometry), st_coordinates(cell_geometry2))

  expect_equal(cell_geometry$X_ind[1:5], c(1, 2, 3, 4, 5))
  expect_equal(cell_geometry$Y_ind[1:5], c(1, 1, 1, 1, 1))

  cell_geometry <- suppressWarnings(
    create_cell_geometry(rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE),
                         rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE),
                         "+init=epsg:4326", geom))

  expect_true(nrow(cell_geometry) == 340)

  cell_geometry <- suppressWarnings(
    create_cell_geometry(rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE),
                         rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE),
                         "+init=epsg:4326", geom, 10000))
  rnz::close_nz(nc)

  expect_true(nrow(cell_geometry) == 700)

  data_source_cells <- st_sf(select(cell_geometry, grid_ids))
  target_polygons <- st_sf(select(geom, CNTY_ID))

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons)

  suppressWarnings(intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                      cell_geometry, nc_coord_vars$X, nc_coord_vars$Y, nc_coord_vars$T))

  expect_equal(intersected$`1832`[intersected$time_stamp == as.POSIXct("2018-09-14 15:00:00", tz = "UTC")],
               0.7737, tolerance = 0.001)

  nc_file <- system.file("extdata/test_stageiv_xyt_borked.nc", package = "intersectr")
  nc <- rnz::open_nz(nc_file)

  cell_geometry <- suppressWarnings(
    create_cell_geometry(rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE),
                         rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE),
                         "+init=epsg:4326", geom, 10000))
  rnz::close_nz(nc)

  data_source_cells <- st_sf(select(cell_geometry, grid_ids))

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons)

  suppressWarnings(intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                                       cell_geometry, nc_coord_vars$X,
                                                       nc_coord_vars$Y, nc_coord_vars$T))

  expect_equal(intersected$`1832`[intersected$time_stamp == as.POSIXct("2018-09-14 15:00:00", tz = "UTC")],
               0.7737, tolerance = 0.001)
})

test_that("curvilinear", {
  nc_file <- system.file("extdata/c201923412.out1_4.nc", package = "intersectr")
  variable_name <- "wvh"
  nc <- rnz::open_nz(nc_file)

  nc_coord_vars <- nc_coord_var(nc_file, variable_name)

  x <- rnz::get_var(nc, nc_coord_vars$X, unpack = TRUE)
  y <- rnz::get_var(nc, nc_coord_vars$Y, unpack = TRUE)

  cell_geometry <- suppressWarnings(
    create_cell_geometry(x, y, "+init=epsg:4326"))

  geom <- cell_geometry[530:600, ]

  data_source_cells <- st_transform(st_sf(select(cell_geometry, grid_ids)), 5070)
  target_polygons <- st_transform(st_sf(select(geom, ids = grid_ids)), 5070)

  area_weights <- calculate_area_intersection_weights(
    data_source_cells,
    target_polygons)

  suppressWarnings(intersected <- execute_intersection(nc_file, variable_name, area_weights,
                                                       cell_geometry, nc_coord_vars$X,
                                                       nc_coord_vars$Y, nc_coord_vars$T))

  expect_equal(ncol(intersected), 72)

  geom["data"] <- as.numeric(intersected[, 2:ncol(intersected)])
  geom$data[geom$data == -99999] <- NA
  geom <- st_transform(geom, 5070)

  check_row <- st_transform(geom[15, ], 4326)
  mid <- st_sfc(st_point(c(mean(st_coordinates(check_row)[1:4, 1]),
                           mean(st_coordinates(check_row)[1:4, 2]))), crs = 4326)

  x_inds <- which(abs((x - st_coordinates(mid)[1,1])) == min(abs((x - st_coordinates(mid)[1,1]))), arr.ind = TRUE)
  y_inds <- which(abs((y - st_coordinates(mid)[1,2])) == min(abs((y - st_coordinates(mid)[1,2]))), arr.ind = TRUE)

  x_ind <- unique(x_inds[which(x_inds[, 1] %in% y_inds[, 1]), 1])
  y_ind <- unique(y_inds[which(x_inds[, 2] %in% y_inds[, 2]), 2])

  val <- rnz::get_var(nc, variable_name, c(x_ind, y_ind, 1), c(1,1,1))

  expect_equal(check_row$data, as.numeric(val))

  # dat <- stars::read_ncdf(nc_file)
  # dat$wvh[dat$wvh == units::set_units(-99999, "m")] <- units::set_units(NA, "m")
  # dat <- st_transform(dat, 5070)
  # dat <- st_as_sf(dat)
  # plot(st_geometry(dat))
  # plot(dat["2019-08-22 13:00:00"], border = NA, add = TRUE)
  #
  # plot(geom["data"], add = TRUE)

})

# geom_data <- select(target_polygons, grid_ids_2)
# g_names <- as.numeric(names(intersected)[2:ncol(intersected)])
#
# plot_fun <- function(row, geom_data, intersected, g_names) {
#   geom_data <- geom_data %>%
#     dplyr::left_join(data.frame(grid_ids_2 = g_names,
#                                 poly_data = as.numeric(intersected[row, 2:ncol(intersected)]),
#                                 stringsAsFactors = FALSE),
#                      by = "grid_ids_2")
#   plot(geom_data["poly_data"], border = NA, breaks = c(0,10,20,30,40,50,60,70,80,90,100))
# }
#
# gifski::save_gif(lapply(1:nrow(intersected), plot_fun, intersected = intersected, geom_data = geom_data,
#                         g_names = g_names), gif_file = "test.gif")
