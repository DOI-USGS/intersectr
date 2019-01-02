#todo... just leaving code here.
nc_file <- "https://cida.usgs.gov/thredds/dodsC/prism_v2"
nc <- RNetCDF::open.nc(nc_file)

(nc_var <- ncmeta::nc_vars(nc))

(nc_coord_vars <- ncmeta::nc_coord_var(nc))

variable_name <- "ppt"
(nc_coord_vars <- nc_coord_vars[nc_coord_vars$variable == variable_name, ])

(geom <- sf::st_transform(sf::read_sf(system.file("shape/nc.shp",
                                                 package = "sf")),
                         5070))

(nc_prj <- ncdfgeom::get_prj(intersectr::get_grid_mapping(nc)))

# Small diffs
col_coords <- RNetCDF::var.get.nc(nc, nc_coord_vars$X)
col_coords <- seq(from = col_coords[1],
                  to = col_coords[length(col_coords)],
                  along.with = col_coords)

row_coords <- RNetCDF::var.get.nc(nc, nc_coord_vars$Y)
row_coords <- seq(from = row_coords[1],
                  to = row_coords[length(row_coords)],
                  along.with = row_coords)

(cell_geometry <-
  create_cell_geometry(col_coords = col_coords,
                       row_coords = row_coords,
                       prj = nc_prj,
                       geom = geom,
                       buffer_dist = 1000))

data_source_cells <- sf::st_sf(dplyr::select(cell_geometry, grid_ids))
target_polygons <- sf::st_sf(dplyr::select(geom, CNTY_ID))
sf::st_agr(data_source_cells) <- "constant"
sf::st_agr(target_polygons) <- "constant"

(area_weights <- calculate_area_intersection_weights(
  data_source_cells,
  target_polygons))

(intersected <- execute_intersection(nc_file,
                                     variable_name,
                                     area_weights,
                                     cell_geometry,
                                     nc_coord_vars$X,
                                     nc_coord_vars$Y,
                                     nc_coord_vars$T))

i <- 1
breaks <- c(0, 10, 20, 40, 60, 80, 100, 150, 300, 600, 800)

gifski::save_gif({
  for(i in 1:nrow(intersected)) {
    geom_data <- select(geom, CNTY_ID) %>%
      left_join(data.frame(CNTY_ID = as.numeric(names(intersected)[2:ncol(intersected)]),
                           poly_data = as.numeric(intersected[i, 2:ncol(intersected)]) / 100),
                by = "CNTY_ID")

    plot(geom_data["poly_data"], breaks = breaks, main = "Demo Prism Animation (ppt mm/month)")
  }
}, gif_file = "test.gif", width = 1024, height = 768, delay = 0.1, progress = TRUE)

plot(intersected$time_stamp, intersected$`1825`, col = NA)
lines(intersected$time_stamp, intersected$`1825`)

