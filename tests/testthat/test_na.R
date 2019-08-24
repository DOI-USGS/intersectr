library(stars)
library(sf)
library(dplyr)

nc_file <- system.file("extdata/prism.nc", package = "intersectr")
stars_dat <- stars::read_ncdf(nc_file)

variable_name <- "ppt"
nc <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))

plot(slice(stars_dat, index = 1, along = "time"), reset = FALSE)

geom <- nc[5,]

test_0 <- test_na(nc_file, variable_name, geom)

st_geometry(geom) <- st_geometry(geom) + c(1.5,0)
st_crs(geom) <- st_crs(nc)

plot(st_geometry(geom), add = TRUE, reset = FALSE, col = NA, border = 'red')

test_1 <- test_na(nc_file, variable_name, geom)

expect_true(test_1$`1832` > 100)

geom <- nc[5,]
st_geometry(geom) <- st_geometry(geom) + c(2,0)
st_crs(geom) <- st_crs(nc)

plot(st_geometry(geom), add = TRUE, reset = FALSE, col = NA, border = 'red')

test_2 <- test_na(nc_file, variable_name, geom)

expect_true(test_2$`1832` > 100)

