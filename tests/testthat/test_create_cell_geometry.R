context("create_cell_geometry")

test_that("area_intersection", {
  row_coords <- c(5, 4, 3, 2, 1)
  col_coords <- c(1, 2, 3, 4)
  prj <- "+init=epsg:5070"
  cells <- create_cell_geometry(col_coords, row_coords, prj)

  expect(all(c("grid_ids", "x_ind", "y_ind") %in% names(cells)))

  expect_equal(cells[cells$grid_ids == 1, ]$x_ind, 1)
  expect_equal(cells[cells$grid_ids == 1, ]$y_ind, 1)

  expect_equal(cells[cells$grid_ids == (5*4), ]$x_ind, 4)
  expect_equal(cells[cells$grid_ids == (5*4), ]$y_ind, 5)

  expect_s3_class(cells$geometry, "sfc_POLYGON")

  expect(st_crs(cells) == st_crs(prj))
})
