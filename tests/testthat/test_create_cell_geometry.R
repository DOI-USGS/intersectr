context("create_cell_geometry")

test_that("basic", {
  row_coords <- c(5, 4, 3, 2, 1)
  col_coords <- c(1, 2, 3, 4)
  prj <- "+init=epsg:5070"
  cells <- create_cell_geometry(col_coords, row_coords, prj)

  expect(all(c("grid_ids", "col_ind", "row_ind") %in% names(cells)))

  expect_equal(cells[cells$grid_ids == 1, ]$col_ind, 1)
  expect_equal(cells[cells$grid_ids == 1, ]$row_ind, 1)

  expect_equal(cells[cells$grid_ids == (5*4), ]$col_ind, 4)
  expect_equal(cells[cells$grid_ids == (5*4), ]$row_ind, 5)

  expect_s3_class(cells$geometry, "sfc_POLYGON")

  expect(st_crs(cells) == st_crs(prj))
})

test_that("daymet_subset", {
  col_coords <- c(1986250, 1987250, 1988250, 1989250, 1990250, 1991250, 1992250,
        1993250, 1994250, 1995250, 1996250, 1997250, 1998250, 1999250,
        2000250, 2001250, 2002250, 2003250, 2004250, 2005250, 2006250,
        2007250, 2008250, 2009250, 2010250, 2011250, 2012250, 2013250,
        2014250, 2015250, 2016250, 2017250, 2018250, 2019250, 2020250,
        2021250, 2022250, 2023250, 2024250, 2025250, 2026250, 2027250,
        2028250, 2029250, 2030250, 2031250, 2032250, 2033250, 2034250,
        2035250, 2036250, 2037250, 2038250, 2039250, 2040250, 2041250,
        2042250, 2043250, 2044250, 2045250, 2046250, 2047250, 2048250,
        2049250, 2050250, 2051250, 2052250, 2053250, 2054250, 2055250,
        2056250, 2057250, 2058250, 2059250, 2060250, 2061250, 2062250,
        2063250, 2064250)
  row_coords <- c(-503500, -504500, -505500, -506500, -507500, -508500, -509500,
         -510500, -511500, -512500, -513500, -514500, -515500, -516500,
         -517500, -518500, -519500, -520500, -521500, -522500, -523500,
         -524500, -525500, -526500, -527500, -528500, -529500, -530500,
         -531500, -532500, -533500, -534500, -535500, -536500, -537500,
         -538500, -539500, -540500, -541500, -542500, -543500, -544500,
         -545500, -546500, -547500, -548500, -549500, -550500, -551500,
         -552500, -553500, -554500, -555500, -556500, -557500)

  prj <- "+proj=lcc +lat_1=25 +lat_2=60 +x_0=0 +y_0=0 +units=m +lat_0=42.5 +lon_0=-100 +a=6378137 +f=0.00335281066474748 +pm=0 +no_defs"

  cells <- create_cell_geometry(col_coords, row_coords, prj)

  expect_equal(cells[cells$grid_ids == 1, ]$col_ind, 1)
  expect_equal(cells[cells$grid_ids == 1, ]$row_ind, 1)

  expect_equal(cells[cells$grid_ids == (length(col_coords) * length(row_coords)), ]$col_ind, length(col_coords))
  expect_equal(cells[cells$grid_ids == (length(col_coords) * length(row_coords)), ]$row_ind, length(row_coords))

  expect_s3_class(cells$geometry, "sfc_POLYGON")

  expect(st_crs(cells) == st_crs(prj))
})
