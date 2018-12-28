#' Area Weighted Intersection
#' @description Returns the fractional percent of each
#' feature in y that is covered by each intersecting feature
#' in x. These can be used as the weights in an area-weighted
#' mean overlay analysis where y is the data source and area-
#' weighted means are being generated for x.
#'
#' @details
#' THIS WILL NOT WORK WITH SELF INTERSECTIONS!!!!!
#'
#' @param x sf data.frame including one geometry column and one identifier column
#' @param y sf data.frame including one geometry column and one identifier column
#' @return data.frame containing fraction of each feature in y that is
#' covered by each feature in x.
#' s
#' @examples
#' b1 = sf::st_polygon(list(rbind(c(-1,-1), c(1,-1),
#'                            c(1,1), c(-1,1),
#'                            c(-1,-1))))
#' b2 = b1 + 2
#' b3 = b1 + c(-0.2, 2)
#' b4 = b1 + c(2.2, 0)
#' b = sf::st_sfc(b1, b2, b3, b4)
#' a1 = b1 * 0.8
#' a2 = a1 + c(1, 2)
#' a3 = a1 + c(-1, 2)
#' a = sf::st_sfc(a1,a2,a3)
#' plot(b, border = 'red')
#' plot(a, border = 'green', add = TRUE)
#'
#' sf::st_crs(b) <- sf::st_crs(a) <- sf::st_crs(5070)
#'
#' b <- sf::st_sf(b, data.frame(idb = c(1, 2, 3, 4)))
#' a <- sf::st_sf(a, data.frame(ida = c(1, 2, 3)))
#'
#' sf::st_agr(a) <- sf::st_agr(b) <- "constant"
#'
#' a_b <- calculate_area_intersection_weights(a, b)
#' b_a <- calculate_area_intersection_weights(b, a)
#'
#' @export
#' @importFrom sf st_intersection st_set_geometry st_area st_crs
#' @importFrom dplyr mutate group_by right_join select ungroup

calculate_area_intersection_weights <- function(x, y) {

  if(st_crs(x) != st_crs(y)) {
    x <- st_transform(x, st_crs(y))
  }

  # Standard evaluation is for chumps.
  id_x <- names(x)[names(x) != attr(x, "sf_column")]
  id_y <- names(y)[names(y) != attr(y, "sf_column")]

  if(length(id_x) != 1 | length(id_y) != 1)
    stop("x and y must have one and only one non-geometry column")

  names(x)[names(x) == id_x] <- "varx"
  names(y)[names(y) == id_y] <- "vary"

  # Get all parts and calculate their individual area
  int <- st_intersection(x, y)
  int <- mutate(int,
                          part_area = as.numeric(st_area(int)))
  int <- group_by(int, varx) # Allow sum parts over each x.
  int <- st_set_geometry(int, NULL)

  # Get the area of x.
  x_area <- st_set_geometry(mutate(x, x_area = as.numeric(st_area(x))), NULL)

  # Join the intersecting area with the all the parts.
  int <- right_join(int, x_area, by = "varx")

  # Join total x area and calculate percent for each sum of intersecting y.
  int <- mutate(int, w = part_area/x_area)
  int <- select(int, varx, vary, w)
  int <- ungroup(int)

  out <- stats::setNames(int, c(id_x, id_y, "w"))

  return(out)
}
