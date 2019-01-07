#' Execute Intersection
#'
#' @description Traverses many time steps of data applying
#' intersection weights for each time step.
#'
#' @param nc_file character path or OPeNDAP URL to NetCDF source
#' @param variable_name character name of variable in NetCDF source to access
#' @param intersection_weights data.frame as output by
#' \code{\link{calculate_area_intersection_weights}}
#' @param cell_geometry sf data.frame with col_ind and row_ind columns as output by
#' \code{\link{create_cell_geometry}}
#' @param x_var character variable with X axis data
#' @param y_var character variable with Y axis data
#' @param t_var character variable with T axis data
#' @param start_datetime character or POSIX character format is
#' \code{\link{strptime}} default e.g. "2010-10-10 10:10:10"
#' @param end_datetime character or POSIX
#' @importFrom dplyr right_join mutate summarise group_by ungroup
#' @importFrom RNetCDF open.nc dim.inq.nc var.inq.nc att.get.nc var.get.nc utcal.nc
#' @export
#' @examples
#'
#' library(dplyr)
#' library(sf)
#'
#' nc_file <- system.file("extdata/metdata.nc", package = "intersectr")
#' nc <- RNetCDF::open.nc(nc_file)
#' ncmeta::nc_vars(nc)
#' variable_name <- "precipitation_amount"
#' cv <- ncmeta::nc_coord_var(nc, variable_name)
#'
#' x <- RNetCDF::var.get.nc(nc, cv$X, unpack = TRUE)
#' y <- RNetCDF::var.get.nc(nc, cv$Y, unpack = TRUE)
#'
#' prj <- ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc))
#' out_prj <- "+init=epsg:5070"
#'
#' geom <- read_sf(system.file("shape/nc.shp", package = "sf")) %>%
#'   st_transform(out_prj)
#'
#' cell_geometry <- create_cell_geometry(x, y, prj, geom, 0)
#'
#' data_source_cells <- st_sf(select(cell_geometry, grid_ids))
#' target_polygons <- st_sf(select(geom, CNTY_ID))
#'
#' st_agr(data_source_cells) <- "constant"
#' st_agr(target_polygons) <- "constant"
#'
#' area_weights <- calculate_area_intersection_weights(
#'   data_source_cells,
#'   target_polygons)
#'
#' intersected <- execute_intersection(nc_file, variable_name, area_weights,
#'                                     cell_geometry, cv$X, cv$Y, cv$T)
#'
execute_intersection <- function(nc_file,
                                 variable_name,
                                 intersection_weights,
                                 cell_geometry,
                                 x_var, y_var, t_var,
                                 start_datetime = NULL,
                                 end_datetime = NULL) {

  names(intersection_weights) <- c("grid_ids", "poly_id", "w")

  nc <- open.nc(nc_file)
  nc_var_info <- var.inq.nc(nc, variable_name)

  if (nc_var_info$ndims != 3) stop("only 3d variables are supported")

  col_var_info <- var.inq.nc(nc, x_var)
  row_var_info <- var.inq.nc(nc, y_var)
  t_var_info <- var.inq.nc(nc, t_var)

  if (col_var_info$ndims > 1 | row_var_info$ndims > 1 | t_var_info$ndims > 1)
    stop("only 1d coordinate variables supported")

  time_steps <- utcal.nc(att.get.nc(nc, t_var_info$name, "units"),
                         var.get.nc(nc, t_var_info$name, unpack = TRUE),
                         "c")
  time_steps <- data.frame(time_stamp = time_steps)

  col_inds <- seq(min(cell_geometry$col_ind), max(cell_geometry$col_ind), 1)
  row_inds <- seq(min(cell_geometry$row_ind), max(cell_geometry$row_ind), 1)

  ids <- get_ids(length(col_inds), length(row_inds))

  if (is.null(start_datetime) & is.null(end_datetime)) {

    out_nrows <- nrow(time_steps)

    t_inds <- seq_len(out_nrows)

  } else {

    if (is.null(start_datetime)) start_datetime <- time_steps[1]
    if (is.character(start_datetime)) {
      start_datetime <- strptime(start_datetime,
                                 format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    }
    if (is.character(end_datetime)) {
      end_datetime <- strptime(end_datetime,
                               format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
    }

    t_inds <- time_steps$time_stamp >= start_datetime &
      time_steps$time_stamp <= end_datetime

    time_steps <- filter(time_steps, t_inds)

    t_inds <- which(t_inds)

    out_nrows <- length(t_inds)
  }
    out_ncols <- length(unique(intersection_weights[, 2][[1]]))

    out <- matrix(nrow = out_nrows, ncol = out_ncols)

    dimid_order <- match(nc_var_info$dimids,
                         c(col_var_info$dimids,
                           row_var_info$dimids,
                           t_var_info$dimids))

    for (i in 1:out_nrows) {
      try_backoff({
        i_data <- var.get.nc(nc, variable_name,
                             start = c(min(col_inds),
                                       min(row_inds), t_inds[i])[dimid_order],
                             count = c(length(col_inds),
                                       length(row_inds), 1)[dimid_order],
                             unpack = TRUE)

        i_data <- data.frame(d = matrix(i_data,
                                        ncol = 1,
                                        byrow = TRUE),
                             grid_ids = matrix(ids, ncol = 1))

        i_data <- right_join(i_data,
                             intersection_weights,
                             by = c("grid_ids"))
        i_data <- group_by(i_data, poly_id)
        i_data <- summarise(i_data,
                            d = (sum( (d * w), na.rm = T) /
                                   sum(w, na.rm = T)))
          i_data <- ungroup(i_data)

        out[i, ] <- i_data$d
      })
    }

  out <- data.frame(out)

  out_names <- as.character(unique(i_data$poly_id))

  out <- select(out, which(!is.na(out_names)))

  out_names <- out_names[which(!is.na(out_names))]

  names(out) <- out_names

  return(cbind(time_steps, out))
}

# Taken from: https://github.com/ramhiser/retry/blob/master/R/try-backoff.r
# Un exported from retry package so including directly.
#' Try/catch with exponential backoff
#'
#' Attempts the expression in \code{expr} up to the number of tries specified
#' in \code{max_attempts}. Each time a failure results, the functions sleeps
#' for a random amount of time before re-attempting the expression. The upper
#' bound of the backoff increases exponentially after each failure.
#'
#' For details on exponential backoff, see:
#' \url{http://en.wikipedia.org/wiki/Exponential_backoff}
#'
#' @param expr an R expression to try.
#' @param silent logical: should the report of error messages be suppressed?
#' @param max_tries the maximum number of times to attempt the expression
#' \code{expr}
#' @param verbose logical: Should detailed messages be reported regarding each
#' attempt? Default: no.
#' @return the value of the expression in \code{expr}. If the final attempt was
#' a failure, the objected returned will be of class try-error".
#' @importFrom stats runif
#' @examples
#' # Example that will never succeed.
#' try_backoff(log("a"), verbose=TRUE, max_attempts=5)
#' @noRd
try_backoff <- function(expr, silent=FALSE, max_attempts=10, verbose=FALSE) {
  for (attempt_i in seq_len(max_attempts)) {
    results <- try(expr = expr, silent = silent)
    if (class(results) == "try-error") {
      backoff <- runif (n = 1, min = 0, max = 2 ^ attempt_i - 1)
      if (verbose) {
        message("Backing off for ", backoff, " seconds.")
      }
      Sys.sleep(backoff)
    } else {
      if (verbose) {
        message("Succeeded after ", attempt_i, " attempts.")
      }
      break
    }
  }
  results
}
