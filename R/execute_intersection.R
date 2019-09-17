#' Execute Intersection
#'
#' @description Traverses many time steps of data applying
#' intersection weights for each time step.
#'
#' @param nc_file character path or OPeNDAP URL to NetCDF source
#' @param variable_name character name of variable in NetCDF source to access
#' @param intersection_weights data.frame as output by
#' \code{\link{calculate_area_intersection_weights}}
#' @param cell_geometry sf data.frame with X_ind and Y_ind columns as output by
#' \code{\link{create_cell_geometry}}
#' @param x_var character variable with X axis data
#' @param y_var character variable with Y axis data
#' @param t_var character variable with T axis data
#' @param start_datetime character or POSIX character format is
#' \code{\link{strptime}} default e.g. "2010-10-10 10:10:10"
#' @param end_datetime character or POSIX
#' @param status boolean print status if TRUE
#' @param writer_fun a function that implements the writer API described in details
#' @param out_file character path for output file to be used with writer_fun
#' @param out_var_meta list passed to writer_fun with ids
#' @importFrom RNetCDF open.nc dim.inq.nc var.inq.nc att.get.nc var.get.nc utcal.nc
#' @importFrom data.table data.table
#' @details The writer_fun input can be a function passed in that supports a basic API.
#' The function should the following inputs:
#' \enumerate{
#'   \item{file_handle: character or otherwise see description}
#'   \item{step: the step to write. Special values are 0 and -1 for
#'   initialization and closing respectively}
#'   \item{size: a length 2 integer vector giving rows and cols of the output
#'   e.g. c(100, 10) for 100 time steps and 10 columns}
#'   \item{var_meta: a named list containing variable name, long_name, units, and
#'   ids. e.g. list(name = "prcp", long_name = "precpiptation in mm", units = "mm",
#'   ids = c("1", "2", "3")) the first three values of this list are taken from the
#'   out_var_meta input.}
#'   \item{timestep_data: a one row data.frame with time in the first column and
#'   output data in 2:ncol(timestep_data)}
#' }
#' The function will be called once with step = 0. This should be used to initialize
#' the output.
#'
#' Initialization should return a file path or open file handle, such as an RNetCDF object.
#'
#' The function will be called once for each timestep with step from 1:timesteps Each
#' row should be written as it is provided. This should return the file path or open
#' file handle.
#'
#' The function will be called once with step = -1 to indicate the file can be closed. The
#' final return value is discarded and the out_file input is returned.
#'
#' @export
#'
#' @examples
#' nc_file <- system.file("extdata/metdata.nc", package = "intersectr")
#' nc <- RNetCDF::open.nc(nc_file)
#'
#' variable_name <- "precipitation_amount"
#' (cv <- ncmeta::nc_coord_var(nc, variable_name))
#'
#' geom <- sf::read_sf(system.file("shape/nc.shp", package = "sf"))[20, ]
#' geom <- sf::st_transform(geom, "+init=epsg:5070")
#'
#' cell_geometry <-
#'   create_cell_geometry(X_coords = RNetCDF::var.get.nc(nc, cv$X, unpack = TRUE),
#'                        Y_coords = RNetCDF::var.get.nc(nc, cv$Y, unpack = TRUE),
#'                        prj = ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc)),
#'                        geom = geom,
#'                        buffer_dist = 2000)
#'
#' plot(sf::st_transform(sf::st_geometry(cell_geometry), "+init=epsg:5070"))
#' plot(sf::st_geometry(geom), add = TRUE)
#'
#' area_weights <- calculate_area_intersection_weights(
#'   sf::st_sf(dplyr::select(cell_geometry, grid_ids), agr = "constant"),
#'   sf::st_sf(dplyr::select(geom, NAME), agr = "constant"))
#'
#' execute_intersection(nc_file, variable_name, area_weights,
#'                      cell_geometry, cv$X, cv$Y, cv$T)
#'

execute_intersection <- function(nc_file,
                                 variable_name,
                                 intersection_weights,
                                 cell_geometry,
                                 x_var, y_var, t_var,
                                 start_datetime = NULL,
                                 end_datetime = NULL,
                                 status = FALSE,
                                 writer_fun = NULL,
                                 out_file = NULL,
                                 out_var_meta = NULL) {

  names(intersection_weights) <- c("grid_ids", "poly_id", "w")

  intersection_weights <- dplyr::filter(intersection_weights, !is.na(poly_id))

  nc <- open.nc(nc_file)
  on.exit(close.nc(nc), add = TRUE)

  ri <- get_request_inds(min(cell_geometry$X_ind), max(cell_geometry$X_ind),
                         min(cell_geometry$Y_ind), max(cell_geometry$Y_ind),
                         nc, variable_name, x_var, y_var,
                         t_var, start_datetime, end_datetime)

  rm(cell_geometry)

  ids <- as.integer(matrix(get_ids(length(ri$X), length(ri$Y)),
                           ncol = 1,
                           byrow = FALSE))

  out_nrows <- length(ri$T)

    out_ncols <- length(unique(intersection_weights[, 2][[1]]))

    out_names <- unique(intersection_weights$poly_id)

    if(!is.null(writer_fun)) {
      size <- c(out_nrows, out_ncols)
      var_meta <- c(out_var_meta, list(ids = out_names))
      file_handle <- writer_fun(out_file, step = 0, size = size, var_meta = var_meta)
    } else {
      out <- matrix(0, nrow = out_nrows, ncol = out_ncols)
    }

    transpose <- FALSE
    if(ri$dimid_order[1] > ri$dimid_order[2]) transpose <- TRUE

    intersection_weights <- data.table(intersection_weights)

    join_indices <- match(intersection_weights$grid_ids, ids)

    for (i in 1:out_nrows) {
      try_backoff({
        timer <- Sys.time()

        i_data <- get_i_data(i, nc, variable_name, ri, transpose,
                             intersection_weights, join_indices)

        i_data<-  as.numeric(i_data[, list(d = sum( (d * w), na.rm = TRUE)),
                                    by = poly_id]$d /
                               i_data[, list(sw = sum(w, na.rm = TRUE)),
                                      by = poly_id]$sw)

        if(!is.null(writer_fun)) {
          file_handle <- writer_fun(file_handle, step = i, size = size,
                                    var_meta = var_meta,
                                    timestep_data = list(ri$time_steps[i, ], i_data))

        } else {
          out[i, ] <- i_data
        }

        if(status) print(paste(i, "of", out_nrows, Sys.time() - timer))
      })
    }

    if(!is.null(writer_fun)) {
      writer_fun(file_handle, step = -1, size = size, var_meta = var_meta)
      return(out_file)
    } else {
      out <- data.frame(out)

      out <- select(out, which(!is.na(out_names)))

      out_names <- out_names[which(!is.na(out_names))]

      names(out) <- out_names

      time_steps <- data.frame(time_stamps = ri$time_steps)

      return(cbind(time_steps, out))
    }
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

get_request_inds <- function(min_x, max_x, min_y, max_y,
                             nc, variable_name, x_var, y_var,
                             t_var, start_datetime, end_datetime) {

  nc_var_info <- var.inq.nc(nc, variable_name)

  if (nc_var_info$ndims != 3) stop("only 3d variables are supported")

  X_var_info <- var.inq.nc(nc, x_var)
  Y_var_info <- var.inq.nc(nc, y_var)
  T_var_info <- var.inq.nc(nc, t_var)

  # (X_inds and Y_inds are the first and second dimensions of 2d coordinate variables)
  X_inds <- seq(min_x, max_x, 1)
  Y_inds <- seq(min_y, max_y, 1)

  if (X_var_info$ndims > 1 | Y_var_info$ndims > 1 | T_var_info$ndims > 1) {
    warning("2d coordinate variable found. Caution, limited testing.")
    warning("assuming dimension order")
    dimid_order <- match(nc_var_info$dimids,
                         c(X_var_info$dimids,
                           T_var_info$dimids))

    X_inds_dim <- X_var_info$dimids[1]
    Y_inds_dim <- X_var_info$dimids[2]
    T_inds_dim <- T_var_info$dimids
    # This whole thing needs to be refactored.
    if(X_inds_dim > Y_inds_dim &
       all(dimid_order[1:2] == c(1, 2)) &
       X_inds_dim < 2) {
      Y_inds <- seq(min_x, max_x, 1)
      X_inds <- seq(min_y, max_y, 1)
    }
  } else {
    dimid_order <- match(nc_var_info$dimids,
                         c(X_var_info$dimids,
                           Y_var_info$dimids,
                           T_var_info$dimids))
  }

  time_steps <- utcal.nc(att.get.nc(nc, T_var_info$name, "units"),
                         var.get.nc(nc, T_var_info$name, unpack = TRUE),
                         "c")

  time_steps <- data.frame(time_stamp = time_steps)

  if (is.null(start_datetime) & is.null(end_datetime)) {

    T_inds <- seq_len(nrow(time_steps))

  } else {

    T_inds <- get_time_inds(time_steps, start_datetime, end_datetime)
    time_steps <- time_steps[T_inds, ,drop = FALSE]

  }

  return(list(time_steps = time_steps, X = X_inds, Y = Y_inds, `T` = T_inds, dimid_order = dimid_order))

}

get_time_inds <- function(time_steps, start_datetime, end_datetime) {
  if (is.null(start_datetime)) start_datetime <- time_steps[1]

  if (is.character(start_datetime)) {
    start_datetime <- strptime(start_datetime,
                               format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }

  if (is.character(end_datetime)) {
    end_datetime <- strptime(end_datetime,
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  }

  T_inds <- time_steps$time_stamp >= start_datetime &
    time_steps$time_stamp <= end_datetime

  T_inds <- which(T_inds)

  if(length(T_inds) == 0) {
    stop("No time steps intersect start and stop datetimes.")
  }

  return(T_inds)
}

get_dap_url <- function(min_x, max_x, min_y, max_y,
                        nc_file, variable_name, x_var, y_var, t_var,
                        start_datetime, end_datetime) {

  nc <- open.nc(nc_file)
  on.exit(close.nc(nc), add = TRUE)

  ri <- get_request_inds(min_x, max_x, min_y, max_y,
                         nc, variable_name, x_var, y_var,
                         t_var, start_datetime, end_datetime)

  get_range <- function(min, max) sprintf("[%s:%s]", min, max)

  X_range <- get_range(min(ri$X), max(ri$X))
  Y_range <- get_range(min(ri$Y), max(ri$Y))
  T_range <- get_range(min(ri$T), max(ri$T))

  # OPeNDAP uses reverse order from RNetCDF so TYX here.
  var_inds <- paste0(c(T_range, Y_range, X_range)[ri$dimid_order], collapse = "")

  sprintf("%s?%s%s,%s%s,%s%s,%s%s",
          nc_file,
          x_var, X_range,
          y_var, Y_range,
          t_var, T_range,
          variable_name, var_inds)
}

get_i_data <- function(i, nc, variable_name, ri, transpose, intersection_weights, join_indices) {
  i_data <- var.get.nc(nc, variable_name,
                       start = c(min(ri$X), min(ri$Y), ri$T[i])[ri$dimid_order],
                       count = c(length(ri$X), length(ri$Y), 1)[ri$dimid_order],
                       unpack = TRUE)

  if(transpose) i_data <- t(i_data)

  i_data <- data.table(d = as.numeric(matrix(i_data,
                                             ncol = 1,
                                             byrow = FALSE)))

  i_data <- cbind(intersection_weights,
                  i_data[join_indices, ])

  i_data <- i_data[, w := ifelse(is.na(d), NA, w)]
}
