

###############################################################################
#########################   Assistant Functions   #############################
###############################################################################

find_dim_lon <- function(ncfile, var = NULL) {
  if(is.character(var)) {
    dim_names <- sapply(ncfile$var[[var]]$dim, function(x)x$name)
  } else {
    dim_names <- names(ncfile$dim)
  }

  for(dim_name in dim_names) {
    units <- ncfile$dim[[dim_name]]$units
    if(!is.null(units) &&
       (tolower(units) == "degrees_east" | tolower(units) == "degree_east"))
      return(dim_name)
    ldn <- tolower(dim_name)
    if(ldn == "long" | ldn == 'lng' | ldn == 'lon' | ldn == 'longitude')
      return(dim_name)
  }
}


find_dim_lat <- function(ncfile, var = NULL) {
  if(is.character(var)) {
    dim_names <- sapply(ncfile$var[[var]]$dim, function(x)x$name)
  } else {
    dim_names <- names(ncfile$dim)
  }

  for(dim_name in dim_names) {
    units <- ncfile$dim[[dim_name]]$units
    if(!is.null(units) &&
       (tolower(units) == "degrees_north" | tolower(units) == "degree_north"))
      return(dim_name)
    ldn <- tolower(dim_name)
    if(ldn == "lat" | ldn == 'la' | ldn == 'latitude')
      return(dim_name)
  }
}

find_dim_time <- function(ncfile, var = NULL) {
  if(is.character(var)) {
    dim_names <- sapply(ncfile$var[[var]]$dim, function(x)x$name)
  } else {
    dim_names <- names(ncfile$dim)
  }

  for(dim_name in dim_names) {
    units <- ncfile$dim[[dim_name]]$units
    if(!is.null(units) && grepl('since', tolower(units)) &&
       (grepl('day', tolower(units)) | grepl('month', tolower(units)) |
        grepl('year', tolower(units)) | grepl('minute', tolower(units)) |
        grepl('hour', tolower(units)) | grepl('second', tolower(units)) ))
      return(dim_name)
    ldn <- tolower(dim_name)
    if(ldn == "time" | ldn == 'date' | ldn == 't')
      return(dim_name)
  }
}
############################################################################


#' Quick build netCDF typical dimensions.
#'
#' @description Quick build typical netCDF dimensions like longitude, latitude and time.
#'
#' @param vals Values of the dimension. For lon and lat can be a series of longitudes and
#'             latitudes, for time can be a series of dates and times. If not provided,
#'             can provide from, csize and len for substitudes.
#' @param from Beginning of the values. Provided if vals was not provided.
#' @param csize Size of the gridcells for lon and lat. Also can be regarded as the interval
#'              of the lon and lat values.
#' @param len Length of the values.
#'
#' @param name Name of the dimensions. Default are lon, lat and time.
#' @param longname Long name of the dimensions.
#' @param time Time or date series for time dimension.
#' @param ... Other parameter for function time2num.
#'
#' @return A netCDF dimension object.
#'
#' @name ncdim
#' @export
ncdim_lon <- function(vals = NULL, from = NULL, csize = NULL, len = NULL, name = 'lon', longname = 'longitute') {
  if(!is.numeric(vals) & (!is.numeric(from) | !is.numeric(csize) | !is.numeric(len)))
    stop('vals and informations of grids must provide at lease one. Or format incorrect.')

  if(!is.numeric(vals)) {
    vals = seq(from = from, by = csize, length.out = len)
  }

  units <- "degrees_east"
  dim_lon <- ncdim_def(name, vals = vals, units = units, longname=longname)
  dim_lon
}

#' @rdname ncdim
#' @export
ncdim_lat <- function(vals = NULL, from = NULL, csize = NULL, len = NULL, name = 'lat', longname = 'latitute') {
  if(!is.numeric(vals) & (!is.numeric(from) | !is.numeric(csize) | !is.numeric(len)))
    stop('vals and informations of grids must provide at lease one. Or format incorrect.')

  if(!is.numeric(vals)) {
    vals = seq(from = from, by = csize, length.out = len)
  }

  units <- "degrees_north"
  dim_lon <- ncdim_def(name, vals = vals, units = units, longname=longname)
  dim_lon
}


#' @rdname ncdim
#' @export
ncdim_time <- function(time, name = 'time', ...) {
  ptime <- time2num(time, ...)
  dim_time <- ncdim_def(name, vals = ptime$val, units = ptime$units)
  dim_time
}

#' Convert numbers to time or convert time to numbers in nc format.
#'
#' @description Convert numbers of netCDF time dimension style to R time series by providing
#'              units or Convert R time series to numbers and units of netCDF time dimension
#'              style.
#'
#' @param time R time or date series.
#' @param since Time datum of the time series after converted to numbers. Default the first
#'              step of the time series.
#' @param freq Time step of the time series. Could select a value when not provided.
#' @param tz Time zone of the time series. Default UTC.
#'
#' @param units NetCDF time units for function num2time.
#'
#' @return A list including the numeric time series and the netCDF style time unit for function
#'         time2num, or an R time series for function num2time.
#'
#' @name nctime
#' @export
time2num <- function(time, since = time[1], freq = NULL, tz = 'UTC') {
  time <- as.POSIXct(time)
  since <- as.POSIXct(since)
  tnum <- as.numeric(time - since)

  if(is.null(freq)) {
    if(length(time) == 1) {
      freq <- 'day'
    } else {
      numintv <- mean(tnum[2:length(time)] - tnum[1:(length(time)-1)])
      if(numintv >= 31536000) {
        freq <- 'year'
      } else if(numintv >= 259200){
        freq <- 'month'
      } else if(numintv >= 86400){
        freq <- 'day'
        intv <- 86400
      } else if(numintv >= 3600){
        freq <- 'hour'
        intv <- 3600
      } else if(numintv >= 60){
        freq <- 'minute'
        intv <- 60
      } else {
        freq <- 'second'
        intv <- 1
      }
    }
  }

  if(freq %in% c('year', 'month')) {
    if(freq == 'month') {
      cumu.mon <- month(time) + (year(time) - 1) * 12
      since.mon <- month(since) + (year(since) - 1) * 12
      time.vals <- cumu.mon - since.mon
    } else {
      time.vals <- year(time) - year(since)
    }
    units <- paste(freq, 's since ', format(since, '%Y-%m-%d'), sep = '')

  } else {
    if(freq == 'day'){
      intv <- 86400
    } else if(freq == 'hour'){
      intv <- 3600
    } else if(freq == 'minute'){
      intv <- 60
    } else {
      intv <- 1
    }
    time.vals <- tnum / intv
    units <- paste(freq, 's since ', format(since, '%Y-%m-%d %H:%M:%S', tz = tz), sep = '')
  }
  return(list(vals = time.vals, units = units))
}

#' @rdname nctime
#' @export
num2time <- function(time, units) {
  units_esm <- strsplit(units, '\\ssince\\s')[[1]]
  freq <- units_esm[1]
  from <- units_esm[2]

  if(freq == 'years' | freq == 'months') {
    from_date <- strsplit(from, '[-\\/]|\\s')[[1]]
    stmonth <- as.numeric(from_date[2])
    styear <- as.numeric(from_date[1])

    if(freq == 'year') {
      out_years <- styear + time
      ts <- as.Date(paste(out_years, stmonth, '01', sep='-'))
    } else {
      orimonth <- stmonth + (styear-1)*12
      tmp_months <- orimonth + time
      out_years <- tmp_months %/% 12 + 1
      out_months <- (tmp_months-1) %% 12 + 1
      ts <- as.Date(paste(out_years, out_months, '01', sep='-'))
    }
    return(ts)
  }

  if(freq == 'minutes') {
    multi <- 60
  } else if(freq == 'hours') {
    multi <- 3600
  } else if(freq == 'days') {
    multi <- 86400
  }
  ts <- as.POSIXct(time * multi, origin = from, tz = "UTC")
  return(ts)
}

#' Get dimension values by one line code.
#'
#' @description Get typical dimension values, like lon, lat and time by one line code.
#'
#' @param file File path of the netCDF file.
#' @param var The name of the variable where to get the dimension.
#' @param conv Convent the time values from number to time or not. Default TRUE.
#'             Specified for function get_nctime.
#'
#' @return The values of the dimension of the netCDF file.
#'
#' @rdname get_ncdim
#' @export
get_nclon <- function(file, var=NULL) {
  ncfile <- nc_open(file)
  tryCatch({
    lon_name <- find_dim_lon(ncfile, var)
    lon_vals <- ncfile$dim[[lon_name]]$vals
  },
  finally = {nc_close(ncfile)})
  return(lon_vals)
}

#' @rdname get_ncdim
#' @export
get_nclat <- function(file, var=NULL) {
  ncfile <- nc_open(file)
  tryCatch({
    lat_name <- find_dim_lat(ncfile, var)
    lat_vals <- ncfile$dim[[lat_name]]$vals
  },
  finally = {nc_close(ncfile)})
  return(lat_vals)
}

#' @rdname get_ncdim
#' @export
get_nctime <- function(file, var = NULL, conv = TRUE) {
  ncfile <- nc_open(file)
  tryCatch({
    time_name <- find_dim_time(ncfile, var)
    time_vals <- ncfile$dim[[time_name]]$vals
    time_units <- ncfile$dim[[time_name]]$units
  },
  finally = {nc_close(ncfile)})

  if(conv)
    time_vals <- num2time(time_vals, time_units)
  return(time_vals)
}

#' Read ArcInfo ASCII grid data.
#'
#' @description Read the grid data from ArcInfo type ASCII file.
#'
#' @param grid_file A ArcInfo ASCII grid file.
#' @return A grid data in list type, including x, y, xllcorner,
#'         yllcorner and cellsize.
#'
#' @export
read_ascgrid <- function(file) {
  meta <- read.table(file, nrows = 6)
  meta[ , 1] <- tolower(meta[ , 1])
  xcor <- meta[3, 2]
  ycor <- meta[4, 2]
  csize <- meta[5, 2]
  ndval <- meta[6, 2]

  if(grepl('center', meta[3, 1]))
    xcor <- xcor - csize/2
  if(grepl('center', meta[4, 1]))
    ycor <- ycor - csize/2

  grid <- read.table(file, skip = 6)
  names(grid) <- NULL
  grid <- as.matrix(grid)
  grid[grid == ndval] <- NA
  grid <- grid[nrows:1, ]
  grid <- t(grid)

  x <- seq(xcor + csize/2, length.out = dim(grid)[1], by = csize)
  y <- seq(ycor + csize/2, length.out = dim(grid)[2], by = csize)

  return(list(x = x, y = y, xcor = xcor, ycor = ycor, csize = csize, grid = grid))
}

#' Read data from nc file by one line code.
#'
#' @description Using one line code to get data from netCDF file.
#'
#' @param file String. File path of the netCDF file.
#' @param var String. Name of the variable to be get data.
#' @param from A vector indicating where to start reading the values at each dimentions.
#'             If byval = TRUE, it's regarded as the dimention value; If byval = FALSE,
#'             it's regarded as the indexes.
#' @param to A vector indicating the end of the values to be reading at each dimentions.
#'           If byval = TRUE, it's regarded as the dimention value; If byval = FALSE,
#'           it's regarded as the indexes.
#' @param dims A vector indicating subset at which dimensions.
#'
#' @param lon A length 2 vector. Indicating the longitude range of the data to be get.
#' @param lat A length 2 vector. Indicating the latitude range of the data to be get.
#' @param time A length 2 vector. Indicating the time range of the data to be get.
#'
#' @param byval If set to TRUE, values of from and to would be regarded as dimension
#'              values; If False, it would be regarded as indexes. Default TRUE.
#' @param msize Maximum size (MB) of the data to be read. This is provide that too much data
#'              read into memory. If < 0 it will regarded as unlimited. Default 2048 MB (2GB).
#'
#' @return The values of the data of the netCDF file.
#' @export
read_ncvar <- function(file, var = NULL, from = NULL, to = NULL, dims = NULL,
                       lon = NULL, lat = NULL, time = NULL, byval = TRUE, msize = 2048) {
  nccon <- nc_open(file)
  val <- NULL
  tryCatch ({
    if(!is.character(var)) {
      var <- nccon$var[[1]]$name
      warning(sprintf("var not provided. Default getting the first variable: \"%s\".", var))
    }
    if(!(var %in% names(nccon$var)))
      stop(sprintf("Variable \"%s\" not exist in the nc file.", var))

    alldims <- sapply(nccon$var[[var]]$dim, function(x) x$name)
    ndims <- nccon$var[[var]]$ndims
    dimlens <- nccon$var[[var]]$size

    # Preset getting range.
    start <- rep(NA, ndims)
    end <- rep(NA, ndims)

    # Find out the site of the indexed dims of all dims.
    if(!is.character(dims)) {
      dims <- alldims
    } else if(any(!(dims %in% alldims))){
      stop(sprintf("Dimention \"%s\" is not exist in the nc file.",
                   dims[!(dims %in% alldims)]))
    }
    dimsites <- sapply(dims, function(x) which(x == alldims))


    if(!is.numeric(from))
      from <- start
    if(!is.numeric(to))
      to <- end
    if(ndims < length(from) | ndims < length(to))
      stop('length of from ,to or dim are incorrect.')

    start[dimsites[1:length(from)]] <- from
    end[dimsites[1:length(to)]] <- to

    if(any(start > end, na.rm = T))
      stop('from could not greater than to.')

    if(byval) {
      dimvals <- sapply(dims, function(x) nccon$dim[[x]]$val)
      start[is.na(start)] <- sapply(dimvals[is.na(start)], min)
      end[is.na(end)] <- sapply(dimvals[is.na(end)], max)

      inrange <- mapply(function(x, y, z)x[x >= y & x <= z],
                 dimvals, start, end)
      getlens <- sapply(inrange, length)
      if(any(getlens <= 0)) return(val)

      getstart <- mapply(function(x, y)which(x[1] == y), inrange, dimvals)
    } else {
      start[is.na(start) | start <= 0] <- 1
      start[start > dimlens] <- dimlens[start > dimlens]
      end[is.na(end)] <- dimlens[is.na(end)]
      end[end <= 0] <- 1
      end[end > dimlens] <- dimlens[end > dimlens]

      getstart <- start
      getlens <- dimlens - getstart + 1
    }

    ####################### Special for lon, lat and time ########################
    if(is.numeric(lon) && length(lon) >= 2) {
      nlon <- find_dim_lon(nccon)
      if(!(nlon %in% alldims))
        stop("Variable \"%s\" does not have longitude dimention.", var)
      slon <- which(nlon == alldims)
      lonval <- nccon$dim[[nlon]]$val
      if(is.na(lon[1])) lon[1] <- min(lonval)
      if(is.na(lon[2])) lon[2] <- max(lonval)
      if(lon[1] > lon[2])
        lon[1:2] <- lon[2:1]
      inrange <- lonval[lonval >= lon[1] & lonval <= lon[2]]
      getlens[slon] <- length(inrange)
      if(any(getlens[slon] <= 0)) return(val)
      getstart[slon] <- which(inrange[1] == lonval)
    }

    if(is.numeric(lat) && length(lat) >= 2) {
      nlat <- find_dim_lat(nccon)
      if(!(nlat %in% alldims))
        stop("Variable \"%s\" does not have latitude dimention.", var)
      slat <- which(nlat == alldims)
      lonval <- nccon$dim[[nlat]]$val
      if(is.na(lon[1])) lon[1] <- min(lonval)
      if(is.na(lon[2])) lon[2] <- max(lonval)
      if(lat[1] > lat[2])
        lat[1:2] <- lat[2:1]
      inrange <- latval[latval >= lat[1] & latval <= lat[2]]
      getlens[slat] <- length(inrange)
      if(any(getlens[slat] <= 0)) return(val)
      getstart[slat] <- which(inrange[1] == latval)
    }

    if(!is.null(time) && length(time) >= 2) {
      ntime <- find_dim_time(nccon)
      if(!(ntime %in% alldims))
        stop("Variable \"%s\" does not have time dimention.", var)
      stime <- which(ntime == alldims)
      timeunits <- nccon$dim[[ntime]]$units
      timeval <- num2time(nccon$dim[[ntime]]$val, timeunits)
      if(is.na(time[1])) {
        sttime <- as.POSIXct(time[1], tz='UTC')
      } else {
        sttime <- min(timeval)
      }
      if(is.na(time[2])) {
        edtime <- as.POSIXct(time[2], tz='UTC')
      } else {
        edtime <- max(timeval)
      }
      inrange <- timeval[timeval >= time[1] & timeval <= time[2]]
      getlens[stime] <- length(inrange)
      if(any(getlens[stime] <= 0)) return(val)
      getstart[stime] <- which(inrange[1] == timeval)
    }

    size <- prod(getlens) / 131072 # Get data size to be read and convent to MB.
    if(msize > 0 & size > msize)
      stop(sprintf('Size of value to be read is %.2f MB that larger than restrect size %.2f MB.
                   This is to avoid reading to much data out of memory.
                   Please set msize (in MB) to enlarge the restrected size or set to -1 for unlimited.',
                   size, msize))

    val <- ncvar_get(nccon, var, getstart, getlens)
  },
  finally = {
    nc_close(nccon)
  })
  return(val)
}

#write_ncvar <- function(val, file, varname, dims = NULL,) {
#  nccon <- nc_open(file, write = TRUE)

#  tryCatch ({


#  },
#  finally = {
#    nc_close(nccon)
#  })
#}

nc_combine <- function(file_list, out_file, combine_dim = 'time', meta_from = file_list[1],
                       from = NULL, to = NULL, compression = 1) {

  # Check params.
  if(is.numeric(from) & is.numeric(to) && from > to) stop('from must not larger than to.')

  # Getting meta data.
  ncf <- nc_open(meta_from)
  global_attr <- ncatt_get(ncf, 0)

  dim_names <- names(ncf$dim)
  dim_attrs <- list()
  for(dim_name in dim_names)
    dim_attrs[[dim_name]] <- ncatt_get(ncf, dim_name)

  names(dim_attrs) <- dim_names

  dim_vals <- list()
  for(dim_name in dim_names)
    dim_vals[[dim_name]] <- ncf$dim[[dim_name]]$vals

  # Getting Var informations.
  var_names <- names(ncf$var)

  var_dims <- list()
  var_attrs <- list()
  for(var_name in var_names) {
    var_dims[[var_name]] <-
      sapply(ncf$var[[var_name]]$dim, function(x) x$name)

    var_attrs[[var_name]] <- ncatt_get(ncf, var_name)
  }

  # Get where is the dimention that to be combine.
  cd_site <- sapply(var_dims, function(x) which(x == combine_dim))

  nc_close(ncf)

  # Check each nc file and get the total length.
  len_get <- c()
  from_get <- c()
  cdim_val <- c()
  from_put <- c()
  ifrom_put <- 1
  for(fn in file_list) {
    ncf <- nc_open(fn)
    # TODO Check the size of each var of each files.

    cd_len <- ncf$dim[[combine_dim]]$len
    ifrom <- from; ito <- to
    if(!is.numeric(ifrom) || ifrom < 1) ifrom <- 1
    if(!is.numeric(ito) || ito > cd_len) ito <- cd_len
    if(ifrom > cd_len) ifrom <- cd_len
    if(ito < 1) ito <- 1
    len_get <- append(len_get, ito - ifrom + 1)
    from_get <- append(from_get, ifrom)
    cdim_val <- append(cdim_val, ncf$dim[[combine_dim]]$vals[ifrom:ito])

    from_put <- append(from_put, ifrom_put)
    ifrom_put <- ifrom_put + ito - ifrom + 1

    nc_close(ncf)
  }

  dim_vals[[combine_dim]] <- cdim_val

  # Create assemble nc file.
  adims <- list()
  for(dim_name in dim_names) {
    adim <- ncdim_def(dim_name, 'None', dim_vals[[dim_name]])
    adims[[dim_name]] <- adim
  }

  adatas <- list()
  for(var_name in var_names) {
    adata <- ncvar_def(var_name, "None", adims[ var_dims[[var_name]] ],
                       compression = compression)
    adatas[[var_name]] <- adata
  }

  print(sprintf('Creating file %s', out_file))

  ncfa <- nc_create(out_file, adatas)

  # Set each attributes of meta data.
  for(attr_name in names(global_attr))
    ncatt_put(ncfa, 0, attr_name, global_attr[[attr_name]])

  for(var_name in var_names) {
    for(attr_name in names(var_attrs[[var_name]]))
      ncatt_put(ncfa, var_name, attr_name, var_attrs[[var_name]][[attr_name]])
  }

  for(dim_name in dim_names) {
    for(attr_name in names(dim_attrs[[dim_name]]))
      ncatt_put(ncfa, dim_name, attr_name, dim_attrs[[dim_name]][[attr_name]])
  }


  # Write in data from each original nc file.
  for(i in 1:length(file_list)) {
    fn <- file_list[i]
    ncf <- nc_open(fn)
    print(paste('Writing', fn))
    for(var_name in var_names) {

      pstart <- ostart <- rep(1, length(var_dims[[var_name]]))
      ostart[cd_site[[var_name]]] = from_get[i]
      pstart[cd_site[[var_name]]] = from_put[i]
      olen <- rep(-1, length(var_dims[[var_name]]))
      olen[cd_site[[var_name]]] = len_get[i]

      # Write in data.
      odata <- ncvar_get(ncf, var_name, ostart, olen)
      ncvar_put(ncfa, var_name, odata, pstart, olen)
    }
    nc_close(ncf)
  }

  nc_close(ncfa)
}
