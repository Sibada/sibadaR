

#########################   Assistant Functions   #############################

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
  return('Null')
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
  return('Null')
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
  return('Null')
}
################################### depart ################################




#' Show the variables or their dimensions of a netCDF file.
#'
#' @description Show the variables or their dimensions of a netCDF file.
#'
#' @param file Path of the nc file.
#' @param show.dim If is TRUE, also list out the dimensions of each variables.
#'                 Default FALSE.
#'
#' @return A vector of the variable names.
#'
#' @export
list_ncvar <- function(file, show.dim = FALSE) {
  nccon <- nc_open(file)
  tryCatch({
    vars <- names(nccon$var)
    if(show.dim) {
      varinfos <- c()
      maxlen <- max(nchar(vars))
      for(var in vars) {
        dims <- sapply(nccon$var[[var]]$dim, function(x) x$name)
        diminfo <- paste(dims, collapse = '\t')
        if(nchar(var) < maxlen - 5)
          var <- paste(var, '   ')
        varinfo <- paste(var, ':', diminfo, sep='\t')
        varinfos <- c(varinfos, varinfo)
        message(varinfo)
      }
    }
  }, finally = {
    nc_close(nccon)
  })
  vars
}



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
ncdim_time <- function(time, name = 'time', calendar = 'proleptic_gregorian', ...) {
  ptime <- time2num(time, ...)
  dim_time <- ncdim_def(name, vals = ptime$val, units = ptime$units, calendar = calendar)
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
  if(is.character(since) && grepl('since', since)) {
    units_esm <- strsplit(since, '\\ssince\\s')[[1]]
    since <- units_esm[2]
  }
  since <- as.POSIXct(since, tz = "UTC")

  time <- as.POSIXct(time, tz = 'UTC')
  tnum <- as.numeric(time) - as.numeric(since)

  if(is.null(freq)) {
    if(length(time) == 1) {
      freq <- 'day'
    } else {
      numintv <- mean(tnum[2:length(time)] - tnum[1:(length(time)-1)])
      if(numintv >= 31536000) {
        freq <- 'year'
      } else if(numintv >= 259200){
        freq <- 'month'
      } else if(numintv >= 86399){
        freq <- 'day'
        intv <- 86400
      } else if(numintv >= 3599){
        freq <- 'hour'
        intv <- 3600
      } else if(numintv >= 59){
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

    if(freq == 'years') {
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
#'
#' @param lon A length 2 vector. Indicating the longitude range of the data to be get.
#' @param lat A length 2 vector. Indicating the latitude range of the data to be get.
#'
#' @param msize Maximum size (MB) of the data to be read. This is provide that too much data
#'              read into memory. If < 0 it will regarded as unlimited. Default 2048 MB (2GB).
#'
#' @return A list, including the data of the netCDF file and the dimensions values.
#' @export
read_ncvar <- function(file, var = NULL, from = NULL, to = NULL,
                       lon = NULL, lat = NULL, msize = 2048) {
  # file = 'E:/meteodatas/CCI/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-monthly-fv02.2.nc'
  # var = NULL; from = NULL; to = NULL; dims = NULL; lon = NULL; lat = NULL; time = NULL
  # tmp = read_ncvar('E:/meteodatas/CCI/ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-monthly-fv02.2.nc',
  # lon=c(70,140), lat=c(15,60))
  nccon <- nc_open(file)
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
    if(is.null(from)) from <- rep(1, ndims)
    if(is.null(to)) to <- dimlens

    if(length(from) != ndims || length(to) != ndims) {
      message(sprintf("The dimensions of variable \"%s\" are:", var))
      print(data.frame(dim = alldims, length = dimlens))
      stop("Lengths of from or to are incorrect.")
    }

    from[is.na(from) | from == -1] <- 1
    from[is.na(to) | to == -1] <- dimlens[is.na(to) | to == -1]

    if(any(from > to, na.rm = T))
      stop('from could not greater than to.')

    from[from < 1] <- 1; from[from > dimlens] <- dimlens[from > dimlens]
    to[to < 1] <- 1; to[to > dimlens] <- dimlens[to > dimlens]

    ####################### Special for lon, lat and time ########################
    if(is.numeric(lon)) {
      nlon <- find_dim_lon(nccon)
      if(!(nlon %in% alldims))
        stop("Variable \"%s\" does not have longitude dimention.", var)
      slon <- which(nlon == alldims)
      lonval <- nccon$dim[[nlon]]$val
      if(is.na(lon[1])) lon[1] <- min(lonval)
      if(is.na(lon[2])) lon[2] <- max(lonval)
      if(lon[1] > lon[2])
        stop("Min longtitude could not larger than max longitude.")

      inrange <- which(lonval >= lon[1] & lonval <= lon[2])
      rangelen <- length(inrange)
      if(rangelen <= 0)
        stop("Not any grid included in the given range of longitude.")
      lonfrom <- inrange[1]
      lonto <- inrange[rangelen]
      from[slon] <- lonfrom
      to[slon] <- lonto
    }

    if(is.numeric(lat)) {
      nlat <- find_dim_lat(nccon)
      if(!(nlat %in% alldims))
        stop("Variable \"%s\" does not have latitude dimention.", var)
      slat <- which(nlat == alldims)
      latval <- nccon$dim[[nlat]]$val
      if(is.na(lat[1])) lat[1] <- min(latval)
      if(is.na(lat[2])) lat[2] <- max(latval)
      if(lat[1] > lat[2])
        stop("Min latitude could not larger than max latitude.")

      inrange <- which(latval >= lat[1] & latval <= lat[2])
      rangelen <- length(inrange)
      if(rangelen <= 0)
        stop("Not any grid included in the given range of latitude.")
      latfrom <- inrange[1]
      latto <- inrange[rangelen]
      from[slat] <- latfrom
      to[slat] <- latto
    }
    ############################# end of lon, lat, time ########################

    output <- list()
    for(i in 1:ndims) {
      dimvals <- nccon$dim[[alldims[i]]]$vals[from[i]:to[i]]
      attr(dimvals, 'units') <- nccon$dim[[alldims[i]]]$units
      if(!is.null(nccon$dim[[alldims[i]]]$calendar))
        attr(dimvals, 'calendar') <- nccon$dim[[alldims[i]]]$calendar
      output[[alldims[i]]] <- dimvals
    }

    getlens <- to - from + 1

    size <- prod(getlens) / 131072 # Get data size to be read and convent to MB.
    if(msize > 0 & size > msize)
      stop(sprintf('Size of value to be read is %.2f MB that larger than restrect size %.2f MB.
                   This is to avoid reading to much data out of memory.
                   Please set msize (in MB) to enlarge the restrected size or set to -1 for unlimited.',
                   size, msize))

    val <- ncvar_get(nccon, var, from, getlens)
    output$data <- val
  },
  finally = {
    nc_close(nccon)
  })
  return(output)
}


#' Write data to netCDF file.
#'
#' @description Write data to netCDF file.
#'
#' @param val Values to be write to the file. Should be matrix or array.
#' @param file Path of the netCDF file to be write.
#' @param var Name of the variable to be write. If the variable is not exist and
#'            the parameter "newvar" set to TRUE then the variable would be create.
#' @param dims A vector of string to determine the dimensions of the variable
#'             when create a new variable. If some dimensions was not determine,
#'             it would atuomatically find a dimension that fit it.
#' @param newvar If is TRUE, if the variable is not exist, it would create a
#'               new variable. Default FALSE.
#'
#' @param units Units of the new variable if the new variable is created.
#'
#' @export
write_ncvar <- function(val, file, var, dims = c(), newvar = FALSE, units = 'NULL') {
  if(!is.character(var))
    stop("Name of variable not provide or incorrect.")
  ncf <- nc_open(file, write = TRUE)
  tryCatch({
    vars <- names(ncf$var)
    if(!(var %in% vars) & !newvar)
      stop(sprintf('Variable "%s" is not exist in the nc file.', var))

    if(!(var %in% vars) & newvar) {
      vdims <- rep(NA, length(dim(val)))
      dimlens <- sapply(ncf$dim, function(x)x$len)
      dimnames <- names(dimlens)
      vdimlens <- dim(val)

      for(vdim in dims) {
        if(!(vdim %in% dimnames))
          stop(sprintf('Dimension "%s" is not exist in the nc file.', vdim))
        dimsite <- which(vdimlens == dimlens[vdim])
        if(length(dimsite) == 0)
          stop(sprintf('Dimension "%s" (length %d) can not match any dim of val.',
                       vdim, dimlens[vdim]))
        vdims[dimsite][which(is.na(vdims[dimsite]))[1]] <- vdim
      }
      for(i in 1:length(vdims)) {
        if(!is.na(vdims[i])) next
        dimlen <- dim(val) [i]
        dimsite <- which(dimlen == dimlens)
        if(length(dimsite) == 0)
          stop(sprintf("The No.%d of the dim of the val (length %d) does not match any dimension of the nc file.", i, dimlen))
        if(length(dimsite) > 1) {
          dimsite <- dimsite[1]
          warning(sprintf('The No.%d of the dim of the val (length %d) matchs more than 1 dimension.
                          Now select dimension "%s" to match it.', i, dimlen, dimnames[dimsite]))
        }
        vdims[i] <- dimnames[dimsite]
      }
      newdim <- list()
      for(vdim in vdims) {
        newdim[[vdim]] <- ncf$dim[[vdim]]
      }
      newvar <- ncvar_def(var, units, newdim)
      ncvar_add(ncf, newvar)

      # The variable could not be write just after be added.
      # Should close the file and reopen.
      nc_close(ncf)
      ncf <- nc_open(file, write = TRUE)
    }
    ncvar_put(ncf, var, val)

  }, finally = {
    nc_close(ncf)
  })
}


#' Combine several netCDF file by time.
#'
#' @description Combine several netCDF file by time. Each netCDF file should have the
#'              same dimensions and variables.
#'
#' @param fs A vector of string. File paths of each nc file to be combine.
#' @param out_file Output path of the new combined nc file.
#' @param from A vector indicating the start of the time range to be combine of each nc file.
#' @param to A vector indicating the end of the time range to be combine of each nc file.
#' @param datum_nc Which nc file to be regard as the datum, i.e. the attributes, dimensions
#'                 and variable names would inherit from it.
#' @param slip Provided if the time datum ("XXX" of the "days since XXX") of each file are not
#'             equal. It appoints the times interval (in time step) between each nc file to be
#'             combine. Default 30.
#'
#' @param compression The same parameters of function ncvar_def. Define the compression level
#'                    of the nc file.
#'
#' @export
nc_combine_time <- function(fs, out_file, from = NULL, to = NULL, datum_nc = 1, compression = NA, slip = 30) {

  if(length(fs) < 1)
    stop('Number of input files should more than one.')

  if(is.null(from))
    from <- rep(1, length(fs))
  if(is.null(to))
    to <- from * 0

  message('Collecting infos...')

  tlen <- c()
  tunits <- c()
  calendars <- c()
  tvals <- list()

  if(length(slip) < length(fs) - 1)
    slip <- rep(slip[1], length(fs))

  for(i in 1:length(fs)) {
    ncf <- nc_open(fs[i])
    tdname <- find_dim_time(ncf)
    timeinfo <- ncf$dim[[tdname]]
    tunits[i] <- timeinfo$units
    tvals[[i]] <- timeinfo$vals
    tlen[i] <- length(tvals[[i]])
    calendars[i] <- timeinfo$calendar
    nc_close(ncf)
    if(to[i] <= 0 | to[i] > tlen[i]) to[i] <- tlen[i]
    if(from[i] <= 0 | from[i] > tlen[i]) from[i] <- tlen[i]
  }
  allfrom <- cumsum(to - from + 1) - (to - from)

  # Define values of the new time dimension.
  timev <- tvals[[1]][from[1]:to[1]]
  for(i in 2:length(fs)) {
    if (tunits[i] == tunits[1]) {
      timev <- append(timev, tvals[[i]][from[i]:to[i]])
    } else {
      movev <- timev[length(timev)] - tvals[[i]][1] + slip[i]
      itimev <- tvals[[i]] + movev
      timev <- append(timev, itimev[from[i]:to[i]])
      message(sprintf('Time move by %f', movev))
    }
  }

  if(is.null(calendars[1])) calendars[1] <- 'proleptic_gregorian'
  dtime <- ncdim_def('time', tunits[1], timev, calendar = calendars[1])

  ncf <- nc_open(fs[datum_nc])
  tdname <- find_dim_time(ncf)
  newdims <- ncf$dim
  newdims[[tdname]] <- dtime

  vars <- names(ncf$var)
  varinfos <- list()
  attrs <- list()
  newvars <- list()

  # Write in the data and attributes of the old nc files.
  for(var in vars) {
    varinfo <- ncf$var[[var]]
    dims <-  sapply(varinfo$dim, function(x)x$name)
    size <- varinfo$size
    units <- varinfo$units
    attr <- ncatt_get(ncf, var)
    attr[['_FillValue']] <- NULL
    attrs[[var]] <- attr
    varinfos[[var]] <- list(dims = dims, size = size)

    vardims <- list()
    if(length(dims) != 0)
      for (i in 1:length(dims)) {
        vardims[[i]] <- newdims[[dims[i]]]
      }
    newvar <- ncvar_def(var, 'NONE',
                        vardims,
                        compression = compression)
    newvars[[var]] <- newvar
  }
  global_attrs <- ncatt_get(ncf, 0)

  nc_close(ncf)


  nncf <- nc_create(out_file, newvars, force_v4 = TRUE)
  for(var in vars)
    for(attr in names(attrs[[var]])){
      ncatt_put(nncf, var, attr, attrs[[var]][[attr]])
    }

  for(i in 1:length(fs)) {
    f <- fs[i]
    message(paste("Combing file ", f, '...', sep=''))
    ncf <- nc_open(f)

    tdname <- find_dim_time(ncf)
    for(var in vars) {
      message(paste("Combing variable ", var, '...', sep=''))
      if(length(ncf$var[[var]]$dim) <= 0)
        next

      ifrom <- rep(1, length(varinfos[[var]]$size))
      iwfrom <- ifrom
      ito <- varinfos[[var]]$size
      timesite <- which(tdname == varinfos[[var]]$dims)
      if(length(timesite) > 0) {
        ifrom[timesite] <- from[i]
        ito[timesite] <- to[i]
        iwfrom[timesite] <- allfrom[i]
      }
      ilen <- ito - ifrom + 1

      val <- ncvar_get(ncf, var, ifrom, ilen)
      ncvar_put(nncf, var, val, iwfrom, ilen)
      rm(val)
    }

    nc_close(ncf)
  }
  nc_close(nncf)
  message('Complete.')
}

#' Write ArcInfo ASCII grid data.
#'
#' @description Write the grid data to the file as a ArcInfo type ASCII file by
#'              offering a list of grid data like the output of points2grid()
#'              or create_flow_direction(), or a matrix and some necessary
#'              informations like cell size, x and y corner of the grid data and
#'              so on.
#'
#' @param grid A list of grid data or a matrix. When it's a matrix, xcor, ycor and
#'             csize are necessary.
#' @param out_file Output file path.
#' @param xcor X corner of the grid. Necessary when grid is a matrix.
#' @param ycor Y corner of the grid. Necessary when grid is a matrix.
#' @param csize Cell size of the grid. Necessary when grid is a matrix.
#' @param NA_value A value to reprent the NA value.
#'
#' @export
write_arc_grid <- function(grid, out_file, xcor=NULL, ycor=NULL, csize=NULL, NA_value=-9999) {
  if(is.list(grid)) {
    if(is.null(grid$xllcorner) | is.null(grid$yllcorner) | is.null(grid$cellsize)) {
      stop("grid should have xllcorner, yllconer and cellsize.")
    }
    ncols=nrow(grid$grid)
    nrows=ncol(grid$grid)
    xllcorner <- grid$xllcorner
    yllcorner <- grid$yllcorner
    cellsize <- grid$cellsize
    griddata <- grid$grid
  } else {
    if(!is.matrix(grid)) {
      stop("grid should be a list or a matirx.")
    }
    if(is.null(xcor) | is.null(ycor) | is.null(csize)) {
      stop("xcor, ycor and csize must provided when grid is a matrix.")
    }
    ncols=nrow(grid)
    nrows=ncol(grid)
    xllcorner <- xcor
    yllcorner <- ycor
    cellsize <- csize
    griddata <- grid
  }
  griddata <- t(griddata)
  griddata <- griddata[ncol(griddata):1, ]
  meta <- c(paste('ncols         ', ncols),
            paste('nrows         ', nrows),
            paste('xllcorner     ', round(xllcorner, 6)),
            paste('yllcorner     ', round(yllcorner, 6)),
            paste('cellsize      ', round(cellsize, 6)),
            paste('NODATA_value  ', NA_value))
  writeLines(meta, out_file)
  write.table(griddata, out_file, append=TRUE, row.names=FALSE, col.names=FALSE, na=paste(NA_value))
}


#' Read ArcInfo ASCII grid data.
#'
#' @description Read the grid data from ArcInfo type ASCII file.
#'
#' @param grid_file A ArcInfo ASCII grid file.
#' @return A grid data in list type, including nrows, ncols, xllcorner, yllcorner and cellsize.
#'
#' @export
read_arc_grid <- function(grid_file) {
  params <- strsplit(readLines(grid_file, 6), '\\s+')
  xcor <- as.numeric(params[[3]][2])
  ycor <- as.numeric(params[[4]][2])
  csize <- as.numeric(params[[5]][2])
  NA_value <- as.numeric(params[[6]][2])

  grid <- read.table(grid_file, skip=6)
  grid <- as.matrix(grid)
  ncols <- ncol(grid)
  nrows <- nrow(grid)
  grid <- t(grid)
  row.names(grid) <- NULL

  grid <- grid[ , nrows:1]
  grid[grid == NA_value] <- NA

  in_grid <- list(ncols=ncols,
                  nrows=nrows,
                  xllcorner=xcor,
                  yllcorner=ycor,
                  cellsize=csize,
                  grid=grid)

  return(in_grid)
}


#' Clip a part from a netCDF file to another netCDF file.
#'
#' @description Clip a part from a netCDF file to another netCDF file by providing dimensions, from an to.
#'
#' @param file A netCDF file to be clip.
#' @param out_file The netCDF file to store the cliped part of the original netCDF file.
#' @param from A vector of the border that where to begin to clip the variables at each
#'             dimensions. The length should be as same as dims. If byval is TRUE, or
#'             time should be the "date-like" integer in the format of YYYYMMDD, such as
#'             that 19600101 means "1960-01-01".
#' @param to   A vector of the border that where to stop to clip the variables at each
#'             dimensions. The length should be as same as dims, details see from.
#' @param dim Names of the dimensions to be clip. Default are longitude, latitude and time.
#' @param vars Variables to be clip. Default all the variables of the origin nc file.
#' @param byval Determine the from and to are values of the dimensions (e.g. longitude
#'              values) or the subscript of the dimensions (e.g. the N th value of the
#'              variable). Default TRUE, use the values of the dimensions.
#' @param compression The same agrument of function ncvar_def.
#'
#' @export
nc_intercept <- function(file, out_file, from = NULL, to = NULL, dims = NULL,
                         vars = NULL, byval = TRUE, compression = NA) {
  nc_in <- nc_open(file)
  tryCatch({
    if(is.null(dims)) {
      #### Default longitude, latitude and time.
      dims <- c(find_dim_lon(nc_in),
                find_dim_lat(nc_in),
                find_dim_time(nc_in))
    }
    if(length(from) < length(dims))
      from <- c(from, rep(NA, length(dims) - length(from)))
    if(length(from) > length(dims))
      from <- from[1:length(dims)]
    if(length(to) < length(dims))
      to <- c(to, rep(NA, length(dims) - length(to)))
    if(length(to) > length(dims))
      to <- to[1:length(dims)]

    if(any(from > to, na.rm=T))
      stop('From could not larger than to.')

    dimlens <- c()
    for(idim in dims) {
      ilen <- nc_in$dim[[idim]]$len
      dimlens <- c(dimlens, ilen)
    }
    dimunits <- c()
    for(idim in dims) {
      iu <- nc_in$dim[[idim]]$units
      dimunits <- c(dimunits, iu)
    }
    time_site <- which(dims == 'time' | grepl('since', dimunits))

    dimvals <- list()
    for(idim in dims) {
      ival <- nc_in$dim[[idim]]$vals
      dimvals[[idim]] <- ival
    }

    if(byval) {
      minvs <- sapply(dimvals, min)
      maxvs <- sapply(dimvals, max)

      if(length(time_site) > 0) {
        for(time_s in time_site) {
          timevals <- dimvals[[time_s]]
          timests <- num2time(timevals, dimunits[[time_s]])
          if(!is.na(from[time_s])){
            fromtime <- paste(from[time_s])
            sttime <- as.POSIXct(fromtime, format = "%Y%m%d", tz = 'UTC')
            from[time_s] <- min(timevals[sttime <= timests])
          }
          if(!is.na(to[time_s])){
            totime <- paste(to[time_s])
            sttime <- as.POSIXct(totime, format = "%Y%m%d", tz = 'UTC')
            to[time_s] <- max(timevals[sttime >= timests])
          }
        }
      }
      from[is.na(from)] <- minvs[is.na(from)]
      to[is.na(to)] <- maxvs[is.na(to)]
      from[from < minvs] <- minvs[from < minvs]
      to[to > maxvs] <- maxvs[to > maxvs]

      sns <- list()
      for(d in 1:length(dims)) {
        dimval <- dimvals[[d]]
        sn <- which(dimval <= to[d] & dimval >= from[d])
        sns[[dims[d]]] <- sn
      }
    } else {
      from[is.na(from)] <- 1
      to[is.na(to)] <- dimlens[is.na(to)]
      from[from < 1] <- 1
      to[to > dimlens] <- dimlens[to > dimlens]

      sns <- list()
      for(d in 1:length(dims)) {
        sns[[dims[d]]] <- floor(from[d]) : ceiling(to[d])
      }
    }

    ########### Build new dims and vars ###########
    newdims <- list()
    for(idim in names(nc_in$dim)) {
      if(idim %in% dims) {
        vals <- nc_in$dim[[idim]]$vals[sns[[idim]]]
      } else {
        vals <- nc_in$dim[[idim]]$vals
      }
      calendar <- NA
      if(!is.null(nc_in$dim[[idim]]$calendar))
        calendar <- nc_in$dim[[idim]]$calendar
      newdim <- ncdim_def(idim, nc_in$dim[[idim]]$units, vals, calendar = calendar)
      newdims[[idim]] <- newdim
    }

    if(is.null(vars))
      vars <- names(nc_in$var)
    delvars <- c()

    newvars <- list()
    for(ivar in vars) {
      idimns <- sapply(nc_in$var[[ivar]]$dim, function(x) x$name)
      if(length(idimns) == 0){
        delvars <- c(ivar)
        next
      }
      idims <- list()
      for(idimn in idimns) {
        idims[[idimn]] <- newdims[[idimn]]
      }
      newvar <- ncvar_def(ivar, 'TEMP', idims, compression = compression)
      newvars[[ivar]] <- newvar
    }
    vars <- vars[!(vars %in% delvars)]

    nc_out <- nc_create(out_file, newvars)

    ########### Copy attributes. ##############
    for(ivar in vars) {
      attrs <- ncatt_get(nc_in, ivar)
      attrs$`_FillValue` <- NULL
      for(attr in names(attrs)) {
        ncatt_put(nc_out, ivar, attr, attrs[[attr]])
      }
    }
    attrs <- ncatt_get(nc_in, 0)
    for(attr in names(attrs)) {
      ncatt_put(nc_out, 0, attr, attrs[[attr]])
    }

    ########### Write in variables. ##########
    for(ivar in vars) {
      message(sprintf("Writing variable %s ...", ivar))
      idimns <- sapply(nc_in$var[[ivar]]$dim, function(x) x$name)
      ifrom <- rep(1, length(idimns))
      ilen <- rep(-1, length(idimns))
      for(i in 1:length(idimns)) {
        if(idimns[i] %in% dims) {
          ilen[i] = length(sns[[idimns[i]]])
          ifrom[i] = sns[[idimns[i]]][1]
        }
      }
      val <- ncvar_get(nc_in, ivar, ifrom, ilen)
      ncvar_put(nc_out, ivar, val)
      rm(val)
    }
    nc_close(nc_out)
  }, finally = {
    nc_close(nc_in)
  })
}

#' Make table-like array to 2d geo-like array
#'
#' @description Make the table-like data to grid-like data. For example,
#'              an arrray has three dimensions, station, height and time.
#'              Then it transfrom the data into a four dimensions
#'              data, with longitude, latitude, height and time as
#'              dimensions. If the original data is a 2d data (grid and
#'              time) then it would turn it to a 3d array. the coordinates
#'              of the grids must be grid-like.
#' @param arr The array of the data to be transform.
#' @param x,y Coordinates of the grids.
#'
#' @param csize Cell size of the grids.
#' @param y Coordinates of the y dimension.
#' @param dimpnt Site of the dimension difine the location of grids.
#'
#' @return An list, including the array tranformed, and the corresponding
#'         coordinates.
#'
#' @export
table2geo <- function(arr, x, y, csize = NULL, dimpnt = NULL, dgt = 5) {
  dims <- dim(arr)
  ldim <- length(dims)
  if(is.null(dimpnt)) {
    dimpnt <- which(dims == length(x))[1]
    if(is.na(dimpnt)) stop('All dimensions of arr not fit to x and y.')
  }
  if(length(dimpnt > 1)) dimpnt <- dimpnt[1]
  if(dims[dimpnt] != length(x) | length(x) != length(y))
    stop('Dimension of points not fit to x and y.')
  if(is.null(csize)) {
    csize <- min(c(diff(sort(unique(round(x, dgt)))),
                   diff(sort(unique(round(y, dgt))))))
  }

  row <- x/csize - min(x/csize) + 1
  col <- y/csize - min(y/csize) + 1
  nrow <- max(row)
  ncol <- max(col)

  newdim <- dims
  newdim <- c(dims[-dimpnt], nrow, ncol)
  nldim <- length(newdim)

  arr <- aperm(arr, c((1:ldim)[-dimpnt], dimpnt))
  newarr <- array(dim = newdim)
  blocksize <- prod(dims[-dimpnt])

  for(i in 1:dims[dimpnt]) {
    ir <- row[i]; ic <- col[i]
    spsite <- (ic-1) * nrow + ir
    newarr[(1:blocksize) + (spsite - 1)*blocksize] <-
      arr[(1:blocksize) + (i-1)*blocksize]
  }

  newarr <- aperm(newarr, c(nldim-1, nldim, (1:nldim)[1:(nldim-2)]))
  list(arr = newarr, x = sort(unique(x)), y = sort(unique(y)))

}


#' Make 2d grid-like data to station-like table data
#'
#' @description Make the grid-like data to table-like data. For example,
#'              an nc data has four dimensions, longitude, latitude, height
#'              and time. Then it transfrom the data into a three dimensions
#'              data, with grid, height and time as dimensions.
#'              If the original data is a three dimension data (x, y and
#'              time) then it would turn it to a table, each row is a time
#'              step while each column is a grid.
#' @param arr The array of the data to be transform.
#' @param dimxy Which dimensions are the x and y. Default c(1, 2), which
#'              means the first and second dimensions are x and y.
#'
#' @param x Coordinates of the x dimension.
#' @param y Coordinates of the y dimension.
#' @param mask An matrix with the same dim of x and y, marking which the
#'             grids to be transform (NA means not be transform).
#'
#' @return An list, including the array tranformed, and the corresponding
#'         coordinates.
#'
#' @export
geo2table <- function(arr, dimxy = NULL, x = NULL, y = NULL, mask = NULL) {
  dims <- dim(arr)

  if(is.null(x) | is.null(y)) {
    if(is.null(dimxy)) {
      dimxy <- c(1, 2)
    } else {
      x <- 1:dims[dimxy[1]]
      y <- 1:dims[dimxy[2]]
    }
  } else {
    if (length(which(dims == length(x))) < 1 |
        length(which(dims == length(y))) < 1) {
      stop("Length of x and y are not fit to arr.")
    }
    dimxy <- 1:2
    dimxy[1] <- which(dims == length(x))[1]
    dimxy[2] <- which(dims == length(y))[1]
  }

  lenxy <- dims[dimxy]

  if(is.null(mask)) {
    mask <- apply(arr, dimxy, function(ia) {
      ifelse(any(!is.na(ia)), 1, NA)
    })
  } else {
    if(any(dim(mask)[1:2] != lenxy))
      stop("Dim of mask is not corresponded to arr.")
  }

  mask[!is.na(mask)] <- 1
  ng <- sum(mask, na.rm = T)

  exitdim <- dims[-dimxy]
  flat.arr <- array(dim = c(exitdim, ng))
  blocksize <- prod(exitdim)

  arr <- aperm(arr, c((1:length(dims))[-dimxy], dimxy))

  lon <- c()
  lat <- c()
  sn <- 1
  for(i in 1:lenxy[1]) {
    for(j in 1:lenxy[2]) {
      if(!is.na(mask[i, j])) {
        spsite <- (j-1) * lenxy[1] + i
        flat.arr[1:blocksize + (sn-1)*blocksize] <-
          arr[(1:blocksize) + (spsite - 1)*blocksize]
        lon[sn] <- x[i]
        lat[sn] <- y[j]
        sn <- sn + 1
      }
    }
  }
  return(list(arr = flat.arr, x = lon, y = lat))
}
