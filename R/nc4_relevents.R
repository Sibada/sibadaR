


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

nc_time2num <- function(time, since = time[1], freq = NULL, tz = 'UTC') {
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
      time.val <- cumu.mon - since.mon
    } else {
      time.val <- year(time) - year(since)
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
    time.val <- tnum / intv
    units <- paste(freq, 's since ', format(since, '%Y-%m-%d %H:%M:%S', tz = tz), sep = '')
  }
  return(list(val = time.val, units = units))
}

nc_num2time <- function(time, units) {
  units_esm <- strsplit(units, '\\ssince\\s')[[1]]
  freq <- units_esm[1]
  from <- units_esm[2]

  if(freq == 'years' | freq == 'months') {
    from_date <- strsplit(from, '[-\\/ ]')[[1]]
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


find_dim_lon <- function(ncfile, var = NULL) {
  if(is.character(var)) {
    dim_names <- sapply(ncfile$var[[var]]$dim, function(x)x$name)
  } else {
    dim_names <- names(ncfile$dim)
  }

  for(dim_name in dim_names) {
    ldn <- tolower(dim_name)
    if(ldn == "long" | ldn == 'lng' | ldn == 'lon' | ldn == 'longitude')
      return(dim_name)
    units <- ncfile$dim[[dim_name]]$units
    if(!is.null(units) &&
       (tolower(units) == "degrees_east" | tolower(units) == "degree_east"))
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
    ldn <- tolower(dim_name)
    if(ldn == "lat" | ldn == 'la' | ldn == 'latitude')
      return(dim_name)
    units <- ncfile$dim[[dim_name]]$units
    if(!is.null(units) &&
       (tolower(units) == "degrees_north" | tolower(units) == "degree_north"))
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
    ldn <- tolower(dim_name)
    if(ldn == "time" | ldn == 'date' | ldn == 't')
      return(dim_name)
    units <- ncfile$dim[[dim_name]]$units
    if(!is.null(units) && grepl('since', tolower(units)) &&
       (grepl('day', tolower(units)) | grepl('month', tolower(units)) |
        grepl('year', tolower(units)) | grepl('minute', tolower(units)) |
        grepl('hour', tolower(units)) | grepl('second', tolower(units)) ))
      return(dim_name)
  }
}


nc_get_lon <- function(file, var=NULL) {
  ncfile <- nc_open(file)
  tryCatch({
    lon_name <- find_dim_lon(ncfile, var)
    lon_vals <- ncfile$dim[[lon_name]]$vals
  },
  finally = {nc_close(ncfile)})
  return(lon_vals)
}

nc_get_lat <- function(file, var=NULL) {
  ncfile <- nc_open(file)
  tryCatch({
    lat_name <- find_dim_lat(ncfile, var)
    lat_vals <- ncfile$dim[[lat_name]]$vals
  },
  finally = {nc_close(ncfile)})
  return(lat_vals)
}

nc_get_time <- function(file, var = NULL, conv = TRUE) {
  ncfile <- nc_open(file)
  tryCatch({
    time_name <- find_dim_time(ncfile, var)
    time_vals <- ncfile$dim[[time_name]]$vals
    time_units <- ncfile$dim[[time_name]]$units
  },
  finally = {nc_close(ncfile)})

  if(conv)
    time_vals <- nc_num2time(time_vals, time_units)
  return(time_vals)
}
