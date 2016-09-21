# Base functions.
#
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Inner function
mk <- function (x) {
  x <- as.numeric(x)
  y <- x[!is.na(x)]
  n <- length(y)
  r <- sapply(1:n, function(i) sum(y[1:i] < y[i]))
  s <- sum(r)
  E <- n * (n + 1) / 4
  Var <- n * (n - 1) * (2*n + 5) / 72
  MK.value <- (s - E) / sqrt(Var)
  p.value <- abs(pnorm(MK.value)-pnorm(-MK.value))
  c(MK.value = MK.value, p.value = p.value)
}

endpoints.seasons <- function(x, on = "spring") {
  if (timeBased(x)) {
    NR <- length(x)
    x <- xts(NULL, order.by = x)
  }
  else NR <- NROW(x)
  if (!is.xts(x))
    x <- try.xts(x, error = "must be either xts-coercible or timeBased")

  posixltindex <- as.POSIXlt(.POSIXct(.index(x)), tz = indexTZ(x))$mon
  if (on == "winter") {
    tocal <- c(11, 0, 1)
  }
  else if (on == "spring") {
    tocal <- c(2, 3, 4)
  }
  else if (on == "summer") {
    tocal <- c(5, 6, 7)
  }
  else if (on == "autumn") {
    tocal <- c(8, 9, 10)
  }

  xi <- rep(0, NR)
  xi[posixltindex %in% tocal] <- 1
  if(xi[1] == 1) {
    ep <- as.integer(c(0, which(diff(xi) != 0)))
  }else {
    ep <- as.integer(which(diff(xi) != 0))
  }
  if(xi[NR] == 1) {
    ep[length(ep) + 1] <- NR
  }
  ep
}

season.apply <- function(x, INDEX, FUN, ...)
{
  x <- try.xts(x, error = FALSE)
  FUN <- match.fun(FUN)

  re <- sapply(1:(length(INDEX)/2), function(y) {
    FUN(x[(INDEX[y*2 - 1] + 1):(INDEX[y*2])], ...)
  })

  if (is.vector(re))
    re <- t(re)
  re <- t(re)
  if (is.null(colnames(re)) && NCOL(x) == NCOL(re))
    colnames(re) <- colnames(x)

  SINDEX <- INDEX[seq(2, length(INDEX), by = 2)]
  reclass(re, x[SINDEX])
}

####################################################

#' Mann-Kendall trend test
#'
#' @param x A time series (vector or matrix) to be test.
#' @return A list include M-K static value(MK.value) and significance value(p.value)
#'         if x is a vector, or a matrix to record MK.value and p.value for each column of x.
#'
#' @export
MK_test <- function (x) {
  if (is.null(dim(x)))
    return(as.list(mk(x)))
  if (ncol(x) == 1)
    return(as.list(mk(x[, 1])))
  return(t(apply(x, 2, mk)))
}


#' Apply function over seasons for xts object.
#'
#' @description The seasons is conclute as that March, April and May as spring,
#'              June, July and August as summer, September, October and November
#'              as autumn, December and January, Feburary of next year as winter.
#' @param x An xts time-series object
#' @param FUN An R function
#' @param ... Additional arguments to FUN.
#'
#' @return The result of the R function apply to the xts object.
#'
#' @name apply.season
#' @export
apply.seasonally <- function(x, FUN, ...) {
  if (timeBased(x)) {
    NR <- length(x)
    x <- xts(NULL, order.by = x)
  }
  else NR <- NROW(x)
  if (!is.xts(x))
    x <- try.xts(x, error = "must be either xts-coercible or timeBased")

  posixltindex <- as.POSIXlt(.POSIXct(.index(x)), tz = indexTZ(x))
  xi <- ((posixltindex$mon+2)%/%3)
  xi[xi == 4] <- 0
  ep <- as.integer(c(0, which(diff(xi) != 0), NR))

  period.apply(x, ep, FUN, ...)
}

#' @rdname apply.season
#' @export
apply.spring <- function(x, FUN, ...) {
  ep <- endpoints.seasons(x, "spring")
  season.apply(x, ep, FUN, ...)
}

#' @rdname apply.season
#' @export
apply.summer <- function(x, FUN, ...) {
  ep <- endpoints.seasons(x, "summer")
  season.apply(x, ep, FUN, ...)
}

#' @rdname apply.season
#' @export
apply.autumn <- function(x, FUN, ...) {
  ep <- endpoints.seasons(x, "autumn")
  season.apply(x, ep, FUN, ...)
}

#' @rdname apply.season
#' @export
apply.winter <- function(x, FUN, ...) {
  ep <- endpoints.seasons(x, "winter")
  season.apply(x, ep, FUN, ...)
}


#' Some template of ggplot2 theme may be offten used.
#'
#' @description Just as title.
#'
#' @return A ggplot2 theme.
#'
#' @name mytheme
#' @export
blank.theme <- function() {
  th <- theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              panel.border = element_rect(color = 'black', fill = NA, size = 0.5),
              strip.background = element_blank(), axis.line = element_blank(),
              axis.text = element_text(color = 'black'),
              text = element_text(color='black'))
  th
}

#' Quickly create a date sequence
#'
#' @description  Quickly create a date sequence without write a mountain of
#'               "as.Date" and write the long argument "length.out = ".
#' @param from A string of the initiate date, "1949", "1949-01", "1949-01-01",
#'             "1949/01", "1949/01/01", "1949.01", "1949.01.01" are acceptabel.
#' @param to A string of the end date, format is as same as from.
#' @param len The length of the sequence.
#'
#' @return A date sequence.
#'
#' @name seq.date
#' @export
seq_day <- function(from, to = NULL, len = NULL) {
  ncfrom <- length(strsplit(from, "[-./\\s]")[[1]])
  from <- sub("[./\\s]","-",from)
  if(ncfrom == 2) {
    from <- paste(from, "01", sep = "-")
  }else if(ncfrom == 1) {
    from <- paste(from, "01", "01", sep = "-")
  }else if(ncfrom == 0 | ncfrom > 3)
    stop('Error: fromat of "from" incorrect.')

  if(!is.null(to)) {
    ncto <- length(strsplit(to, "[-./\\s]")[[1]])
    to <- sub("[./\\s]","-",to)
    if(ncto == 2) {
      to <- paste(to, "01", sep = "-")
    }else if(ncto == 1) {
      to <- paste(to, "01", "01", sep = "-")
    }else if(ncto == 0 | to > 3)
      stop('Error: fromat of "to" incorrect.')

    sq <- seq(from = as.Date(from), to = as.Date(to), by = 'day')

  }else if(!is.null(len)){
    sq <- seq(from = as.Date(from), length.out = len, by = 'day')

  }else{
    stop("Error: to and len shouldn't be NULL at same time.")
  }
  sq
}

#' @rdname seq.date
#' @export
seq_month <- function(from, to = NULL, len = NULL) {
  ncfrom <- length(strsplit(from, "[-./\\s]")[[1]])
  from <- sub("[./\\s]","-",from)
  if(ncfrom == 2) {
    from <- paste(from, "12", sep = "-")
  }else if(ncfrom == 1) {
    from <- paste(from, "12", "01", sep = "-")
  }else if(ncfrom == 0 | ncfrom > 3)
    stop('Error: fromat of "from" incorrect.')

  if(!is.null(to)) {
    ncto <- length(strsplit(to, "[-./\\s]")[[1]])
    to <- sub("[./\\s]","-",to)
    if(ncto == 2) {
      to <- paste(to, "12", sep = "-")
    }else if(ncto == 1) {
      to <- paste(to, "12", "01", sep = "-")
    }else if(ncto == 0 | to > 3)
      stop('Error: fromat of "to" incorrect.')

    sq <- seq(from = as.Date(from), to = as.Date(to), by = 'month')

  }else if(!is.null(len)){
    sq <- seq(from = as.Date(from), length.out = len, by = 'month')

  }else{
    stop("Error: to and len shouldn't be NULL at same time.")
  }
  sq
}

#' @rdname seq.date
#' @export
seq_year <- function(from, to = NULL, len = NULL) {
  ncfrom <- length(strsplit(from, "[-./\\s]")[[1]])
  from <- sub("[./\\s]","-",from)
  if(ncfrom == 2) {
    from <- paste(from, "01", sep = "-")
  }else if(ncfrom == 1) {
    from <- paste(from, "01", "01", sep = "-")
  }else if(ncfrom == 0 | ncfrom > 3)
    stop('Error: fromat of "from" incorrect.')

  if(!is.null(to)) {
    ncto <- length(strsplit(to, "[-./\\s]")[[1]])
    to <- sub("[./\\s]","-",to)
    if(ncto == 2) {
      to <- paste(to, "01", sep = "-")
    }else if(ncto == 1) {
      to <- paste(to, "01", "01", sep = "-")
    }else if(ncto == 0 | to > 3)
      stop('Error: fromat of "to" incorrect.')

    sq <- seq(from = as.Date(from), to = as.Date(to), by = 'year')

  }else if(!is.null(len)){
    sq <- seq(from = as.Date(from), length.out = len, by = 'year')

  }else{
    stop("Error: to and len shouldn't be NULL at same time.")
  }
  sq
}
