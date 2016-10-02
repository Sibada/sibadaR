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

mk_U <- function(x) {
  n <- length(x)
  r <- s <- 0
  for(i in 2:n){
    r[i] <- length((1:i)[x[1:i] < x[i]])
    s[i] <- s[i-1] + r[i]
  }
  k <- 1:n
  U <- (s - k*(k+1)/4) / sqrt(k*(k-1)*(2*k+5)/72)
  U[1] <- 0
  return(U)
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
#' @name sq.date
#' @export
sq.day <- function(from, to = NULL, len = NULL) {
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

#' @rdname sq.date
#' @export
sq.month <- function(from, to = NULL, len = NULL) {
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

#' @rdname sq.date
#' @export
sq.year <- function(from, to = NULL, len = NULL) {
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

#' Mann-Kendall mutation test
#' @description Make Mann-Kendall mutation test and draw plot. Where the UF and UB curve cross means a mutation may be happen.
#' @param x A single time series (can be xts object) to do that test.
#' If x has multi column, it will only choose the first column to do that test.
#' @param plot TRUE or FALSE. If it will draw a plot.
#' @param out.value TRUE or FALSE. If it return the UF and UB value series or plot.
#' @param index Index for x axis of plot.
#' @param p.size
#' @param l.size Point and line size of plot.
#' @return Plot or values of UF and UB series from M-K mutation test.
#' @export
mk_mut_test <- function(x, plot = TRUE, out.value = FALSE, index=NULL, p.size=3, l.size=0.6) {
  if(!is.null(dim(x)) && ncol(x) > 1){
    x <- x[ , 1]
    warning("x is not a single series, and now only use the first column.")
  }
  n <- length(x)
  isxts <- FALSE
  if(is.xts(x)){
    t.ind <- index(x)
    x <- as.numeric(x)
    isxts <- TRUE
  }

  if(!is.null(index)){
    ind <- index
  }else if(isxts){
    ind <- t.ind
  }else ind <- 1:n

  UF <- mk_U(x)
  UB <- mk_U(rev(x))
  UB <- -rev(UB)

  mut.result <- data.frame(UF, UB)
  if(isxts) mut.result <- xts(mut.result, t.ind)

  if(!plot & out.value) {
    return(mut.result)
  }

  conf.bound <- c(0)
  if(max(mut.result) > 1.5) conf.bound <- append(blef.bond, 1.645)
  if(min(mut.result) < -1.5) conf.bound <- append(blef.bond, -1.645)
  if(max(mut.result) > 1.85) conf.bound <- append(blef.bond, 1.96)
  if(min(mut.result) < -1.85) conf.bound <- append(blef.bond, -1.96)

    mut.melt <- melt(data.frame(mut.result))
    mut.melt$index <- ind

  mut.plot <- ggplot(data=mut.melt, aes(x=index, y=value, group=variable)) +
    geom_line(aes(linetype=variable), size=l.size) +
    geom_point(aes(shape=variable), size=p.size) +
    geom_hline(yintercept = conf.bound, linetype=3) +
    geom_hline(yintercept = 0) +
    ylab("M-K statistic") +
    theme(axis.title.x=element_blank(), legend.title=element_blank())

  if(plot & out.value){
    plot(mkm.plot)
    return(mut.result)
  }else  return(mut.plot)
}

#' Quickly get point shapes.
#' @description Quickly get what shapes of the point shapes look like in
#' ggplot2.
#' @param Size parameter of points.
#' @return A chart showing point shapes with its ID.
#' @exportB
show_pointshapes <- function(size=7) {
  xy <- merge(1:5, 1:6)[1:26, ]
  xy$s <- (xy$y-1)*5 + xy$x - 1
  p <- ggplot(data=xy,aes(x=x,y=y))+geom_text(aes(label=s),vjust=2.6, size=5) + scale_y_reverse(limits=c(6.5,1))
  p <- p +
    # To the begin I had ever try to use loop to make that, however I got a shit.
    # Maybe Hadley hate loop. Finally I use Excel to generate those below.
    geom_point(inherit.aes=F, aes(x=1, y=1), shape=0, size=size)+
    geom_point(inherit.aes=F, aes(x=2, y=1), shape=1, size=size)+
    geom_point(inherit.aes=F, aes(x=3, y=1), shape=2, size=size)+
    geom_point(inherit.aes=F, aes(x=4, y=1), shape=3, size=size)+
    geom_point(inherit.aes=F, aes(x=5, y=1), shape=4, size=size)+
    geom_point(inherit.aes=F, aes(x=1, y=2), shape=5, size=size)+
    geom_point(inherit.aes=F, aes(x=2, y=2), shape=6, size=size)+
    geom_point(inherit.aes=F, aes(x=3, y=2), shape=7, size=size)+
    geom_point(inherit.aes=F, aes(x=4, y=2), shape=8, size=size)+
    geom_point(inherit.aes=F, aes(x=5, y=2), shape=9, size=size)+
    geom_point(inherit.aes=F, aes(x=1, y=3), shape=10, size=size)+
    geom_point(inherit.aes=F, aes(x=2, y=3), shape=11, size=size)+
    geom_point(inherit.aes=F, aes(x=3, y=3), shape=12, size=size)+
    geom_point(inherit.aes=F, aes(x=4, y=3), shape=13, size=size)+
    geom_point(inherit.aes=F, aes(x=5, y=3), shape=14, size=size)+
    geom_point(inherit.aes=F, aes(x=1, y=4), shape=15, size=size)+
    geom_point(inherit.aes=F, aes(x=2, y=4), shape=16, size=size)+
    geom_point(inherit.aes=F, aes(x=3, y=4), shape=17, size=size)+
    geom_point(inherit.aes=F, aes(x=4, y=4), shape=18, size=size)+
    geom_point(inherit.aes=F, aes(x=5, y=4), shape=19, size=size)+
    geom_point(inherit.aes=F, aes(x=1, y=5), shape=20, size=size)+
    geom_point(inherit.aes=F, aes(x=2, y=5), shape=21, size=size)+
    geom_point(inherit.aes=F, aes(x=3, y=5), shape=22, size=size)+
    geom_point(inherit.aes=F, aes(x=4, y=5), shape=23, size=size)+
    geom_point(inherit.aes=F, aes(x=5, y=5), shape=24, size=size)+
    geom_point(inherit.aes=F, aes(x=1, y=6), shape=25, size=size)+
    geom_point(inherit.aes=F, aes(x=2, y=6), shape=26, size=size)
  p <- p + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
  return(p)
}

#' Quickly get line types.
#' @description Quickly get what line types look like in ggplot2.
#' @return A chart showing line types with its ID.
#' @param size Size parameter of lines.
#' @export
shwo_linetypes <- function(size=0.5) {
  g <- ggplot()
  for(i in 1:6) g <- g + geom_hline(yintercept = i, linetype = i,size = size)
  g <- g + scale_y_reverse(breaks = 6:1) +
    theme(axis.title = element_blank(), axis.text = element_text(size = 15),
          axis.ticks = element_blank())
  plot(g)
}
