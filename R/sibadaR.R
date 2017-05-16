# Base functions.
#
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

# Inner function

mk <- function (x, sen.slope = TRUE, p.value = TRUE) {
  x <- as.numeric(x)
  t <- 1:length(x)
  t <- t[!is.na(x)]
  x <- x[!is.na(x)]
  n <- length(x)
  S <- sum(sapply(1:n, function(k) sum(sign(x[k:n] - x[k])) ) )
  ties <- rle(sort(x))$lengths
  tfix <- sum(ties * (ties - 1) * (2*ties + 5) )
  varS <- (n * (n - 1) * (2*n + 5) - tfix) / 18
  Z <- (S - sign(S)) / sqrt(varS)
  mk.results <- c(Z = Z)
  if(p.value)
    mk.results <- c(mk.results, p = 1 - 2 * pnorm(-abs(Z)))

  if(sen.slope){
    sen <- c()
    #senr <- c()
    for(k in 2:n) {
      kslo <- (x[k] - x[1:(k-1)]) / (t[k] - t[1:(k-1)])
      sen <- append(sen, kslo)
      #senr <- append(senr, kslo/x[1:(k-1)])
    }
    mk.results <- c(mk.results, sen.slope = median(sen))
  }
  mk.results
}

mk_U <- function(x) {
  n <- length(x)
  r <- sapply(1:n, function(k)sum(x[1:k] < x[k]))
  s <- cumsum(r)
  k <- 1:n
  var <- sapply(k, function(ki) {
    tf <- rle(sort(x[1:ki]))$lengths
    (ki * (ki-1) * (2*ki+5) - sum(tf * (tf-1) * (2*tf+5)))/72
  })
  E <- sapply(k, function(ki) {
    tf <- rle(sort(x[1:ki]))$lengths
    (ki*(ki - 1) - sum(tf * (tf-1)))/4
  })
  U <- (s - E) / sqrt(var)
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
#' @param sen.slope Output Sen's slope or not. Default TRUE.
#' @param p.value Output significance confidence value or not. Default TRUE.
#' @return A named vector include M-K static value(Z), Sen's slope value
#'         (sen.slope) and significance confidence value (p) when x is
#'         a vector, or a matrix to record Z, sen.slope and p for
#'         each column of x.
#'
#' @export
MK_test <- function (x, sen.slope = TRUE, p.value = TRUE) {
  mk
  if (is.null(dim(x)))
    return(mk(x))
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

# Assistant function for sq.date
todate <- function(ori.date, is.end = FALSE) {
  if(is.numeric(ori.date)) ori.date <- as.character(ori.date)

  date.cpns <- strsplit(ori.date, "[-./\\s\\\\]")[[1]] # Split can be "-", ".", "/", "\" and space.
  num.cpns <- length(date.cpns)

  if(num.cpns == 1) {
    strlen <- nchar(date.cpns)
    if(strlen >= 7) {
      date.cpns <- c(substr(date.cpns, 1, strlen-4),
                     substr(date.cpns, strlen-3, strlen-2),
                     substr(date.cpns, strlen-1, strlen))
      num.cpns <- 3
    } else if(strlen >=5) {
      date.cpns <- c(substr(date.cpns, 1, strlen-2),
                     substr(date.cpns, strlen-1, strlen))
      num.cpns <- 2
    } else {
      if(!is.end)
        return(paste(date.cpns[1], "01", "01", sep = '-'))
      return(paste(date.cpns[1], "12", "31", sep = '-'))
    }
  }

  if(num.cpns >= 3)
    return(paste(date.cpns[1], date.cpns[2], date.cpns[3], sep = '-'))
  if(num.cpns == 2) {
    date.start <- as.Date(paste(date.cpns[1], date.cpns[2], '01', sep = '-'))
    if(!is.end)
      return(date.start)
    return(paste(date.cpns[1], date.cpns[2], days_in_month(date.start), sep = '-'))
  }
  if(num.cpns == 0)
    stop(sprintf('Error: fromat of date "%s" is incorrect.', ori.date))
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
  from <- todate(from)
  if(!is.null(to)) {
    to <- todate(to, is.end = TRUE)
    sq <- seq(from = as.Date(from), to = as.Date(to), by = 'day')
  }else if(!is.null(len)){
    sq <- seq(from = as.Date(from), length.out = len, by = 'day')
  }else{
    stop("Error: to and len should provide at lease one.")
  }
  sq
}

#' @rdname sq.date
#' @export
sq.month <- function(from, to = NULL, len = NULL) {
  from <- todate(from)
  if(!is.null(to)) {
    to <- todate(to, is.end = TRUE)
    sq <- seq(from = as.Date(from), to = as.Date(to), by = 'month')
  }else if(!is.null(len)){
    sq <- seq(from = as.Date(from), length.out = len, by = 'month')
  }else{
    stop("Error: to and len should provide at lease one.")
  }
  sq
}

#' @rdname sq.date
#' @export
sq.year <- function(from, to = NULL, len = NULL) {
  from <- todate(from)
  if(!is.null(to)) {
    to <- todate(to, is.end = TRUE)
    sq <- seq(from = as.Date(from), to = as.Date(to), by = 'year')
  }else if(!is.null(len)){
    sq <- seq(from = as.Date(from), length.out = len, by = 'year')
  }else{
    stop("Error: to and len should provide at lease one.")
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
#' @param p.size Point  size of plot.
#' @param l.size Line size of plot.
#' @return Plot or values of UF and UB series from M-K mutation test.
#' @export
MK_mut_test <- function(x, plot = TRUE, out.value = FALSE, index = NULL, p.size = 3, l.size = 0.6) {
  if(!is.null(dim(x))){
    if(ncol(x) > 1)
      warning("x is not a single series, and now only use the first column.")
    x <- x[ , 1]
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
  if(max(mut.result) > 1.5) conf.bound <- append(conf.bound, 1.645)
  if(min(mut.result) < -1.5) conf.bound <- append(conf.bound, -1.645)
  if(max(mut.result) > 1.85) conf.bound <- append(conf.bound, 1.96)
  if(min(mut.result) < -1.85) conf.bound <- append(conf.bound, -1.96)

  mut.melt <- melt(data.frame(index = ind, mut.result), id.vars = 'index')

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

#' Quickly get point shapes of ggplot2.
#' @description Quickly get what shapes of the point shapes look like in
#' ggplot2.
#' @param Size parameter of points.
#' @return A chart showing point shapes with its ID.
#' @export
gg_pointshapes <- function(size=7) {
  xy <- merge(1:5, 1:6)[1:26, ]
  xy$s <- (xy$y-1)*5 + xy$x - 1
  p <- ggplot(data=xy,aes(x=x,y=y)) +
    geom_text(aes(label=s),vjust=2.6, size=5) +
    scale_y_reverse(limits=c(6.5,1))
  p <- p +
    # At first I had ever try to use loop to make that but I got a shit in the end.
    # Maybe Hadley hate loop. Finally I use Excel to generate those below.
    geom_point(inherit.aes=F, aes(x=1, y=1), shape=0, size=size) +
    geom_point(inherit.aes=F, aes(x=2, y=1), shape=1, size=size) +
    geom_point(inherit.aes=F, aes(x=3, y=1), shape=2, size=size) +
    geom_point(inherit.aes=F, aes(x=4, y=1), shape=3, size=size) +
    geom_point(inherit.aes=F, aes(x=5, y=1), shape=4, size=size) +
    geom_point(inherit.aes=F, aes(x=1, y=2), shape=5, size=size) +
    geom_point(inherit.aes=F, aes(x=2, y=2), shape=6, size=size) +
    geom_point(inherit.aes=F, aes(x=3, y=2), shape=7, size=size) +
    geom_point(inherit.aes=F, aes(x=4, y=2), shape=8, size=size) +
    geom_point(inherit.aes=F, aes(x=5, y=2), shape=9, size=size) +
    geom_point(inherit.aes=F, aes(x=1, y=3), shape=10, size=size) +
    geom_point(inherit.aes=F, aes(x=2, y=3), shape=11, size=size) +
    geom_point(inherit.aes=F, aes(x=3, y=3), shape=12, size=size) +
    geom_point(inherit.aes=F, aes(x=4, y=3), shape=13, size=size) +
    geom_point(inherit.aes=F, aes(x=5, y=3), shape=14, size=size) +
    geom_point(inherit.aes=F, aes(x=1, y=4), shape=15, size=size) +
    geom_point(inherit.aes=F, aes(x=2, y=4), shape=16, size=size) +
    geom_point(inherit.aes=F, aes(x=3, y=4), shape=17, size=size) +
    geom_point(inherit.aes=F, aes(x=4, y=4), shape=18, size=size) +
    geom_point(inherit.aes=F, aes(x=5, y=4), shape=19, size=size) +
    geom_point(inherit.aes=F, aes(x=1, y=5), shape=20, size=size) +
    geom_point(inherit.aes=F, aes(x=2, y=5), shape=21, size=size) +
    geom_point(inherit.aes=F, aes(x=3, y=5), shape=22, size=size) +
    geom_point(inherit.aes=F, aes(x=4, y=5), shape=23, size=size) +
    geom_point(inherit.aes=F, aes(x=5, y=5), shape=24, size=size) +
    geom_point(inherit.aes=F, aes(x=1, y=6), shape=25, size=size) +
    geom_point(inherit.aes=F, aes(x=2, y=6), shape=26, size=size)
  p <- p + theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank())
  return(p)
}

#' Quickly get line types of ggplot2.
#' @description Quickly get what line types look like in ggplot2.
#' @return A chart showing line types with its ID.
#' @param size Size parameter of lines.
#' @export
gg_linetypes <- function(size=0.5) {
  g <- ggplot()
  for(i in 1:6) g <- g + geom_hline(yintercept = i, linetype = i,size = size)
  g <- g + scale_y_reverse(breaks = 6:1) +
    theme(axis.title = element_blank(), axis.text = element_text(size = 15),
          axis.ticks = element_blank())
  plot(g)
}


#' Quickly get point shapes of base plot.
#' @description Quickly get what shapes of the point shapes look like in
#' base plot function.
#' @return A chart showing point shapes with its ID.
#' @export
p_pointshapes <- function() {
  xy <- merge(1:5, 1:5)
  s <- 1:25
  ps <- data.frame(xy,s)
  names(ps)  <-  c("x", "y", "s")
  plot(x = c(), xlim = c(0.5,5.5), ylim = c(0.5,5.5), xlab = "type: pch, size: cex", ylab = "", asp = 1)
  for(i in 1:25) {
    points(ps$x[i], ps$y[i], pch = ps$s[i], cex = 2)
    text(ps$x[i] + 0.1, ps$y[i] - 0.3, labels = ps$s[i])
  }
}


#' Quickly get line types of base plot.
#' @description Quickly get what the line types look like in
#' base plot function.
#' @return A chart showing line types with its ID.
#' @export
p_linetypes <- function() {
  plot(x = c(), xlim =c(-1, 1), ylim = c(0.5,6.5), xlab = "type: lty, size: lwd", ylab = "")
  for(lt in 1:6) {
    lines(c(-1, 1), c(lt, lt), lty = lt, lwd = 2)
  }
}

#' Root mean square error
#' @description Calculate root mean square error of two data series.
#' @param sim Simulated data series or data to be evaluated.
#' @param obs Observed data series or data as reference.
#' @return Root mean square error of sim and obs series.
#' @export
RMSE <- function(sim, obs) {
  valid <- !is.na(sim) & !is.na(obs)
  sim <- sim[valid]
  obs <- obs[valid]
  return(sqrt(mean((sim - obs)**2)))
}

#' Normailzed mean square error
#' @description Calculate normailzed mean square error of two data series.
#' @param sim Simulated data series or data to be evaluated.
#' @param obs Observed data series or data as reference.
#' @return Normailzed mean square error of sim and obs series.
#' @export
NMSE <- function(sim, obs) {
  valid <- !is.na(sim) & !is.na(obs)
  sim <- sim[valid]
  obs <- obs[valid]
  return(mean((obs - sim)**2) / mean((obs - mean(obs))**2) )
}

#' Nash-Sutcliffe coefficient of efficiency
#' @description Calculate Nash-Sutcliffe coefficient of efficiency of two data series.
#' @param sim Simulated data series or data to be evaluated.
#' @param obs Observed data series or data as reference.
#' @return Nash-Sutcliffe coefficient of efficiency of sim and obs series.
#' @export
NSCE <- function(sim, obs) {
  return(1 - NMSE(sim, obs))
}

#' Relative bias
#' @description Calculate relative bias (fraction, not percentage) of two data series.
#' @param sim Simulated data series or data to be evaluated.
#' @param obs Observed data series or data as reference.
#' @return Relative bias of efficiency of sim and obs series.
#' @export
BIAS <- function(sim, obs) {
  valid <- !is.na(sim) & !is.na(obs)
  sim <- sim[valid]
  obs <- obs[valid]
  return(mean(sim)/mean(obs) - 1)
}

#' Statistics of each respective column
#' @description Returns the statistic values, such as sum, mean, max,
#' min of each column of dataframe or matrix.
#' @param x a matrix of dataframe
#' @param na.rm	logical. Should missing values (including NaN) be removed?
#' @return A series of statistic value of each column of x.
#' @name statistic_col
#' @export
sum_col <- function(x, ...) {
  apply(x, 2, sum, ...)
}


#' @rdname statistic_col
#' @export
mean_col <- function(x, ...) {
  apply(x, 2, mean, ...)
}


#' @rdname statistic_col
#' @export
max_col <- function(x, ...) {
  apply(x, 2, max, ...)
}


#' @rdname statistic_col
#' @export
min_col <- function(x, ...) {
  apply(x, 2, min, ...)
}


#' @rdname statistic_col
#' @export
var_col <- function(x, ...) {
  apply(x, 2, var, ...)
}


#' @rdname statistic_col
#' @export
sd_col <- function(x, ...) {
  apply(x, 2, sd, ...)
}


#' write.table without col names and row names
#' @description The original write.table set col.names and row.names TRUE,
#' and when we want to merely write out a pure dataset or matrix we should
#' write in a long line: "col.names = F, row.names = F" and it so inconvenient
#' and troublesome. Now have this function we could not write it again.
#' Other usage is as same as write.table.
#' @export
wt <- function(x, file = "", ...) {
  write.table(x, file, row.names = FALSE, col.names = FALSE, ...)
}


#' Get benzi from nh zhan.
#' @param aurl Lianjie of benzi.
#' @return Nai
#' @description Hei hei hei
#' @import rvest downloader
#' @export
get_nh <- function(aurl, dir = "", mustnewdir=FALSE) {
  hasrvest <- require('rvest')
  hasdwer <- require('downloader')
  if(!hasrvest | !hasdwer) stop('rvesr or downloader was not installed.')

  lendir = nchar(dir)
  if(dir != "" & substr(dir, lendir, lendir) != "/")
    dir = paste(dir, "/", sep="")

  apg=read_html(aurl)

  title <- apg %>% html_node('div[id="info-block"]')%>% html_node('div[id="info"]') %>% html_node('h1') %>% html_text()
  print(paste("Collecting", title, "..."))

  pgnum <- length(apg %>% html_node('div[id="content"]') %>% html_node('div[id="thumbnail-container"]') %>% html_nodes('div[class="thumb-container"]'))
  print(paste("Total", pgnum, "pgs."))

  title <- gsub('[:\\*\\\\\\"\\?\\|\\/<>]', ' ', x=title)
  ifmkdir <- dir.create(paste(dir, title, sep = ""), recursive=F)
  if(!ifmkdir){
    if(mustnewdir){
      return
    } else {
      warning("Warning: dir exist.")
    }
  }

  i <- 1
  while(i <= pgnum){
    suburl=paste(aurl,i,sep = '')

    getsucceed=F
    while(!getsucceed){
      tryCatch({
        subpg <- read_html(suburl,timeout=1000)

        pgurl <-  subpg %>% html_node('div[id="content"]') %>% html_node('div[id="page-container"]') %>% html_node('section[id="image-container"]') %>% html_node('img') %>% html_attr('src')
        if(substr(pgurl, 1, 6) != "https:")
          pgurl <- paste('https:', pgurl, sep = '')

        ntmp <-  strsplit(pgurl, '/')[[1]]
        pgname <- ntmp[length(ntmp)]
        pg_path <- paste(dir, title, '/', pgname, sep='')
        download(pgurl, pg_path, mode='wb')
        print(paste("Saving pg to", pg_path))

        getsucceed <- T
      },
      error = function(e){
        print(paste('Error:', e))
        getsucceed <- F
      })
    }

    i <- i + 1
  }

}

#' Get the deficiency of the time series.
#' @param ts Time series to be find the lose time.
#' @return The time stamp lost in the time series.
#' @export
loss_time <- function(ts){
  ts <- sort(ts)
  len <- length(ts)
  int_ts <- as.integer(as.POSIXct(ts))
  diss <- as.numeric(int_ts[2:len] - int_ts[1:(len - 1)])
  mdis <- as.numeric(names(which.max(table(diss))))
  lp <- which(diss != mdis)
  ll <- diss[diss != mdis]
  loss <- c()

  if(length(ll) == 0)
    return(ll)
  for (i in 1:length(lp)) {
    if(ll[i] == 0)break
    stime <- int_ts[lp[i]]
    for (g in seq(mdis, ll[i], by=mdis)) {
      ls <- stime + g
      loss <- append(loss, ls)
    }
  }
  loss <- as.POSIXct(loss, origin="1970-01-01")
  return(loss)
}

#' Read IMERG data from netCDF4 file.
#' @description Only for the nc4 format IMERG data. Offering a series of path
#'              of the IMERG files, or single file path, and the coordinates
#'              of the data.
#' @param files Paths of the IMERG files to be read.
#' @param x Longitudes of the points where get the precipitation data.
#' @param y Latitudes of the points where get the precipitation data. Length
#'          should be as same as x.
#' @param varsn Order of which variable to be read. Default the first variable.
#'
#' @param verbose If TRUE, it will print the files read. Default FALSE.
#' @return A data frame of the IMERG data, for each column means a point, and
#'         each row means a time step. If only offered one point or one file,
#'         will return a vector of the data.
#' @export
readIMERG_nc4 <- function(files, x, y, varsn=1, verbose=FALSE) {
  if(length(x) != length(y)){
    stop("Lengths of x and y are not equal.")
  }

  l <- min(x)
  r <- max(x)
  b <- min(y)
  t <- max(y)

  f <- files[1]
  nc <- nc_open(f)
  lon <- nc$dim$lon$vals
  lat <- nc$dim$lat$vals
  nc_close(nc)

  rows <- rep(0, length(x))
  cols <- rep(0, length(y))
  for(i in 1:length(x)){
    col <- which.min(abs(x[i]-lon))
    row <- which.min(abs(y[i]-lat))
    cols[i] <- col
    rows[i] <- row
  }

  lc <- which.min(abs(l-lon))
  rc <- which.min(abs(r-lon))
  br <- which.min(abs(b-lat))
  tr <- which.min(abs(t-lat))

  ncol <- rc-lc+1
  nrow <- tr-br+1

  srows <- rows-br+1
  scols <- cols-lc+1
  dat <- data.frame(matrix(0, ncol=length(x), nrow=length(files)))

  for(d in 1:length(files)) {
    if(verbose) {
      print(paste("Reading", files[d]))
    }
    nc=nc_open(files[d])
    pvar <- nc$var[[varsn]]
    value <- ncvar_get(nc, pvar, c(br, lc), c(nrow, ncol))

    for(i in 1:length(x)) {
      dat[d, i] <- value[srows[i], scols[i]]
    }
    nc_close(nc)
  }
  if(nrow(dat) == 1)
    dat <- as.numeric(dat)
  if(ncol(dat) == 1)
    dat <- as.numeric(dat)

  return(dat)
}

#' Extract the sub strings of a single string or a string vector by regex.
#' @description Extrect the sub strings by offering a pattern, which can
#'              be regular expression in R.
#' @param pattern A string or regular expression. The same parameter in function
#'                gregexpr.
#' @param text A character vector to be extract sub string.
#' @param ... Other parameters of function gregexpr.
#'
#' @return A character vector when offering a single text, or a list containing
#'         the matched sub strings of each texts. If each texts has the same
#'         number of matched sub strings, it would be a multi-demantion vector.
#' @export
group <- function(pattern, text, ...) {
  rgx <- gregexpr(pattern, text, ...)
  mapply(function(rg, txt) {
    start <- as.integer(rg)
    end <- start + attr(rg, 'match.length') - 1
    mapply(function(s, e) {
      substr(txt, s, e)
    }, start, end)
  }, rgx, text)
}

#' Probability of the detection
#' @description Probability of the detection, like precipitation.
#' @param sim Series to be evaluated.
#' @param obs Reference series.
#' @param threshold Threshold to determind the event if happened. Default 0.
#' @return POD of the evaluated series and the reference series.
#' @export
POD <- function(sim, obs, th=0) {
  valid <- !is.na(sim) & !is.na(obs)
  sim <- sim[valid]
  obs <- obs[valid]
  n11 <- length(sim[sim > th & obs > th])
  n <- length(sim[obs > th])
  return(n11 / n)
}

#' False alarm ratio
#' @description False alarm ratio of the estimation, like precipitation.
#' @param sim Series to be evaluated.
#' @param obs Reference series.
#' @param threshold Threshold to determind the event if happened. Default 0.
#' @return FAR of the evaluated series and the reference series.
#' @export
FAR <- function(sim, obs, th=0) {
  valid <- !is.na(sim) & !is.na(obs)
  sim <- sim[valid]
  obs <- obs[valid]
  n01 <- length(sim[sim > th & obs <= th])
  n <- length(sim[sim > th])
  return(n01 / n)
}

#' Save urls to be download as meta4 format.
#' @description Save the url list of a series of file to be downloaded as
#'              .meta4 format for downthemall. No one can doubt that downthemall
#'              is a good add-on to make download more eaiser, convinient and
#'              managable. However, there are also some shortcommings to
#'              make it shit, for example, it would automatically remove the contents
#'              after square brackets without any consults when batch import
#'              url list, but those contents are necessary for some file urls such as
#'              OpenDAP. This function is to save the urls as meta4 file so
#'              that downthemall could not amend the urls when import them to
#'              the queue.
#' @param urls Url list of the file to be downloaded, such as OpenDAP.
#' @param file either a character string naming a file or a connection open
#'             for writing. "" indicates output to the console.
#'
#' @export
write.meta4 <- function(urls, file = "") {
  urls <- paste(urls)
  filenames <- sapply(strsplit(urls,'[/?]'), function(x) {
    pfn <- x[grepl("\\.",x)]
    if(length(pfn) >= 2) {
      pfn[2]
    } else {
      pfn[1]
    }
  })
  mtlg <- paste('<file name="',
               filenames,'" a0:num="',1:length(urls),
               '" a0:startDate="1488501103335"><url priority="100" a0:usable="',
               urls, '">', urls, '</url></file>')
  mtlh <- paste('<?xml version="1.0"?><metalink xmlns="urn:ietf:params:xml:ns:metalink" ',
                'version="4.0" a0:version="3.0.8" ',
                'xmlns:a0="http://www.downthemall.net/properties#"><generator>',
                'DownThemAll!/3.0</generator><published>',
                'Sat, 01 October 1949 00:08:00 GMT</published>', sep = '')
  mtl <- c(mtlh, mtlg, '</metalink>')
  write.table(mtl, file, quote = FALSE, col.names = F, row.names = F)
}


#' Calculate shortwave radiation using sunshine duration.
#' @description Calculate shortwave radiation using sunshine duration using
#'              Angstrom formula at daily scale.
#' @param ssd Sunshine duration (hour) in a day.
#' @param lat Latitude of the station location.
#' @param ts Date or date series of the ssd data. Would regard the first
#'           day in a year as the first ssd data if not provided.
#' @param a Coefficient of the Angstrom formula. Determine the relationship
#'          between ssd and radiation. Default 0.25.
#' @param b Coefficient of the Angstrom formula as a.
#'
#' @references Martinezlozano J A, Tena F, Onrubia J E, et al. The historical
#'             evolution of the Angstrom formula and its modifications: review
#'             and bibliography.[J]. Agricultural & Forest Meteorology, 1985,
#'             33(2):109-128.
#' @export
radiation <- function(ssd, lat, ts = NULL, a = 0.25, b=0.5) {
  lat <- lat*pi/180
  if(is.null(ts)) {
    J <- rep_len(c(1:365, 1:365, 1:365, 1:366), length(ssd))
  } else {
    J <- yday(ts)
  }
  dr <- 1 + 0.033*cos(pi*J/182.5)
  delta <- 0.408 * sin(pi*J/182.5 - 1.39)
  ws <- acos(-tan(lat) * tan(delta))

  Ra <- 1367 * dr/pi*(ws*sin(lat)*sin(delta)+cos(lat)*cos(delta)*sin(ws))

  N <- ws*24/pi
  Rs <- (a + b * ssd / N) * Ra
  Rs
}
