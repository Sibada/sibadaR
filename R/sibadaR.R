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

#' Mann-Kendall趋势分析
#'
#' @param x 待分析的时间序列，可以是数据框或者是矩阵形式
#' @return 如果x为单个序列就返回一个列表，其中MK.value为M-K统计值，p.value为显著性达到了多少置信度。
#'         如果x是个数据框或者矩阵就挨列计算，返回一个矩阵记录每一列的MK.value和p.value。
#'
#' @export
MK_test <- function (x) {
  if (is.null(dim(x)))
    return(as.list(mk(x)))
  if (ncol(x) == 1)
    return(as.list(mk(x[, 1])))
  return(t(apply(x, 2, mk)))
}


#' 对xts时间序列按季节（如按春夏秋冬或者全部季节）应用函数。
#'
#' @description 按3~5月为春季，6～8月为夏季，9~11月为秋季，12～次年2月为冬季的方法算。
#' @param x xts时间序列
#' @param FUN 某R函数
#'
#' @return 对该时间序列按季节应用FUN产生的结果
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
