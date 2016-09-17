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
