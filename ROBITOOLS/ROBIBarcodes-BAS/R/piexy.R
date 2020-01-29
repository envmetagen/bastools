#'@include ROBIBarcodes.R
NULL

#' @export
pie.xy = function (x,y,data, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
          init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
          col = NULL, border = NULL, lty = NULL, ...) 
{
  if (!is.numeric(data) || any(is.na(data) | data < 0)) 
    stop("'data' values must be positive.")
  if (is.null(labels)) 
    labels <- as.character(seq_along(data))
  else labels <- as.graphicsAnnot(labels)
  data <- c(0, cumsum(data)/sum(data))
  dx <- diff(data)
  nx <- length(dx)
#  plot.new()
#  pin <- par("pin")
  xlim <- ylim <- c(-1, 1)
#   if (pin[1L] > pin[2L]) 
#     xlim <- (pin[1L]/pin[2L]) * xlim
#   else ylim <- (pin[2L]/pin[1L]) * ylim
#  dev.hold()
#  on.exit(dev.flush())
#  plot.window(xlim, ylim, "", asp = 1)
  if (is.null(col)) 
    col <- if (is.null(density)) 
      c("white", "lightblue", "mistyrose", "lightcyan", 
        "lavender", "cornsilk")
  else par("fg")
  if (!is.null(col)) 
    col <- rep_len(col, nx)
  if (!is.null(border)) 
    border <- rep_len(border, nx)
  if (!is.null(lty)) 
    lty <- rep_len(lty, nx)
  angle <- rep(angle, nx)
  if (!is.null(density)) 
    density <- rep_len(density, nx)
  twopi <- if (clockwise) 
    -2 * pi
  else 2 * pi
  t2xy <- function(t) {
    t2p <- twopi * t + init.angle * pi/180
    list(x = radius * cos(t2p) , y = radius * sin(t2p))
  }
  for (i in 1L:nx) {
    n <- max(2, floor(edges * dx[i]))
    P <- t2xy(seq.int(data[i], data[i + 1], length.out = n))
    if (nx>1)
      polygon(c(P$x + x, x), c(P$y + y , y), 
              density = density[i], angle = angle[i], 
              border = border[i], col = col[i], lty = lty[i])
    else
      polygon(P$x + x, P$y + y, 
              density = density[i], angle = angle[i], 
              border = border[i], col = col[i], lty = lty[i])
    P <- t2xy(mean(data[i + 0:1]))
    lab <- as.character(labels[i])
    if (!is.na(lab) && nzchar(lab)) {
      lines(c(1 , 1.05) * P$x + x, c(1, 1.05) * P$y + y)
      text(1.1 * P$x + x , 1.1 * P$y + y, labels[i], xpd = TRUE, 
           adj = ifelse(P$x < 0, 1, 0), ...)
    }
  }
#  title(main = main, ...)
  invisible(NULL)
}