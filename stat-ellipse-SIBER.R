#' Compute normal confidence ellipses
#'
#' The method for calculating the ellipses has been modified from
#' `car::ellipse` (Fox and Weisberg, 2011)
#'
#' @references John Fox and Sanford Weisberg (2011). An {R} Companion to
#'   Applied Regression, Second Edition. Thousand Oaks CA: Sage. URL:
#'   \url{http://socserv.socsci.mcmaster.ca/jfox/Books/Companion}
#' @param level The confidence level at which to draw an ellipse (default is 0.95),
#'   or, if `type="euclid"`, the radius of the circle to be drawn.
#' @param type The type of ellipse.
#'   The default `"t"` assumes a multivariate t-distribution, and
#'   `"norm"` assumes a multivariate normal distribution.
#'   `"euclid"` draws a circle with the radius equal to `level`,
#'   representing the euclidean distance from the center.
#'   This ellipse probably won't appear circular unless `coord_fixed()` is applied.
#' @param segments The number of segments to be used in drawing the ellipse.
#' @inheritParams layer
#' @inheritParams geom_point
#' @export
#' @examples
#' ggplot(faithful, aes(waiting, eruptions)) +
#'   geom_point() +
#'   stat_ellipse()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3)) +
#'   geom_point() +
#'   stat_ellipse()
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3)) +
#'   geom_point() +
#'   stat_ellipse(type = "norm", linetype = 2) +
#'   stat_ellipse(type = "t")
#'
#' ggplot(faithful, aes(waiting, eruptions, color = eruptions > 3)) +
#'   geom_point() +
#'   stat_ellipse(type = "norm", linetype = 2) +
#'   stat_ellipse(type = "euclid", level = 3) +
#'   coord_fixed()
#'
#' ggplot(faithful, aes(waiting, eruptions, fill = eruptions > 3)) +
#'   stat_ellipse(geom = "polygon")
stat_ellipseSIBER <- function(mapping = NULL, data = NULL,
                         geom = "path", position = "identity",
                         ...,
                         type = "SEAc",
                         level = 0.95,
                         segments = 100,
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatEllipseSIBER,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      level = level,
      segments = segments,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatEllipseSIBER <- ggproto("StatEllipseSIBER", Stat,
                       required_aes = c("x", "y"),
                       
                       compute_group = function(data, scales, type = "SEAc", 
                                                level = 0.95,
                                                segments = 100, na.rm = FALSE) {
                         calculate_ellipse_SIBER(data = data, 
                                                 vars = c("x", "y"), 
                                                 type = type,
                                           level = level, segments = segments)
                       }
)

calculate_ellipse_SIBER <- function(data, vars, type, level, segments){
  
  
  # sample size
  n <- nrow(data)
  
  # check we have enough data points. At least 3.
  if (n < 3) {
    message("Too few points to calculate an ellipse")
    ellipse <- rbind(as.numeric(c(NA, NA)))}
  
  # create the ellipse here
  
  # the covariance ellipse
  ss <- cov.wt(data[,vars])
  
  sigma  <- ss$cov
  center <- ss$center
  
  
  
  # if ci.mean is F (default) then we are plotting quantiles of the sample, 
  # i.e. prediction ellipses, and so we set c <- 1 so that it has no effect
  # below. Else it divides the radius calculation below by sqrt(m) to include
  # the conversion from standard deviation to standard error of the mean.
  ifelse(type == "mean", 
         c.scale <- n,
         c.scale <- 1
  )
  
  # The small.sample toggles on and off (default) the small sample size 
  # correction to essentially plot the SEAc in place of the SEA. It can be 
  # used inconjuction with any prediction ellipse.
  if(type == "SEAc") q <- (n - 1) / (n - 2)
    
  # Dont apply small sample size correction if only SEA
  if(type == "SEA")  q <- 1
  
  
  # if p is NULL then plot a standard ellipse with r = 1
  # else generate a prediction ellipse that contains
  # approximately proportion p of data by scaling r
  # based on the chi-squared distribution.
  # p defaults to NULL.
  ifelse(is.null(level), 
         r <- 1, 
         r <- sqrt(stats::qchisq(level, df=2))
  )
  
  
  # get the eigenvalues and eigenvectors of sigma
  # if ci.mean = T then the covariance matrix is divided by the sample size
  # so as to produce confidence ellipses for the mean. Else it has no 
  # effect with c.scale = 1.
  e = eigen(sigma / c.scale)
  
  # 
  SigSqrt <- e$vectors %*% diag(sqrt(e$values * q)) %*% t(e$vectors)
  
  # create a unit radius circle to transform
  cc <- SIBER:::genCircle(segments, r)
  
  # transform the unit circle according to the covariance 
  # matrix sigma
  
  # a function to transform the points
  back.trans <- function(x) {
    return(SigSqrt %*% x + center)
  }
  
  # apply the transformation to calculate the ellipse.
  
  ellipse = t(apply(cc, 1, back.trans))
  
  ellipse <- as.data.frame(ellipse)
  colnames(ellipse) <- vars
  ellipse
}