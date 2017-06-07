# Imported from http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
# Accessed on 8 Jun 2016

#' @name spatial
#' @title Great-circle distance calculations in R
#' @author Mario Pineda-Krch
#' @source \url{http://www.r-bloggers.com/great-circle-distance-calculations-in-r/}
#' 
#' @description The functions to calculate the spatial distance (meters) between 
#' a large number of longitude and latitude locations.
#' 
#' @details Converts degrees to radians
#' 
#' @param deg The degrees.
#' @export
#' @keywords spatial
#' @examples 
#' lat1 <- deg2rad(deg)
#' 
#' @rdname spatial
deg2rad <- function(deg) return(deg*pi/180)

#' @details \code{gcd.slc} calculates the geodesic distance (meters) between two points  
#' specified by radian latitude/longitude using the spherical Law of Cosines (slc).
#' 
#' @param long1,lat1,long2,lat2 radian longitude or latitude of two points. 
#' @export
#' @examples 
#' long1 <- deg2rad(-36.8485)
#' lat1 <- deg2rad(174.7633)
#' long2 <- deg2rad(-37.7870)
#' lat2 <- deg2rad(175.2793)
#' # Auckland <=> Hamilton
#' dist.slc <- gcd.slc(long1, lat1, long2, lat2)
#' 
#' @rdname spatial
gcd.slc <- function(long1, lat1, long2, lat2) {
  RAD <- 6371000 # Earth mean radius [meter]
  #d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  # use pmin to avoid acos bug http://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  d <- acos( pmin(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1), 1.0) ) * RAD
  return(d) # Distance in meter
}

#' @details \code{gcd.hf} calculates the geodesic distance (meters) between two points 
#' specified by radian latitude/longitude using the Haversine formula (hf).
#' @export
#' @examples 
#' dist.hf <- gcd.hf(long1, lat1, long2, lat2)
#' 
#' @rdname spatial
gcd.hf <- function(long1, lat1, long2, lat2) {
  RAD <- 6371000 # Earth mean radius [meter]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = RAD * c
  return(d) # Distance in meter
}

#' @details \code{gcd.vif} calculates the geodesic distance (meters) between two points 
#' specified by radian latitude/longitude using Vincenty inverse formula for ellipsoids (vif).
#' @export
#' @examples 
#' dist.vif <- gcd.vif(long1, lat1, long2, lat2)
#' 
#' @rdname spatial
gcd.vif <- function(long1, lat1, long2, lat2) {
  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # ength of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid
  
  L <- long2-long1 # difference in longitude
  U1 <- atan((1-f) * tan(lat1)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)
  
  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL
  
  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                        (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
      (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) - 
                                             B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  s <- b*A*(sigma-deltaSigma)
  
  return(s) # Distance in meter
}

#' @details \code{gcd.slc} calculates the geodesic distance (meters) between two points 
#' specified by degrees latitude/longitude using Haversine formula (hf), 
#' Spherical Law of Cosines (slc) and Vincenty inverse formula for ellipsoids (vif).
#' 
#' @param long1.dg,lat1.dg,long2.dg,lat2.dg degrees longitude or latitude of two points. 
#' @export
#' @examples 
#' dist.vif <- gcd(-36.8485, 174.7633, -37.7870, 175.2793, distance="vif")
#' 
#' @rdname spatial
gcd <- function(long1.dg, lat1.dg, long2.dg, lat2.dg, distance=c("slc", "hf", "vif")) {
  # Convert degrees to radians
  long1 <- deg2rad(long1.dg)
  lat1 <- deg2rad(lat1.dg)
  long2 <- deg2rad(long2.dg)
  lat2 <- deg2rad(lat2.dg)
  
  distance <- match.arg(distance)
  if (distance=="hf") {
    dist.m <- gcd.hf(long1, lat1, long2, lat2)
  } else if (distance=="vif") {
      dist.m <- gcd.vif(long1, lat1, long2, lat2)
  } else {
    dist.m <- gcd.slc(long1, lat1, long2, lat2)
  }
  return(dist.m) # Distance in meter
}


#' @details \code{geodesicDist.slc} is a batch processing function to use \code{gcd.slc}.
#' 
#' @param lati.long.degrees A 2-column data frame contains degrees latitude/longitude coordinates.
#' The 1st col is latitude and 2nd is longitude.
#' @param file Either a character string naming a file or a connection open for writing. 
#' "" indicates output to the console. If NA, as default, then do nothing.
#' @param sep The field separator string. 
#' @export
#' @keywords spatial
#' @examples 
#' dists.btw <- geodesicDist.slc(env.byplot[, c("Latitude","Longitude")])
#' 
#' @rdname spatial
geodesicDist.slc <- function(lati.long.degrees, file = NA, sep = "\t") {
  # latitude/longitude coordinates in a 2-col data frame
  n.coord <- nrow(lati.long.degrees)
  
  dlist <- list()
  n <- 1
  for(i in 1:(n.coord-1)){
    for(j in 2:n.coord){
      lat1 <- deg2rad(lati.long.degrees[i,1])
      long1 <- deg2rad(lati.long.degrees[i,2])
      lat2 <- deg2rad(lati.long.degrees[j,1])
      long2 <- deg2rad(lati.long.degrees[j,2])
      
      d <- gcd.slc(long1, lat1, long2, lat2)
      #cat("i = ", i, "j = ", j, ", dist = ", d, "\n")
      
      res <- list("plot1" = rownames(lati.long.degrees)[i], "plot2" = rownames(lati.long.degrees)[j], 
                  "lat1" = lat1, "lat2" = lat2, "long1" = long1, "long2" = long2, "dist" = d)
      dlist[[n]] <- res
      n <- n + 1
    }
  }
  
  require(data.table)
  dists.btw <- rbindlist(dlist)
  if (!is.na(file))
    write.table(dists.btw, file = file, sep = sep, quote = F, col.names = NA)
  return(dists.btw)
}

#' @details \code{geodesicDist.hf} is a batch processing function to use \code{gcd.hf}.
#'
#' @export
#' @keywords spatial
#' @examples 
#' dists.btw <- geodesicDist.hf(env.byplot[, c("Latitude","Longitude")])
#' 
#' @rdname spatial 
geodesicDist.hf <- function(lati.long.degrees, file = NA, sep = "\t") {
  # latitude/longitude coordinates in a 2-col data frame
  n.coord <- nrow(lati.long.degrees)
  
  # Distance calculations
  RAD <- 6371000 # Earth mean radius [m]
  dlist <- list()
  n <- 1
  for(i in 1:(n.coord-1)){
    for(j in 2:n.coord){
      lat1 <- deg2rad(lati.long.degrees[i,1])
      long1 <- deg2rad(lati.long.degrees[i,2])
      lat2 <- deg2rad(lati.long.degrees[j,1])
      long2 <- deg2rad(lati.long.degrees[j,2])
      
      d = gcd.hf(long1, lat1, long2, lat2)
      #cat("i = ", i, "j = ", j, ", dist = ", d, "\n")
      
      res <- list("plot1" = rownames(lati.long.degrees)[i], "plot2" = rownames(lati.long.degrees)[j], 
                  "lat1" = lat1, "lat2" = lat2, "long1" = long1, "long2" = long2, "dist" = d)
      dlist[[n]] <- res
      n <- n + 1
    }
  }
  
  require(data.table)
  dists.btw <- rbindlist(dlist)
  if (!is.na(file))
    write.table(dists.btw, file = file, sep = sep, quote = F, col.names = NA)
  return(dists.btw)
}
