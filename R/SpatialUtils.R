# Author: Andrew Dopheide
# Accessed on 15 Sep 2016

#' @name spatial
#' @title Utils related to spatial distances
#'
#' @description The functions related to spatial distances.
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

#' @details \code{geodesicDist.slc} calculates the geodesic distance (km) 
#' between two points specified by radian latitude/longitude 
#' using the Spherical Law of Cosines (slc).
#' @source \url{http://www.r-bloggers.com/great-circle-distance-calculations-in-r/}
#' 
#' @param lati.long A 2-column data frame contains latitude/longitude coordinates.
#' @param file Either a character string naming a file or a connection open for writing. 
#' "" indicates output to the console. If NA, as default, then do nothing.
#' @param sep The field separator string. 
#' @export
#' @keywords spatial
#' @examples 
#' dists.btw <- geodesicDist.slc(env.byplot[, c("Latitude","Longitude")])
#' 
#' @rdname spatial
geodesicDist.slc <- function(lati.long, file = NA, sep = "\t") {
  # latitude/longitude coordinates in a 2-col data frame
  n.coord <- nrow(lati.long)
  
  # Distance calculations
  RAD <- 6371000 # Earth mean radius [m]
  dlist <- list()
  n <- 1
  for(i in 1:(n.coord-1)){
    for(j in 2:n.coord){
      lat1 <- deg2rad(lati.long[i,1])
      long1 <- deg2rad(lati.long[i,2])
      lat2 <- deg2rad(lati.long[j,1])
      long2 <- deg2rad(lati.long[j,2])
      # use pmin to avoid acos bug http://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
      d <- acos( pmin(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1), 1.0) ) * RAD
      #cat("i = ", i, "j = ", j, ", dist = ", d, "\n")
      res <- list("plot1" = rownames(lati.long)[i], "plot2" = rownames(lati.long)[j], 
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

#' @details \code{geodesicDist.hf} calculates the geodesic distance (km) 
#' between two points specified by radian latitude/longitude 
#' using the Haversine formula (hf).
#'
#' @export
#' @keywords spatial
#' @examples 
#' dists.btw <- geodesicDist.hf(env.byplot[, c("Latitude","Longitude")])
#' 
#' @rdname spatial 
geodesicDist.hf <- function(lati.long, file = NA, sep = "\t") {
  # latitude/longitude coordinates in a 2-col data frame
  n.coord <- nrow(lati.long)
  
  # Distance calculations
  RAD <- 6371000 # Earth mean radius [m]
  dlist <- list()
  n <- 1
  for(i in 1:(n.coord-1)){
    for(j in 2:n.coord){
      lat1 <- deg2rad(lati.long[i,1])
      long1 <- deg2rad(lati.long[i,2])
      lat2 <- deg2rad(lati.long[j,1])
      long2 <- deg2rad(lati.long[j,2])
      
      delta.long <- (long2 - long1)
      delta.lat <- (lat2 - lat1)
      a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
      c <- 2 * asin(min(1,sqrt(a)))
      d = RAD * c
      #cat("i = ", i, "j = ", j, ", dist = ", d, "\n")
      
      res <- list("plot1" = rownames(lati.long)[i], "plot2" = rownames(lati.long)[j], 
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
