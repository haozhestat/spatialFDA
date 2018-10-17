#distance on earth
earth_dist <- function (location1, location2) 
{
  long1 <- as.numeric(location1$lon) # location1=(longitude,latitude)
  lat1 <- as.numeric(location1$lat)
  long2 <- as.numeric(location2$lon)
  lat2 <- as.numeric(location2$lat)
  rad <- pi/180
  a1 <- lat1 * rad
  a2 <- long1 * rad
  b1 <- lat2 * rad
  b2 <- long2 * rad
  dlon <- b2 - a2
  dlat <- b1 - a1
  a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  R <- 6378.145
  d <- R * c
  return(d)
}