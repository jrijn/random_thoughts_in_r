# This is the pracma/circlefit.R script from GitHub (https://github.com/cran/pracma/blob/master/R/circlefit.R),
# but then adapted to return the RMS error as the fourth element of the output vector.
# i.e.the output is: c(xcoord, ycoord, radius, RMSerror).
circlefits <- function(xp, yp, fast = FALSE) {
  if (!is.vector(xp, mode="numeric") || !is.vector(yp, mode="numeric"))
    stop("Arguments 'xp' and 'yp' must be numeric vectors.")
  if (length(xp) != length(yp))
    stop("Vectors 'xp' and 'yp' must be of the same length.")
  
  n  <- length(xp)
  p <- qr.solve(cbind(xp, yp, 1), matrix(xp^2 + yp^2, ncol = 1))
  r <- c(p[1]/2, p[2]/2, sqrt((p[1]^2 + p[2]^2)/4 + p[3]))
  
  cerr <- function(v)
    sqrt(sum((sqrt((xp - v[1])^2 + (yp - v[2])^2) - v[3])^2)/n)
  
  if (fast) {
    RMS <- cerr(r)
  } else {
    q <- optim(p, cerr)
    RMS <- q$value
    r <- q$par
  }
  
  r <- append(r, RMS)
  return(r)
}

# This function is a wrapper function, which returns NA if the circlefits function fails.
# Now you can use this in a dplyr %>% summarise().
circlefitter <- function(xp, yp) {
  out <- tryCatch(
    {
      circlefits(xp, yp)
    },
    error = function(cond) {
      return(NA)
    }
  )
  return(out)
}

library(tidyverse)

spots <- read_csv("Spots in tracks statistics.csv")

spots <- spots %>%
  group_by(TRACK_ID) %>%
  mutate(X_NORM = POSITION_X - POSITION_X[1],
         Y_NORM = POSITION_Y - POSITION_Y[1])

circles <- spots %>%
  group_by(TRACK_ID) %>%
  summarise(X = circlefitter(X_NORM, Y_NORM)[1],
            Y = circlefitter(X_NORM, Y_NORM)[2],
            RADIUS = circlefitter(X_NORM, Y_NORM)[3],
            RMSERROR = circlefitter(X_NORM, Y_NORM)[4])
