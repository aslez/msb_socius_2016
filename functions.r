#DEFINE CUSTOM FUNCTIONS

gwr2df <- function(formula, spdf, bw = NULL, longlat = NULL) {
  if (is.null(bw))
    bw <- gwr.sel(formula, data = spdf@data, longlat = longlat,
                  coords = coordinates(spdf), verbose = FALSE)
  mod <- gwr(formula, data = spdf@data, longlat = longlat,
             coords = coordinates(spdf), bandwidth = bw)
  gwr2fema(mod)
}

shp2df <- function(shp, slopes) {
  shp@data[, c('coord.x', 'coord.y')] <- coordinates(shp)
  shp@data <- left_join(shp@data %>% select(coord.x, coord.y), slopes)
  shp$id <- as.character(1:NROW(shp@data))
  plot_dat <- select(shp@data, id, coord.x, coord.y, WHO)
  fortify(shp, region = 'id') %>%
    left_join(plot_dat)
}

fema <- function(shp_lst, opt = FALSE, longlat = TRUE) {
  #identify bandwidth
  bw_lst <- lapply(shp_lst, function(z) {
    bw <- gwr.sel(VOTE ~ WHO, data = z@data, longlat = longlat,
                  coords = coordinates(z), verbose = FALSE)
  })
  mean_bw <- mean(do.call('rbind', bw_lst))
  print(paste0('Mean Bandwidth = ', mean_bw))

  #generate field measures
  if (opt) df_lst <- lapply(shp_lst, gwr2df,
                            longlat = longlat,
                            formula = VOTE ~ WHO)
  else df_lst <- lapply(shp_lst, gwr2df,
                        formula = VOTE ~ WHO,
                        longlat = longlat,
                        bw = mean_bw)
  field_lst <- lapply(df_lst, function(z) {
    summary(discorr(z, longlat = longlat))
    })

  measure_lst <- lapply(field_lst, function(x) {
    unlist(x[c('range', 'std.effect', 'strength')])
  })

  measure_df <- data.frame(do.call('rbind', measure_lst)) %>%
    mutate(Year = seq(1890, 1896, 2)) %>%
    gather(measure, Value, -Year) %>%
    mutate(measure = plyr::revalue(measure,
                                   c('range' = 'Range',
                                     'std.effect' = 'Organization',
                                     'strength' = 'Strength')))

  #generate mapping data
  wrk_lst <- shp_lst
  yvec <- seq(1890, 1896, 2)
  fort_lst <- lapply(seq_along(wrk_lst), function(z) {
    result <- shp2df(wrk_lst[[z]], df_lst[[z]])
    result$YEAR <- yvec[z]
    result
  })
  fort_df <- do.call(rbind, fort_lst)

  #generate centroids
  #use 10 and 20 when working with projected coordinates
  cent_df <- fort_df %>%
    select(group, coord.x, coord.y, WHO, YEAR) %>%
    unique() %>%
    mutate(angle = atan(WHO),
           xe = coord.x + .15 * cos(angle),
           ye = coord.y + .15 * sin(angle),
           xs = xe + .3 * cos(angle + pi),
           ys = ye + .3 * sin(angle + pi))

  #results
  list(measures = measure_df,
       shapes = fort_df,
       centroids = cent_df)
}

# summarize FEMA objects
summary.fema <- function(x) {
  total.effect <- sum(x$rmat * x$dmat)
  range <- total.effect / sum(x$rmat)
  std.effect <- total.effect / sum(x$dmat)
  bdis <- apply(x$slopes, 1, function(z) sqrt(sum(z^2)))
  strength <- std.effect * mean(bdis)
  result <- list(coords = x$coords, slopes = x$slopes,
                 rmat = x$rmat, dmat = x$dmat,
                 total.effect = total.effect,
                 range = range,
                 std.effect = std.effect,
                 strength = strength)
  class(result) <- c('summary.fema', 'fema')
  result
}

discorr <- function(df, longlat = FALSE, ctype = 'single1') {
  if (!'fema.df' %in% class(df)) {
    stop('df is not a fema-class data frame')
  }

  #generate distance matrix
  coords <- cbind(df$coord.x, df$coord.y)
  dmat <- spDists(coords, longlat = longlat)

  #generate correlations
  cloc <- sapply(c('coord.x', 'coord.y'), function(x) grep(x, names(df)))
  slopes <- data.frame(df[, -cloc])
  if (NCOL(slopes) > 1) ctype <- 'multi'
  names(slopes) <- names(df)[-cloc]
  rmat <- matrix(0, nrow(dmat), ncol(dmat))
  for (i in 2:nrow(slopes)) {
    for (j in 1:(i - 1)) {
      rmat[i, j] = pwc(slopes[i, ], slopes[j, ], ctype)
    }
  }
  if (ctype == 'multi') {
    rmat <- (max(rmat) - 2 * rmat) / max(rmat)
  }
  rmat <- rmat + t(rmat)

  #compile results
  result <- list(coords = coords, slopes = slopes,
                 longlat = longlat, ctype = ctype,
                 dmat = dmat, rmat = rmat)
  class(result) <- 'fema'
  result
}

pwc <- function(x1, x2, ctype) {
  x1x2 <- data.frame(x1 = x1, x2 = x2)
  if (ctype == 'single1') {
    x1x2 <- rbind(1, x1x2)
    cpx1 <- with(x1x2, crossprod(x1))
    cpx2 <- with(x1x2, crossprod(x2))
    cpx1x2 <- with(x1x2, crossprod(x1, x2))
    r <- cpx1x2 / sqrt(cpx1 * cpx2)
  }
  if (ctype == 'single2') r <- with(x1x2, (pi - 2 * abs(x1 - x2)) / pi)
  if (ctype == 'multi') r <- with(x1x2, sqrt(sum((x1 - x2)^2)))
  r
}

gwr2fema <- function(x, slopes = NULL, intercept = FALSE) {
  if (class(x) != 'gwr') {
    stop('x is not a gwr-class object')
  }
  cnames <- names(coef(x$lm))
  iloc <- grep('(Intercept)', cnames)
  if (!intercept) cnames <- cnames[-iloc]
  if (!is.null(slopes)) {
    sloc <- sapply(slopes, function(z) grep(z, cnames))
    if (intercept) cnames <- cnames[c(iloc, sloc)]
    else cnames <- cnames[sloc]
  }
  df <- as.data.frame(x$SDF[, cnames])
  class(df) <- c('fema.df', 'data.frame')
  df
}


# same as scale() but for some reason I did this instead
std <- function(x) {
  # standardizes a vector
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = T)
  stdv <- (x - m) / s
  return(stdv)
}

# convenience ftn for euclidean distance matrix
eDist <- function(x) {
  # x: vector for calculating distances
  d <- dist(x, diag = T)
  d <- as.matrix(d)
  d[upper.tri(d)] <- NA
  return(d)
}

# break observations into bins of size N, find mean, then plot
binPlot <- function(geo, edist, binsize = 1000, title, color="black", ...) {
  # Plots mean of equal-N bins.
  geov <- geo[lower.tri(geo)]
  edistv <- edist[lower.tri(edist)]
  dfall <- as.data.frame(cbind(geov, edistv))
  dfall <- dfall[order(geov), ]
  # plot by bins
  # define min, max, and increment for bins
  bins <- seq(0, length(geov), binsize)
  mnd <- c()          # mean euclidean distance per bin
  mng <- c()          # mean geographical distance per bin
  for (i in 2:length(bins)) {
    # mean euclidean distance in bin
    mnd <- c(mnd, mean(dfall[bins[i - 1]:bins[i], "edistv"],
      na.rm = TRUE))
    # mean geographic distance in bin
    mng <- c(mng, mean(dfall[bins[i - 1]:bins[i], "geov"],
      na.rm = TRUE))
  }
  plot(mng, mnd,
    type = "l", col = "black", lwd = 2,
    ylab = "Mean difference",
    xlab = "Geographical distance (miles)",
    main = title,
    ...
  )
}

# find distance between two points
geoDist <- function(lat1, lon1, lat2, lon2, units = "mi") {
  # This function ouputs the (great circles) distance between to geographic
  # points, using the Spherical Law of Cosines
  #
  # Args:
  #   lon1: longitude for point 1
  #   lat1: latitutde for point 1
  #   lon2: longitude for point 2 (number or vector)
  #   lat2: latitutde for point 2 (number or vector)
  #
  # Returns:
  #   The distance (in miles) between the two given points.

  # set units
  r <- 3959.873 # Earth mean radius (miles)
  if (units == "km") {
    r <- 6372.798
  }

  # convert decimial degree input to radian degrees :
  z <- (lat1 == lat2 & lon1 == lon2)  # logical vector
  lon1.rad <- lon1 * (pi / 180)
  lat1.rad <- lat1 * (pi / 180)
  lon2.rad <- lon2 * (pi / 180)
  lat2.rad <- lat2 * (pi / 180)

  distance <- (acos(sin(lat1.rad) * sin(lat2.rad) +
    cos(lat1.rad) * cos(lat2.rad) * cos(lon2.rad - lon1.rad)) * r)

  distance[z] <- 0  # Distance = 0 if lat and lon are identical

  return(distance) # Distance
}


# matrix of distances between many points
geoDistMatrix <-
  function(df, group=1, lon="lon", lat="lat", units = "mi") {
  # Calculates the geodesic distance between two points specified by radian
  # latitude/longitude using the Spherical Law of Cosines (slc).
  # Args:
  #   df: dataframe containing coordinates and group names
  #   group: unique identifier for each set of coordinates; column in
  #     dataframe to be used for dimension names.
  #   lon: longitude column name (decimal degree only!)
  #   lat: latitude column name (decimal degree only!)
  #   units: "mi" for miles, or "km" for kilometers
  #
  # Returns:
  #   An n by n matrix of distances between points (lower triangle values only).

  # Check for valid group, lon, and lat parameters
  if (is.null(df[[lat]]) | is.null(df[[lon]])) {
    stop("Invalid coordinate columns specified.")
  }
  if (is.null(df[[group]])) {
    stop("Invalid group column specified.")
  }

  # parse, clean, and prep data
  coords <- unique(cbind(df[c(names(df[group]), lon, lat)])) # no dups
  coords <- coords[order(coords[group]), ] # sort by group
  n <- nrow(coords)
  matrix.names = list(coords[[group]], coords[[group]])
  distance.matrix <- matrix(nrow = n, ncol = n, dimnames = matrix.names)

  # calculate distances
  if (n == 1) {
      return(distance.matrix)
  }
  for (i in 1:n) {
      j = 1:i
      distance.matrix[i, j] = suppressWarnings(geoDist(coords[i, lat], coords[i, lon],
      coords[j, lat], coords[j, lon], units = units))
      # distance.matrix[i, i] = 0  # replaces NA's with 0 for diagonal
  }
  return(distance.matrix)
}


# calculates social distance from df of variables
socDist <- function(df, group = "group", na.cases = "mean") {
  # function for calculating social distance
  # requires dataframe of ONLY the variables on which dist is to be
  #   calculated
  # Returns 'dist' type object

  ##### PREP DATA #####
  df <- unique(df)              # keep only unique combinations
  df <- sapply(df, as.numeric)  # convert factors to numeric

  # standardize each column
  df <- apply(df, 2, scale)

  ##### CALCULATE DISTANCE #####
  sdist <- dist(df)

  return(sdist)
}

# Creates a matrix of correlations for a set of slopes (vectors).
vCorr <- function(df, col = 1, group = NULL) {
  # Args:
  #   df: M by N dataframe or matrix of slopes
  #   col: column to use in df. Can be name (ex: "col") or column index
  #   group: column of group names, used for matrix names
  # Returns:
  #   matrix of vector correlations between slopes

  # Define internal correlation function for a pair
  vCorrPair <- function (x, y = NULL) {
    # Code borrowed from Fridolin Wild's cosine() function in 'lsa' library
    if (is.matrix(x) && is.null(y)) {
      co = array(0, c(ncol(x), ncol(x))) # create empty n x n matrix
      f = colnames(x)
      dimnames(co) = list(f, f) # give matrix cols and rows names
      for (i in 2:ncol(x)) {
        for (j in 1:(i - 1)) {
          co[i, j] = vCorrPair(x[, i], x[, j])
        }
      }
      co = co + t(co)
      diag(co) = NA
      return(as.matrix(co))
    }
    else if (is.vector(x) && is.vector(y)) {
      return(crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)))
    }
    else {
      stop("Argument mismatch. Either one matrix or two vectors needed as input.")
    }
  }

  # convert to dataframe, if needed
  if(!is.data.frame(df)) {
    print("Converting data to dataframe.")
    df <- as.data.frame(df)
  }

  # create 2 x N matrix for cosine similarity, each row (1, x)
  slope.vectors <- rbind(1, df[[col]])
  # add group names, if applicable
  if (!is.null(group)) {
    colnames(slope.vectors) <- df[[group]]  # name matrix rows/cols
  }
  # calculate cosine similarity/vector correlation matrix
  vcorr.matrix <- vCorrPair(slope.vectors)
  # Replace values above diag to NA
  for (i in 1:(nrow(vcorr.matrix)-1)) {
    vcorr.matrix[i, (i+1):ncol(vcorr.matrix)] <- NA
  }
  return(vcorr.matrix)
}

# Calculates FEMA multiple correlation, based on euclidean distance
distCorr <- function (df, group = 1, slopes = c(6:ncol(df)), zslopes = FALSE,
    stddist = TRUE) {
  # Arguments:
  #   df: dataframe of slopes (columns) by groups (rows), or single
  #     vector of slopes
  #   slopes: vector of slopes to use in calculation.
  #   group: Unique identifier for each set of coordinates; column in
  #   zslopes: standardize slopes for uniform scale
  #   stddist: standardize distances by number of valid slope
  #     comparisons (for dealing with missing slopes)
  # Returns:
  #   Matrix of correlations (bottom triangle values only)

  # Define distance function
  stdDist <- function (x, y = NULL, zslopes = FALSE, stddist = TRUE) {
    # custom distance function
    if (is.matrix(x) && is.null(y)) {
      if (zslopes == TRUE) {
        for (i in 1:ncol(x)) {
          m <- mean(x[, i], na.rm = T)
          s <- sd(x[, i], na.rm = T)
          x[,i] <- (x[, i] - m)/s
        }
      }
      co = array(NA, c(nrow(x), nrow(x)))
      f = rownames(x)
      dimnames(co) = list(f, f)
      for (i in 2:nrow(x)) {
        for (j in 1:(i - 1)) {
          co[i, j] = stdDist(x[i,], x[j,], stddist = stddist)
        }
      }
      # set diagonal values to NA
      diag(co) = NA
      return(as.matrix(co))
    }
    else if (is.vector(x) && is.vector(y)) {
      n <- 1
      if (stddist == TRUE) {
        n <- sum(!is.na(x - y))
      }
      ssq <- sum((x - y)^2, na.rm=T)
      d <- sqrt(ssq/n)
      return(d)
    }
    else {
      stop("Argument mismatch. Either one matrix or two vectors needed as input.")
    }
  }

  x <- as.matrix(df[,slopes]) # convert dataframe to matrix of only valid slopes
  rownames(x) <- df[[1]]  # give rownames to matrix

  dmx <- stdDist(x, zslopes = zslopes, stddist = stddist) # calculate standardized distance
  mx <- max(dmx, na.rm=T) # define maximum distance

  # correlation calculation
  corrs <- (mx - (2 * dmx)) / mx

  return(corrs)
}

# regress by group and time
femaRegress <-
  function(df, group="group", model="y ~ x1", timecol=NULL,
    lat=NULL, lon=NULL, ...) {
  # Performs time-aware group-level regressions.
  # Returns dataframe with group, time, coordinates, and coefficients
  #
  # Args:
  #   df: dataframe to be used
  #   group: column in dataframe by which data is nested
  #   model: the actual model for running regression,
  #     for example, 'y ~ x1 + x2'
  #   timecol: column in dataframe specifying time (optional)
  #   lat: latitude column in dataframe (optional)
  #   lon: longitude column in dataframe (optional)
  #   ...: optional arguments passed on to glm()
  #

  # 1. CHECK FOR VALID ARGUMENTS:
  if (!is.element(group, names(df))) {
    stop("Group column not found. Stopping.")
  }
  if (!is.null(lat)) {
    if (!is.element(lat, names(df))) {
      stop("Latitude column not found. Stopping.")
    }
  }

  if (!is.null(lon)) {
    if (!is.element(lon, names(df))) {
      stop("Longitude column not found. Stopping.")
    }
  }

  if (!is.null(timecol)) {
    if (!is.element(timecol, names(df))) {
      stop("Time column not found. Stopping.")
    }
  }

  message("All parameters OK. Continuing.")

  # 2. CLEAN AND PREP DATA
  # Remove rows without group column
  df <-df[!is.na(df[[group]]),]

  # Sort dataframe by group column
  df <- df[order(df[group]),]

  # 3a. REGRESSIONS (NO TIME COLUMN)
  # if no timecol is specified
  if (is.null(timecol)) {
    message("No time variable selected. Conducting group-level regressions only.")
    # Run regular group-level regression
    slopesByGroup <- groupRegress(df, group, model, ...)
    # rename intercept column
    colnames(slopesByGroup) <- c("intercept", colnames(slopesByGroup)[2:length(colnames(slopesByGroup))])

    # Add other important variables
    if (!is.null(slopesByGroup)) {  # if groupRegress does not return NULL
      groupDf <- unique(df[c(timecol, group, lat, lon)])
      groupSize <- as.matrix(table(df[[group]]))
      fullByGroupDf <- data.frame(cbind(groupDf, size = groupSize, slopesByGroup))
      rownames(fullByGroupDf) <- NULL
      return(fullByGroupDf)
    }
  }

  # 3b. REGRESSIONS (TIME COLUMN SPECIFIED)
  else {
    message("Time variable found. Conducting regressions by time and group.")
    listByTime <- split(df, df[[timecol]])  # split df into list of dfs, by timecol

    # create list of slope dfs
    slopesListByTime <- lapply(listByTime, function(x) {
      message(paste('Conducting regressions for time:', x[1, timecol]))
      tryCatch(groupRegress(x, group, model, ...), error = function(e) NA)
    })

    # add time and other variables to each dataframe
    for (t in names(slopesListByTime)) {
      if (!is.null(slopesListByTime[[t]])) {  # if list NOT NULL
        # create working copy of time-specific slope dataframe
        timeDf <- slopesListByTime[[t]]

        # Add other important variables
        groupDf <- unique(listByTime[[t]][c(timecol, group, lat, lon)])
        groupSize <- as.matrix(table(listByTime[[t]][[group]]))
        fullByGroupDf <- cbind(groupDf, size = groupSize, timeDf)
        rownames(fullByGroupDf) <- NULL
        slopesListByTime[[t]] <- fullByGroupDf
      }
      else {  #if regressions for time t are NULL
        slopesListByTime[t] <- NULL # delete empty list entries
      }
    }

    # combine list of slope dataframes into single dataframe
    slopesDfByTime <- data.frame(do.call(rbind, slopesListByTime))
    rownames(slopesDfByTime) <- NULL
    return(slopesDfByTime)
  }
}

# Performs regression on each subset of observations.
# used in femaRegress
groupRegress <-
function(df, group, model, ...) {
  # Returns a dataframe of slopes.
  #
  # Args:
  #   df: dataframe to be used
  #   group: column in dataframe by which data is nested
  #   model: the actual model for running regression,
  #     for example, 'y ~ x1 + x2'
  #   ...: optional arguments passed on to glm()
  # Run regressions by group
  regressList <- by(df, df[group],
    function(x) {
      tryCatch(glm(as.formula(model), data = x, ...),  # Regression
      error = function(e) NA) # Errors return NA
    }
  )
  # Strip out slopes and constants
  slopesList <- lapply(regressList, function(x) {
    # vector of regression coefficients
    tryCatch(x[[1]], error = function(e) NA)
  }
  )
  # Convert to dataframe
  slopesDf <- do.call(rbind, slopesList)
  # if dataframe is NOT empty:
  if (FALSE %in% is.na(slopesDf)) {
    return(slopesDf)
  }
  # if dataframe is empty, return error and NULL
  else {
    message('Data contains too many missing values.')
    return(NULL)
  }
}


# Calculate field range, effect, and strength using social distance.
femaDo.socDist <- function(df, model) {
  # df: dataframe of individual level data to use.
  # Returns a list with range, effect, and strength elements.
  #

  # prepare test data sor soc dist.
  # df1 <- df[ ,c(sdist.vars, field.i, field.d)]   # remove all but group vars
  print("Preparing data. Generating group variable for social distance.")
  df$group <- interaction(df[ , model$sdist.vars])  # create group parameter
  df <- df[which(!is.na(df$group)),]     # drop ungrouped (incomplete) cases
  df$group <- factor(df$group)           # drop unpopulated group levels
  print("Done!")

  # CALCULATE SOCIAL DISTANCE MATRIX
  print("Calculating social distance.")
  sdist <- socDist(df[,model$sdist.vars], na.cases = "mean")
  # convert to matrix, keep lower tri only (no diag)
  sdist <- as.matrix(sdist)
  sdist[upper.tri(sdist, diag = TRUE)] <- NA
  print("Done!")

  # SLOPE CORRELATION MATRIX
  # group-level regressions
  print("Running group-level regressions.")
  slopes <- femaRegress(df, group = "group", model = model$model, family = model$family)
  # rescale slopes to normal distribution
  slopes[,3:length(colnames(slopes))] <- as.data.frame(scale(slopes[,3:length(colnames(slopes))]))
  print('Done!')

  # calculate slope correlations
  print("Calculating slope correlations.")
  slopeCor <- vCorr(slopes, model$field.var, group = "group")
  print("Done!")

  # FIELD MEASURES

  # create vectors from lower triangle of each matrix
  slopeCor.v <- slopeCor[lower.tri(slopeCor)]
  sdist.v <- sdist[lower.tri(sdist)]

  # logical vector for identifying missing data
  complete <- complete.cases(cbind(slopeCor.v, sdist.v))

  # list of field measures
  fm <- list()

  # field range
  fieldRange <- sum(slopeCor.v[complete] * sdist.v[complete]) / sum(slopeCor.v[complete])
  print(paste("Field Effect Range:", round(fieldRange, 4)))
  fm$range <- fieldRange

  # standardized field effect
  stdEffect <- sum(slopeCor.v[complete] * sdist.v[complete]) / sum(sdist.v[complete])
  print(paste("Standardized Field Effect:", round(stdEffect, 4)))
  fm$stdEffect <- stdEffect

  # field strength
  # NOTE: STILL NEED TO CALCULATE FIELD STRENGTH FOR MULTI SLOPE CORRS.
  if (length(model$field.var) == 1) {
    v <- mean(abs(df[[model$field.var]]), na.rm=T)
    strength <- stdEffect * v
    print(paste("Field Strength:", round(strength, 4)))
    fm$strength <- strength
  }

  # assemble results into list
  results = list(
    sumDist = sum(sdist.v, na.rm = TRUE),
    meanDist = mean(sdist.v, na.rm = TRUE),
    # slopes = slopes,
    # slopeCor = slopeCor,
    range = fieldRange,
    effect = stdEffect,
    strength = strength
    )

  return(results)
}

