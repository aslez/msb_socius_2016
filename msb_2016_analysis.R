# NOTE: MUST RUN MSB_2016_CLEAN FIRST!!!
setwd("~/Projects/msb_socius_2016")

source('msb_2016_clean.R')

#NOTE: ASSUMES THAT PLYR IS INSTALLED BUT NOT LOADED

#SETTINGS
proj.old <- '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'

#ADDITIONAL LIBRARIES
#the femaR library can be installed using the following:
# devtools::install_github("aslez/femaR")
library(stringr)
library(femaR)
library(spgwr)
library(ggplot2)
library(grid)
library(maptools)
library(rgdal)

# additional functions
source("functions.r")

#DEFINE FUNCTIONS
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

#FIGURE 2
poly_dat <- rbind(c(1, 0), fig2, c(101, 0))
ggplot() +
  geom_line(data = fig2, aes(x = x, y = y), size = 1) +
  geom_polygon(data = poly_dat, aes(x = x, y = y), alpha = .2) +
  scale_x_continuous("Distance", limits = c(1, 101), breaks = seq(1, 101, 4)) +
  scale_y_continuous("Degree of alignment", limits = c(0, 100)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())

#FIGURE 3
#generate and reproject shapefiles
int_gwr_prj <- lapply(int_pdf_list,
                      function(x) spTransform(x, CRS(proj.old)))
real3_gwr_prj <- lapply(real_pdf_list,
                        function(x) spTransform(x, CRS(proj.old)))
st_shp <- unionSpatialPolygons(shplist[[1]], shplist[[1]]$STATE_TERR)
st_shp_prj <- spTransform(st_shp, CRS(proj.old))

#analysis
res_lst <- fema(real3_gwr_prj)
int_lst <- fema(int_gwr_prj)

#vector maps
ggplot(data = res_lst$shapes, aes(long, lat, group = group)) +
  geom_polygon(data = fortify(st_shp_prj), fill = 'grey90') +
  geom_polygon(fill = 'white', colour = 'grey50') +
  geom_segment(data = res_lst$centroids,
               aes(x = xs, xend = xe, y = ys, yend = ye),
               arrow = arrow(angle = 30,
                             length = unit(0.1, 'cm'),
                             type = 'closed')) +
  geom_path(data = fortify(st_shp_prj), size = 1) +
  facet_wrap(~ YEAR, ncol = 2) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.key = element_rect(colour = 'black'))

#FIGURE 4
line_df <- data.frame(x = c(0, 0),
                      xend = c(1, 1),
                      y = c(0, 0),
                      yend = c(.75, 2.5),
                      vec = 1:2)
ggplot() +
  geom_segment(data = line_df,
               aes(x = x, xend = xend,
                   y = y, yend = yend,
                   colour = factor(vec)),
               arrow = arrow(),
               size = 2) +
  scale_x_continuous('Independent variable', limits = c(0, 3)) +
  scale_y_continuous('Dependent variable', limits = c(0, 3)) +
  scale_colour_grey(guide = FALSE) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())

#FIGURE 5
ggplot(int_lst$measures, aes(x = Year, y = Value)) +
  geom_line() +
  facet_wrap(~ measure, scale = 'free_y', ncol = 1) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())


#########################################################################
#########################################################################
#
## ---- figure 6 ----------------------------------------------
#
# Code to produce figure 6. Assumes CCES08 data is loaded/cleaned
#
# Requires :
# - femaRegress
# - geoDistMatrix
# - distCorr
#
#########################################################################
#########################################################################

# generate district-level slopes
m1 <- femaRegress(
  cces,
  group = "fips_dist",
  model = "party ~ edu + inc + race",
  lon = "lon", lat = "lat",
  family = gaussian
)

# FIGURE 7 -- ALTERNATE VERSION (LINE ONLY)
geo <- geoDistMatrix(m1, group = "fips_dist", units = "mi")
corr <- distCorr(
  m1,
  group = "fips_dist",
  slopes = c("inc", "edu", "race"),
  zslopes = T
)

# create vectors from lower triangle of each matrix
corrv <- corr[lower.tri(corr)]
geov <- geo[lower.tri(geo)]

# logical vector for identifying missing data
complete <- complete.cases(cbind(corrv, geov))

# plot points and fitted line
smSpline <- smooth.spline(geov[complete], corrv[complete], spar = .75)

# convert to dataframe
lineDf <- data.frame(x = smSpline$x, y = smSpline$y)

###############
# PLOT FIGURE 6
###############

ggplot(lineDf, aes(x = x, y = y)) +
    geom_line(size = 1.25) +
    labs(
      x = "Geographical distance (miles)",
      y = "Correlation"
    ) +
    theme_bw() +
    theme(
      aspect.ratio = 2 / 3,
      text = element_text(size = 20),
      axis.text = element_text(size = 13),
      panel.grid = element_blank()
    )


#########################################################################
#########################################################################
#
## ---- figure 7  ----------------------------------------------
#
#########################################################################
#########################################################################

# DEFINE CUSTOM FUNCTIONS

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


# CALCULATE GEOGRAPHICAL DISTANCE
geo <- geoDistMatrix(acsSoc$data)

# CALCULATE EUCLIDEAN DISTANCE
highschool <- eDist(acsSoc$data$V264)
oneunit <- eDist(acsHouse$data$V30 + acsHouse$data$V34)
poverty <- eDist(acsEcon$data$V444)
inc <- eDist(acsEcon$data$V252)

# set plot parameters
par(
  mfrow = c(2, 2),
  mar = c(4.2, 4, 1.5, 2),
  mgp = c(2.1, 1, 0)
)

# set label size
labSize <- 1.15
# set bin size
bin <- 1000

# plots
binPlot(geo, oneunit, bin, "Percent One-Family Structures", cex.lab = labSize)
binPlot(geo, inc, bin, "Median Household Income", cex.lab = labSize)
binPlot(geo, poverty, bin, "Poverty Rate", cex.lab = labSize)
binPlot(geo, highschool, bin, "High School Graduation Rate", cex.lab = labSize)


#########################################################################
#########################################################################
#
## ---- figure 8 ----------------------------------------------
#
#########################################################################
#########################################################################

ggplot(fig8, aes(x = Year, y = Range)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0, 1200)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())

#########################################################################
#########################################################################
#
## ---- figure 9 ----------------------------------------------
#
#########################################################################
#########################################################################

#### DEFINE MODEL
models <- list(
  m1 = list(
    sdist.vars = c("agecat", "sex", "srcbelt"),
    field.var = "prest",
    model = as.formula("pvote ~ prest"),
    family = "binomial"
  )
)

#### COMPUTE MEASURES

# DEFINE PERIODS
periods <- list(
  per1 = seq(1981, 1984, 1),
  per2 = seq(1985, 1988, 1),
  per3 = seq(1989, 1992, 1),
  per4 = seq(1993, 1996, 1),
  per5 = seq(1997, 2000, 1),
  per6 = seq(2001, 2004, 1),
  per7 = seq(2005, 2008, 1),
  per8 = seq(2009, 2012, 1)
)

# empty list to add model(s)
fe8 <- list(
  )

for (model in names(models)) {
  # temp list for models
  feTemp <- list(
    range = c(),
    effect = c(),
    strength = c()
    )

  # for all periods
  for (per in  periods) {

    # logical vector of range
    perLog <- gss$year %in% per

    fe <- femaDo.socDist(gss[perLog, ], models[[model]])
    for (el in names(fe)) {
      feTemp[[el]] <- append(feTemp[[el]], fe[[el]])
    }
  }
  fe8[[model]] <- feTemp
}


# x-axis labels
period <- c("1981-\n1984", "1985-\n1988", "1989-\n1992", "1993-\n1996",
 "1997-\n2000", "2001-\n2004", "2005-\n2008", "2009-\n2012")

# convert lists to dfs
fe8 <- lapply(fe8, function(x) do.call(cbind.data.frame, x))

fe8$m1$period <- period

# plot
ggplot(fe8$m1, aes(x = period, y = effect, group = NA)) +
  geom_line(size = 1.25) +
  labs(
    x = "Period",
    y = "Measure"
  ) +
  theme_bw() +
  theme(
  aspect.ratio = 2 / 3,
  text = element_text(size = 20),
  axis.text = element_text(size = 13),
  panel.grid = element_blank()
  )


