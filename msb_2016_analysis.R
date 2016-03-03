# Script assumes that all files are at the top level of the
#   working directory.
# setwd("SET/DIRECTORY/HERE/IF/NECESSARY")

# MUST RUN MSB_2016_CLEAN FIRST!!!
source("msb_2016_clean.R")

# NOTE: ASSUMES THAT PLYR IS INSTALLED BUT NOT LOADED

#ADDITIONAL LIBRARIES
library(stringr)
library(spgwr)
library(ggplot2)
library(grid)
library(maptools)
library(rgdal)

# load custom functions
source("functions.R")

# SETTINGS
proj.old <-
  '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'


#########################################################################
#########################################################################
#
## ---- figure 2 ----------------------------------------------
#
#########################################################################
#########################################################################

poly_dat <- rbind(c(1, 0), fig2, c(101, 0))
ggplot() +
  geom_line(data = fig2, aes(x = x, y = y), size = 1) +
  geom_polygon(data = poly_dat, aes(x = x, y = y), alpha = .2) +
  scale_x_continuous("Distance", limits = c(1, 101), breaks = seq(1, 101, 4)) +
  scale_y_continuous("Degree of alignment", limits = c(0, 100)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())


#########################################################################
#########################################################################
#
## ---- FIGURE 3 ----------------------------------------------
#
#########################################################################
#########################################################################

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

#########################################################################
#########################################################################
#
## ---- FIGURE 4 ----------------------------------------------
#
#########################################################################
#########################################################################

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

#########################################################################
#########################################################################
#
## ---- FIGURE 5 ----------------------------------------------
#
#########################################################################
#########################################################################

ggplot(int_lst$measures, aes(x = Year, y = Value)) +
  geom_line() +
  facet_wrap(~ measure, scale = 'free_y', ncol = 1) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())

#########################################################################
#########################################################################
#
## ---- FIGURE 6 ----------------------------------------------
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

# PLOT FIGURE 6
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
## ---- FIGURE 7  ----------------------------------------------
#
#########################################################################
#########################################################################

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
## ---- FIGURE 8 ----------------------------------------------
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
## ---- FIGURE 9 ----------------------------------------------
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

#### CALCULATE MEASURES

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


