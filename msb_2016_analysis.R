#NOTE: MUST RUN MSB_2016_CLEAN FIRST!!!
#NOTE: ASSUMES THAT PLYR IS INSTALLED BUT NOT LOADED
#SETTINGS
proj.old <- '+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0'

#ADDITIONAL LIBRARIES
#the aiR library can be installed using the following:
#devtools::install_github("aslez/femaR")
library(stringr)
library(femaR)
library(spgwr)
library(ggplot2)
library(grid)
library(maptools)
library(rgdal)

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

#FIGURE 8
ggplot(fig8, aes(x = Year, y = Range)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0, 1200)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        panel.grid = element_blank())

