#SETTINGS
#set working directory to replication folder

#LIBRARIES
library(dplyr)
library(tidyr)

#######################
# CLEAN ECOLOGICAL DATA
#######################
#LOAD DATA
load("msb_2016.Rdata")

#BUILD TABULAR DATA
#project data on to intersections using areal weighting
int_df <- left_join(int_dat@data, clst) %>%
  left_join(census, by = c("FIPS1890" = "FIPS")) %>%
  gather(YEAR, FIPS, -(AREA:WHEAT)) %>%
  mutate(YEAR = as.numeric(substring(YEAR, 5, 8))) %>%
  left_join(voting) %>%
  group_by(FIPS, YEAR) %>%
  mutate(WGT = AREA / sum(AREA),
         GOVV = GOVV * WGT, 
         GOVT = GOVT * WGT,
         WHEAT = WHEAT * WGT,
         WHEATAC = WHEATAC * WGT) %>%
  group_by(comp) %>%
  filter(all(GOVT > 0) & all(WHEATAC > 0)) %>%
  ungroup() %>%
  select(FIPS, YEAR, STATE,
         EDGE = edge, COMP = comp, 
         GOVV, GOVT, WHEAT, WHEATAC)

#aggregate up to real counties
real_df <- int_df %>%
  group_by(FIPS, YEAR, STATE, COMP) %>%
  summarise_each(funs(sum)) %>%
  ungroup() %>%
  mutate(WHO = WHEAT / WHEATAC,
         VOTE = (GOVV / GOVT) * 100) %>%
  select(FIPS, YEAR, WHO, VOTE)

#update intersection data
int_df <- int_df %>%
  mutate(WHO = WHEAT / WHEATAC,
         VOTE = (GOVV / GOVT) * 100) %>%
  select(EDGE, FIPS, YEAR, WHO, VOTE) %>%
  left_join(clst %>% select(-comp), by = c("EDGE" = "edge"))

#JOIN TABULAR DATA TO SHAPEFILES
int_pdf_list <- lapply(seq(1890, 1896, 2), function(x) {
  shp <- int_dat
  shp$ROWNUM <- 1:NROW(shp@data)
  d_old <- shp@data %>% select(contains("FIPS"), ROWNUM)
  d <- left_join(d_old, subset(int_df, YEAR == x)) %>%
    arrange(ROWNUM)
  rownames(d) <- d$ROWNUM
  shp@data <- d
  shp[!is.na(shp$VOTE), ]
})

real_pdf_list <- lapply(seq(1890, 1896, 2), function(x) {
  shp <- shplist[[match(x, seq(1890, 1896, 2))]]
  shp$ROWNUM <- 1:NROW(shp@data)
  d_old <- shp@data %>% select(FIPS, ROWNUM)
  d <- left_join(d_old, subset(real_df, YEAR == x)) %>%
    arrange(ROWNUM)
  rownames(d) <- d$ROWNUM
  shp@data <- d
  shp[!is.na(shp$VOTE), ]
})
