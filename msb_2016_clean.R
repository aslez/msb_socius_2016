#SETTINGS

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

#######################
# CLEAN CCES 2008 DATA
#######################

# CREATE NEW VARIABLES

# generate age from birth year
cces$age <- 2008 - as.integer(as.character(cces$V207))

# REMOVE UNWANTED VARIABLES

# var list required for FEMA analysis
vars <- c("fips_dist" = "fips_dist", "lat" = "lat", "lon" = "lon",
  "race" = "V211", "edu" = "V213", "inc" = "V246", "age" = "age",
  "party" = "CC424"
)

# keep only relevant variables, and rename
cces <- cces[, vars]
names(cces) <- names(vars)


### RECODE AND CONVERT FROM FACTORS TO NUMERIC

# replace unused reponses with NA
for (i in colnames(cces)) {
  nas <- which(cces[[i]] %in% c("Not sure", "Skipped", "Not Asked",
    "Missing", "Prefer not to say", "Don't know"))
  print (paste(i, ":", length(nas)))
  cces[nas, i] <- NA
}


# recode all ordinals to numeric
for (i in names(cces)) {
  cces[[i]] <- as.numeric(cces[[i]])
}


# 7-point party id
cces[["party"]][cces["party"] > 7] <- NA

# race (0 = white, 1 = black)
cces$race <- replace(cces$race, cces$race > 2, NA)
cces$race <- cces$race - 1

# income (remove missing/refused)
cces$inc <- replace(cces$inc, cces$inc > 14, NA)

# CLEANUP
rm(nas, i)


#######################
# CLEAN GSS DATA
#######################

# AGE CATEGORICAL, broken up into quintiles
gss$agecat <- cut(gss$age,
  breaks = quantile(gss$age, probs = seq(0, 1, .2), na.rm = T),
  include.lowest = TRUE
  )

# PRESTIGE
# combine pre and post 1980 scores
gss$prestige <- gss$prestige
gss$prest <- gss$prestg80
gss$prest[is.na(gss$prest)] <- gss$prestige[is.na(gss$prest)]


# PRESIDENTIAL VOTE
# create single variable for presidential vote (by party)
gss$pvote <- NA
for (y in rev(c("68", "72", "76", "80", "84", "88", "92", "96",
  "00", "04", "08"))) {

    # variable names
    sy <- as.character(y)
    pv <- paste("pres", sy, sep = "")      # pres vote
    iv <- paste("if", sy, "who", sep = "") # 'if' vote

    # vector of votes for given year
    vote <- as.integer(gss[, pv])   # vector of pres votes
    vote[is.na(vote)] <- as.integer(gss[is.na(vote), iv]) # add 'if' votes

    # add to df
    gss[is.na(gss$pvote), "pvote"] <- vote[is.na(gss$pvote)]
}

# remove all but two major parties
gss[which(gss$pvote > 3), "pvote"] <- NA
# make binary: 0 = dem, 1 = rep
gss$pvote <- gss$pvote - 2


# remove unused factor levels for all factor cols
for (col in colnames(gss)) {
  if (is.factor(gss[[col]])) {
    gss[[col]] <- factor(gss[[col]])
  }
}

# cleanup
rm(col)

