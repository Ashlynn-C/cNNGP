# code from Andy Finley
rm(list=ls())
library(tidyverse)
library(viridis)
library(sf)
library(stars)
library(patchwork)
library(scico)

## ~1x1 km GEDI L4B version 2 https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2299
## 2019-04-18 to mission week 223 ending on 2023-03-16,
## mean aboveground biomass density (AGBD) (Mg/ha)

bio <- read_stars("data/me_GEDI04_B_MW019MW223_02_002_02_R01000M_MU.tif")

## Remove some crazy high biomass values that are just not believable for ME.
bio[bio > 350] <- NA

## ~1x1 m tree canopy cover (%),
## USFS Science TCC 2021 https://data.fs.usda.gov/geodata/rastergateway/treecanopycover/
## I resampled this from 30x30 m to ~1x2 km via
## gdalwarp -ts 397 396 -r average me_science_tcc_conus_2020_v2021-4_epsg6933.tif me_science_tcc_conus_2020_v2021-4_epsg6933_1km.tif

tcc <- read_stars("data/me_science_tcc_conus_2020_v2021-4_epsg6933_1km.tif")
tcc[tcc > 100] <- NA ##Values should only be between 0 and 100, all else are NA.

## Shapefile for me.
me <- st_read("data", "me_cb_2018_us_state_20m_epsg6933")

st_crs(tcc) <- st_crs(me)
st_crs(bio) <- st_crs(me)

## Take a look.

## There are a few large biomass values that make it hard to see the variation,
## hence the bias argument in the color ramp.
bio.col <- colorRampPalette(scico(100, palette = 'roma'), bias = 3)

p.bio <- ggplot() +
    geom_stars(data = bio) +
    geom_sf(data = me, fill = NA, color = "orange") +
    scale_fill_gradientn("Biomass (Mg/ha)", colors = bio.col(100), na.value = NA)

p.tcc <- ggplot() +
  geom_stars(data = tcc) +
    geom_sf(data = me, fill = NA, color = "orange") +
    scale_fill_gradientn("TCC (%)", colors = scico(100, palette = 'bamako',  direction = -1), na.value = NA)

p.bio + p.tcc

p.bio + ylab("") + xlab("") + theme_bw() +
theme(axis.text = element_text(size = 12),
      legend.text=element_text(size = 12))
ggsave("plots/bio.png", width = 7, height = 6, dpi = 600)

p.tcc + ylab("") + xlab("") + theme_bw() +
    theme(axis.text = element_text(size = 12),
          legend.text=element_text(size = 12))
ggsave("plots/tcc.png", width = 7, height = 6, dpi = 600)

## Get TCC at each gedi pixel.
dat <- bio
dat.df <- as.data.frame(dat, xy = TRUE)
colnames(dat.df)[3] <- "bio"
dat.df$indx <- 1:nrow(dat.df)
dat.sf <- st_as_sf(dat.df, coords = c("x","y"), crs = st_crs(me))

## Grab only those pixels in me.
dat.sf <- st_intersection(dat.sf, me) %>% select(bio, indx)
head(dat.sf)
dat.sf$tcc <- st_extract(tcc, dat.sf) %>% st_drop_geometry()

dat.df <- as.matrix(cbind(st_coordinates(dat.sf), st_drop_geometry(dat.sf)))

colnames(dat.df) <- c("x","y","bio","bio_pixel_indx","tcc")

dat.df <- as.data.frame(dat.df[,c("x","y","bio_pixel_indx","bio","tcc")])

head(dat.df)

plot(dat.df[,c("tcc","bio")])

plot(dat.df[,c("x","y")])

##Just remove rows with missing TCC.
dat.df <- dat.df[!is.na(dat.df$tcc),]

write_csv(dat.df, file = "me_coords_bio_tcc.csv")
