# load packages
library(sf)
library(GPvecchia)
library(ggplot2)
library(dplyr)
library(blockCV)
library(stars)

# load data
data_all = read.csv("./data/me_coords_bio_tcc.csv")
me = st_read("data", "me_cb_2018_us_state_20m_epsg6933")

# partition data into fitting/testing with polygons
data.sf = sf::st_as_sf(data_all, coords = c("x", "y"))
st_crs(data.sf) = st_crs(me)

set.seed(123)
folds = cv_spatial(x = data.sf,
                            k = 5,
                            rows_cols = c(30, 30),
                            selection = "random",
                            iteration = 1,
                            biomod2 = F,
                            plot = T,
                            report = T)

gedi.training = data_all[folds$folds_list[[1]][[1]],] %>% tidyr::drop_na()
save(gedi.training, file = "gedi.training.Rdata")

gedi.holdouts = data_all[folds$folds_list[[1]][[2]],] %>% tidyr::drop_na()
save(gedi.holdouts, file = "gedi.holdouts.Rdata")

gedi.missing = data_all[is.na(data_all$bio),]
save(gedi.missing, file = "gedi.missing.Rdata")


folds
ggsave("plots/blockcv_grid.png", width = 5, height = 4)

data.raster = stars::st_rasterize(data.sf %>% select(bio, geometry))
fold_plot = cv_plot(cv = folds,
        x = data.sf,
        r = data.raster,
        raster_colors = terrain.colors(10, alpha = 0.5),
        label_size = 4,
        num_plots = 1)


sp.auto = cv_spatial_autocor(r = data.raster,
    x = data.sf %>% tidyr::drop_na(),
                   column = 'bio',
                   num_sample = 500)


data.sf$col = "training"
data.sf$col[folds$folds_list[[1]][[2]]] = 'holdout'
data.sf$col[is.na(data.sf$bio)] = 'missing'

data.raster$col = "training"
data.raster$col[folds$folds_list[[1]][[2]]] = 'holdout'
data.raster$col[is.na(data.sf$bio)] = 'missing'

library(wesanderson)
pal = wes_palette("Royal1", 4, type = "discrete")[c(4,3,1)]

ggplot() + 
    geom_sf(data = data.sf, aes(color = col), alpha = 0.8, size = 0.001) +
    geom_sf(data = folds$blocks, fill = NA, color = "black") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          legend.text=element_text(size = 12)) +
    scale_color_manual(name="", values = pal) +
    guides(color = guide_legend(override.aes = list(size=4)))

ggsave("plots/gedi_holdouts.png", width = 7, height = 6, dpi = 600)
