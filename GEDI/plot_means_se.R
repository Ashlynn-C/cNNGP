library(stars)
library(sp)
library(scico)
library(ggplot2)
library(dplyr)
library(ggh4x)

#### Plot means and standard error #############################################
load("data/NNGP_fitted_means_se.Rdata")
colnames(nngp.means.se)[colnames(nngp.means.se) == 'post.mean'] = 'NNGP.mean'
colnames(nngp.means.se)[colnames(nngp.means.se) == 'post.sd'] = 'NNGP.se'

load("data/cNNGP_fitted_means_se.Rdata")
colnames(cnngp.means.se)[colnames(cnngp.means.se) == 'post.mean'] = 'cNNGP.mean'
colnames(cnngp.means.se)[colnames(cnngp.means.se) == 'post.sd'] = 'cNNGP.se'

gedi.fitting = dplyr::left_join(cnngp.means.se,
                                nngp.means.se,
                                by = c("x","y")) %>%
  mutate(x = x*1000, y = y*1000)

load("data/missing_holdout_predictions.Rdata")

gedi.data = rbind(gedi.fitting,gedi.predictions)

means.df = gedi.data %>% tidyr::pivot_longer(cols = contains("mean"),
                                          names_to = "Model",
                                          values_to = "mean") %>%
                            mutate(Model = gsub(".mean", "", Model))

me = st_read("data", "me_cb_2018_us_state_20m_epsg6933")
means.sf = st_as_sf(means.df, coords = c('x', 'y'))
st_crs(means.sf) = st_crs(me)

bio.col = colorRampPalette(scico(100, palette = 'roma'), bias = 3)

ggplot(means.df) +
  geom_sf(data = means.sf,fill=NA,alpha=0) +
  geom_raster(data = means.df, aes(x = x, y = y, fill = mean), interpolate = T) +
  theme_bw() +
  facet_grid(~Model) +
  theme_minimal() +
  theme(strip.text.x = element_text(colour = 'black', size = 14),
        axis.text = element_text(size = 10, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.background = element_blank(),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12)) +
  xlab("") + ylab("") +
  scale_fill_gradientn("Predicted \nBiomass \n(Mg/ha)",
                        colors = bio.col(100), na.value = NA)

ggsave("plots/gedi-means.png", dpi = 600, width = 9, height = 5)

se.df = gedi.data %>% tidyr::pivot_longer(cols = contains("se"),
                                             names_to = "Model",
                                             values_to = "mean") %>%
  mutate(Model = gsub(".se", "", Model))

se.sf = st_as_sf(se.df, coords = c('x', 'y'))
st_crs(se.sf) = st_crs(me)

ggplot() +
  geom_sf(data = se.sf,fill=NA,alpha=0) +
  geom_raster(data = se.df, aes(x = x, y = y, fill = mean), interpolate = T) +
  theme_bw() +
  facet_grid(~Model) +
  theme_minimal() +
  theme(strip.text.x = element_text(colour = 'black', size = 14),
        axis.text = element_text(size = 10, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.background = element_blank(),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12)) +
  xlab("") + ylab("") +
  scale_fill_gradientn("Standard \nerror",
                       colors = bio.col(100), na.value = NA)

ggsave("plots/gedi-se.png", dpi = 600, width = 9, height = 5)

################################################################################
### Plot holdout prediction scatterplots #####
################################################################################
load("data/gedi.holdouts.Rdata")
load("data/gedi.training.Rdata")

library(scales)

pal = pal_viridis(option = "G")(20)
show_col(pal_viridis(option = "G")(20))

plt.data = gedi.holdouts %>% left_join(gedi.predictions, by = c("x", "y")) %>%
  select(!contains("se")) %>%
  tidyr::pivot_longer(cols = contains("mean"), names_to = "Model", values_to = "post.mean") %>%
  mutate(Model = gsub(".mean", "", Model)) %>%
  mutate(dist = abs(post.mean-bio)+1)

plt.data %>% group_by(Model) %>% summarize(cor = cor(post.mean, bio))

ggplot(plt.data, aes(y = bio, x = post.mean, color = dist)) +
  geom_point(alpha = 0.2, size = 0.7, pch = 19) +
  facet_nested(~Model) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text.x = element_text(colour = 'black', size = 13),
        #strip.text.y = element_text(colour = 'black', size = 12),
        axis.text = element_text(size = 8, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        #panel.grid = element_blank(),
        panel.background = element_blank()) +
  geom_abline(slope=1, intercept=0,lty = "dashed",lwd = 0.8) +
  xlim(c(0,355)) + ylim(c(0,355)) +
  scale_color_gradientn(
    colours = c(pal[20], pal[10], pal[6]),
    transform = "sqrt"
    ) +
  ylab("Observed biomass (Mg/ha)") + xlab("Predicted biomass (Mg/ha)")

ggsave("gedi holdout predictions.png", dpi = 300, width = 7, height = 3)

################################################################################
### Plot holdout prediction differences #####
################################################################################

diff.df = gedi.data %>%
  mutate(diff = cNNGP.mean - NNGP.mean) %>%
  filter(diff >= quantile(diff,probs = 0.01) & diff <= quantile(diff, probs = 0.99))

diff.sf = st_as_sf(diff.df, coords = c('x', 'y'))
st_crs(diff.sf) = st_crs(me)

ggplot() +
  geom_sf(data = diff.sf,fill=NA,alpha=0) +
  geom_raster(data = diff.df,aes(x = x, y = y, fill = diff), interpolate = T) +
  theme_bw() +
  #facet_grid(~Model) +
  theme_minimal() +
  theme(strip.text.x = element_text(colour = 'black', size = 14),
        axis.text = element_text(size = 9, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.background = element_blank(),
        legend.text=element_text(size=9),
        legend.title=element_text(size=9)) +
  xlab("") + ylab("") +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name = "cNNGP - NNGP\nPredicted\nBiomass")

ggsave("plots/gedi model mean diff.png", dpi = 300, width = 7, height = 3)

se.diff.df = gedi.data %>%
  mutate(diff = cNNGP.se - NNGP.se) %>%
  filter(diff >= quantile(diff,probs = 0.01) & diff <= quantile(diff, probs = 0.99))

se.diff.sf = st_as_sf(se.diff.df, coords = c('x', 'y'))
st_crs(se.diff.sf) = st_crs(me)

ggplot() +
  geom_sf(data = se.diff.sf,fill=NA,alpha=0) +
  geom_raster(data = se.diff.df,aes(x = x, y = y, fill = diff), interpolate = T) +
  theme_bw() +
  #facet_grid(~Model) +
  theme_minimal() +
  theme(strip.text.x = element_text(colour = 'black', size = 14),
        axis.text = element_text(size = 9, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.background = element_blank(),
        legend.text=element_text(size=9),
        legend.title=element_text(size=9)) +
  xlab("") + ylab("") +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name ="cNNGP - NNGP\nStandard\nError")

ggsave("plots/gedi model se diff.png", dpi = 300, width = 9, height = 5)

###############################################################################
# Plot difference to true biomass
###############################################################################
load("data/gedi.missing.Rdata")

mean.true.df = gedi.data %>%
      left_join(gedi.holdouts, by = c("x", "y")) %>%
      left_join(gedi.training, by = c("x", "y")) %>%
      anti_join(gedi.missing, by = c("x", "y")) %>%
      mutate(bio = ifelse(!is.na(bio.x),bio.x,bio.y)) %>%
      select(x,y,contains("mean"), bio) %>%
      tidyr::pivot_longer(cols = contains("mean"),
                          names_to = "Model", values_to = "post.mean") %>%
      mutate(Model = gsub(".mean", "", Model),
             diff = post.mean - bio) %>%
      filter(diff >= quantile(diff,probs = 0.01) & diff <= quantile(diff, probs = 0.99))

ggplot() +
  geom_sf(data = se.diff.sf,fill=NA,alpha=0) +
  geom_raster(data = mean.true.df,aes(x = x, y = y, fill = diff), interpolate = T) +
  theme_bw() +
  facet_grid(~Model) +
  theme_minimal() +
  theme(strip.text.x = element_text(colour = 'black', size = 14),
        axis.text = element_text(size = 10, color = "black"),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.background = element_blank(),
        #panel.spacing.x = unit(1, "lines"),
        legend.text=element_text(size=11),
        legend.title=element_text(size=12)) +
  xlab("") + ylab("") +
  scale_fill_gradient2(
    low = "#4575B4", mid = "white", high = "#D73027",
    midpoint = 0, name ="Estimated - true\nbiomass\n(Mg/ha)")

ggsave("plots/gedi estimated-true.png", dpi = 300, width = 9, height = 5)
