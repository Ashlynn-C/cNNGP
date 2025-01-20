library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

data_GEDI <- read_csv("data/me_coords_bio_tcc.csv")

data_long = data_GEDI %>% filter(!is.na(bio)) %>%
  pivot_longer(cols = bio:tcc, names_to = "var")

ggplot(data_long, aes(x = value, fill = var)) +
  geom_density(alpha = 0.7) +
  facet_grid(~var, scales = "free",
             labeller = as_labeller(c("bio" = "Biomass (Mg/ha)", "tcc" = "TCC (%)"))) +
  theme_minimal() +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        strip.text = element_text(colour = 'black', size = 11)) +
  scale_fill_manual(values = c("bio" = "#4575b4", "tcc" = "#66bd63")) +
  xlab("Value") + ylab("Density")
ggsave("bio_tcc_dist.png", width = 5, height = 2.5)

ggplot(data_GEDI, aes(x = tcc, y = bio, color = bio)) +
  geom_point(pch=19,alpha=0.35,size=0.8) +
  theme_minimal() +
  scale_color_viridis_c(option = "G",direction = 1) +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        strip.text = element_text(colour = 'black', size = 11)) +
  xlab("TCC (%)") + ylab("Biomass (Mg/ha)")
ggsave("bio_vs_tcc.png", width = 5, height = 3)

