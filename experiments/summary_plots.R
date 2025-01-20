library(tidyr)
library(dplyr)
library(ggplot2)
library(ggh4x)


load("sim.results.Rdata")

reps = 30
true.sigma.sq = 1
true.tau.sq = 0.1
true.beta.0 = 1
true.beta.1 = 5
n.total = c("n == 2500","n == 10000") 
################################################################################
# parameter estimates 
################################################################################
true.param.vals = data.frame(var = factor(rep(c("beta[0]","beta[1]","tau^2","sigma^2","phi"), times = 3),
                             levels = c("beta[0]", "beta[1]", "tau^2", "sigma^2", "phi")),
                             phi = rep(c("phi == 2.88","phi == 11.51"),each = 5),
                             n = rep(n.total,each = 15),
                             value = c(1,5,0.1,1,2.88,1,5,0.1,1,11.51))

plot.data = sim.results %>% 
  dplyr::select((c("n", "phi", "rep", "id", "model") | 
                             ends_with(c("sigma.sq", "phi", "tau.sq", "beta0", "beta1")) &
                             !contains("ESS"))) %>% 
  pivot_longer(cols = !c("n", "phi", "rep", "id", "model"), names_to = "var") %>% 
  mutate(var = case_when(var == "sigma.sq" ~ "sigma^2",
                         var == "est.phi" ~ "phi",
                         var == "beta0" ~ "beta[0]",
                         var == "beta1" ~ "beta[1]",
                         var == "tau.sq" ~ "tau^2",
                         .default = var)) %>%
  mutate(phi = factor(paste0("phi == ", round(phi,2)),
                              levels = c("phi == 11.51", 
                                         "phi == 2.88")),
         n = factor(paste0("n == ", n), levels = c("n == 10000", "n == 2500"))) %>% 
  mutate(model = gsub("_10", "(10)",model) %>% 
           gsub("_20", "(20)",.) %>% 
           gsub("_30", "(30)",.)) %>% 
  mutate(model = factor(model, levels = c("meshedGP", "blockNNGP", "NNGP(30)", "NNGP(20)",  "NNGP(10)",
                                          "cNNGP(30)", "cNNGP(20)", "cNNGP(10)"))) %>% 
  mutate(var = factor(var, levels = c("beta[0]", "beta[1]", "tau^2", "sigma^2", "phi")))

ggplot(plot.data, aes(x=value, y=model, fill = model)) + 
    geom_boxplot() +
    geom_vline(aes(xintercept = value), data = true.param.vals, 
               color = gray(0.2,0.6), linewidth = 1) +
    facet_nested(n+phi~var, scales = "free_x", 
               labeller = label_parsed, as.table = F,
               nest_line = element_line(color = "gray", linewidth = 1)) +
    theme_minimal() +
    theme(legend.position = 'none',
          #strip.background =element_rect(fill="#D5D5D5"),
          strip.text.x = element_text(colour = 'black', size = 14),
          strip.text.y = element_text(colour = 'black', size = 12),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 40),
          axis.line.x = element_blank(),
          axis.line.y.right = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
          panel.grid = element_blank(),
          panel.background = element_blank()) +
    scale_fill_viridis_d(name="",option = "G",direction = 1) +
  ylab("") + xlab("Parameter Values")

ggsave("plots/param_estimates.png", width = 10, height = 6)

# credible interval widths -----------------------------------------------------
# parameter interval widths ----
param.widths = sim.results %>% 
  mutate(phi.width = phi.ub - phi.lb, 
         sigma.sq.width = sigma.sq.ub - sigma.sq.lb,
         tau.sq.width = tau.sq.ub - tau.sq.lb,
         beta0.width = beta0.ub - beta0.lb,
         beta1.width = beta1.ub - beta1.lb) %>% 
  select((c("n", "phi", "rep", "id", "model") | 
            ends_with("width"))) %>% 
  pivot_longer(cols = contains("width"), names_to = "var") %>% 
  mutate(var = gsub(".width", "", var)) %>% 
  mutate(var = case_when(var == "sigma.sq" ~ "sigma^2",
                         var == "est.phi" ~ "phi",
                         var == "beta0" ~ "beta[0]",
                         var == "beta1" ~ "beta[1]",
                         var == "tau.sq" ~ "tau^2",
                         .default = var)) %>%
  mutate(phi = factor(paste0("phi == ", round(phi,2)),
                      levels = c("phi == 11.51", 
                                 "phi == 2.88")),
         n = factor(paste0("n == ", n), levels = c("n == 10000", "n == 2500"))) %>% 
  mutate(model = gsub("_10", "(10)",model) %>% 
           gsub("_20", "(20)",.) %>% 
           gsub("_30", "(30)",.)) %>% 
  mutate(model = factor(model, levels = c("meshedGP", "blockNNGP", "NNGP(30)", "NNGP(20)",  "NNGP(10)",
                                          "cNNGP(30)", "cNNGP(20)", "cNNGP(10)"))) %>% 
  mutate(var = factor(var, levels = c("beta[0]", "beta[1]", "tau^2", "sigma^2", "phi"))) %>% 
  ungroup() 

# max widths ----
max.widths = param.widths %>% 
  group_by(model, var, n, phi) %>%
  summarize(max = max(value)) %>% 
  ungroup() 

param.widths %>% 
  filter(var == "phi") %>% 
  group_by(model, n, phi) %>%
  summarize(med = median(value)) %>% 
  ungroup() %>% 
  group_by(n,phi) %>% 
  filter(med == max(med))

# coverage counts ----
cov.counts = as_tibble(sim.results) %>% 
  mutate(phi.cov = (phi >= phi.lb & phi <= phi.ub),
         sigma.sq.cov = (true.sigma.sq >= sigma.sq.lb & true.sigma.sq <= sigma.sq.ub),
         tau.sq.cov = (true.tau.sq >= tau.sq.lb & true.tau.sq <= tau.sq.ub),
         beta0.cov = (true.beta.0 >= beta0.lb & true.beta.0 <= beta0.ub),
         beta1.cov = (true.beta.1 >= beta1.lb & true.beta.1 <= beta1.ub)) %>% 
  select((c("n", "phi", "rep", "id", "model") | 
            contains(c("cov","width")) & !c("cov","med.widths"))) %>%
  pivot_longer(cols = contains("cov"), names_to = "var", values_to = "cov") %>% 
  group_by(phi,n,model,var) %>%
  summarize(cov.count = sprintf("%#.2f",sum(cov)/reps)) %>% 
  ungroup() %>%
  mutate(var = gsub(".cov", "", var)) %>% 
  mutate(var = case_when(var == "sigma.sq" ~ "sigma^2",
                         var == "est.phi" ~ "phi",
                         var == "beta0" ~ "beta[0]",
                         var == "beta1" ~ "beta[1]",
                         var == "tau.sq" ~ "tau^2",
                         .default = var)) %>%
  mutate(phi = factor(paste0("phi == ", round(phi,2)),
                      levels = c("phi == 11.51",
                                 "phi == 2.88")),
         n = factor(paste0("n == ", n), levels = c("n == 10000", "n == 2500"))) %>%
  mutate(model = gsub("_10", "(10)",model) %>%
           gsub("_20", "(20)",.) %>%
           gsub("_30", "(30)",.)) %>%
  mutate(model = factor(model, levels = c("meshedGP", "blockNNGP", "NNGP(30)", "NNGP(20)",  "NNGP(10)",
                                          "cNNGP(30)", "cNNGP(20)", "cNNGP(10)"))) %>% 
  mutate(var = factor(var, levels = c("beta[0]", "beta[1]", "tau^2", "sigma^2", "phi"))) %>% 
  right_join(max.widths, by = c("model","var","n","phi")) 

ggplot(data = param.widths, aes(x=value, y=model, fill = model)) + 
  geom_boxplot() +
  geom_text(data = cov.counts,
            aes(x = max, y = model, label = cov.count), size = 3, hjust = -0.4) +
  geom_point(data = cov.counts, aes(x = max*1.5, y = model), alpha = 0) +
  facet_nested(n+phi~var, scales = "free_x", 
               labeller = label_parsed, as.table = F,
               nest_line = element_line(color = "gray", linewidth = 1)) +
  theme_minimal() +
  theme(legend.position = 'none',
        #strip.background =element_rect(fill="#D5D5D5"),
        strip.text.x = element_text(colour = 'black', size = 14),
        strip.text.y = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, angle = 40),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_fill_viridis_d(name="",option = "G",direction = 1) +
  ylab("") + xlab("Credible Interval Widths")

ggsave("plots/param_widths.png", width = 8, height = 7)


ggplot(aes(x=cov.count, y=model, fill = model, label = cov.count)) + 
  geom_col() +
  geom_text(position = position_dodge(width = 1), vjust = 0.45, hjust = -0.2) +
  facet_nested(n+phi~var, scales = "free_x", 
               labeller = label_parsed, as.table = F,
               nest_line = element_line(color = "gray", linewidth = 1)) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text.x = element_text(colour = 'black', size = 14),
        strip.text.y = element_text(colour = 'black', size = 12),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y.right = element_blank(),
        panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_fill_viridis_d(name="",option = "G",direction = 1) +
  xlim(0,40) +
  ylab("") + xlab("Counts of Parameters Contained in Credible Intervals")

ggsave("plots/param_cov.png", width = 8, height = 6)

################################################################################
# prediction metrics ####
################################################################################
plot.data = sim.results %>%
  dplyr::select(c("n", "phi", "rep", "id", "model", "crps", "cov", "med.widths", "mae", "rmspe")) %>% 
  rename(Coverage = cov, 'Widths' = med.widths,CRPS = crps, MAE = mae, RMSPE = rmspe) %>% 
  mutate(phi = factor(paste0("phi == ", round(phi,2)),
                                levels = c("phi == 11.51", 
                                              "phi == 2.88")),
         n = factor(paste0("n == ", n), levels = c("n == 10000", "n == 2500"))) %>% 
  pivot_longer(cols = !c("n", "phi", "rep", "id", "model"), names_to = "var") %>% 
  mutate(model = gsub("_10", "(10)",model) %>% 
           gsub("_20", "(20)",.) %>% 
           gsub("_30", "(30)",.)) %>% 
  mutate(model = factor(model, levels = c("meshedGP", "blockNNGP", "NNGP(30)", "NNGP(20)",  "NNGP(10)",
                                          "cNNGP(30)", "cNNGP(20)", "cNNGP(10)")))  

cov.df = data.frame(var = "Coverage",
                    phi = rep(c("phi == 2.88","phi == 11.51"),each = 2),
                    n = rep(n.total,times = 2),
                    value = 0.95)

ggplot(plot.data, aes(x=value, y=model, fill = model)) + 
  geom_boxplot() +
  geom_vline(aes(xintercept = value), data = cov.df, 
             color = gray(0.2,0.6), linewidth = 1) +
  facet_nested(n+phi~factor(var, 
                            levels = c("Coverage", "Widths", "CRPS", "MAE", "RMSPE")), scales = "free_x", #independent = "x",
               labeller = label_parsed, as.table = F,
               nest_line = element_line(color = "gray", linewidth = 1)) +
    theme_minimal() +
    theme(legend.position = 'none',
          strip.text = element_text(colour = 'black', size = 12),
          axis.text.y = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle = 40),
          axis.line.x = element_blank(),
          axis.line.y.right = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
          panel.grid = element_blank(),
          panel.background = element_blank()) +
    scale_fill_viridis_d(name="",option = "G",direction = 1) +
  ylab("") + xlab("")

ggsave("prediction_metrics.png", width = 10, height = 6)

################################################################################
# WAIC
################################################################################
plot.data = sim.results %>%
  dplyr::select(c("n", "phi", "rep", "id", "model", "waic")) %>% 
  rename(WAIC = waic) %>% 
  mutate(phi = factor(paste0("phi == ", round(phi,2)),
                                levels = c("phi == 2.88", 
                                              "phi == 11.51")),
         n = factor(paste0("n == ", n), levels = c("n == 10000", "n == 2500"))) %>% 
  mutate(model = gsub("_10", "(10)",model) %>% 
           gsub("_20", "(20)",.) %>% 
           gsub("_30", "(30)",.)) %>% 
  mutate(model = factor(model, levels = c("meshedGP", "blockNNGP", "NNGP(30)", "NNGP(20)",  "NNGP(10)",
                                          "cNNGP(30)", "cNNGP(20)", "cNNGP(10)")), var = "WAIC")  

ggplot(plot.data, aes(x=WAIC, y=model, fill = model)) + 
  geom_boxplot() +
  facet_nested(n~phi, scales = "free_x", #independent = "x",
               labeller = label_parsed, as.table = F,
               nest_line = element_line(color = "gray", linewidth = 1)) +
    theme_minimal() +
    theme(legend.position = 'none',
          strip.text = element_text(colour = 'black', size = 12),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 8, color = "black", angle = 0),
          axis.line.x = element_blank(),
          axis.line.y.right = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
          panel.grid = element_blank(),
          panel.background = element_blank()) +
    scale_fill_viridis_d(name="",option = "G",direction = 1) +
  ylab("") + xlab("WAIC")

ggsave(filename = "plots/waic.png", width = 6, height = 4, dpi = 300)

################################################################################
###fitting time
################################################################################
plot.data = sim.results %>%
  dplyr::select(c("n", "phi", "rep", "id", "model", matches("fit|pred"))) %>% 
  mutate(phi = factor(paste0("phi == ", round(phi,2)),
                                levels = c("phi == 2.88",
                                           "phi == 11.51")),
         n = factor(paste0("n == ", n), levels = c("n == 2500", "n == 10000")),
         fit.time = fit.user.self + fit.sys.self + fit.sys.child + fit.user.child,
         pred.time = pred.user.self + pred.sys.self + pred.sys.child + pred.user.child) %>% 
  mutate(total.time = fit.time+pred.time)  %>% 
  select("n", "phi", "rep", "id", "model", "fit.time", "pred.time", "total.time")  %>% 
  pivot_longer(cols = !c("n", "phi", "rep", "id", "model"), names_to = "var") %>% 
  mutate(model = gsub("_10", "(10)",model) %>% 
           gsub("_20", "(20)",.) %>% 
           gsub("_30", "(30)",.)) %>% 
  mutate(model = factor(model, levels = c("meshedGP", "blockNNGP", "NNGP(30)", "NNGP(20)",  "NNGP(10)",
                                          "cNNGP(30)", "cNNGP(20)", "cNNGP(10)"))) %>% 
  filter(var == "fit.time")  %>% 
  mutate(value = value/(60*60)) 

med.df = plot.data %>% 
  group_by(model,n,phi) %>% 
  summarize(med.time = median(value),
    med.time.label = sprintf(median(value), fmt = '%#.1f'),
    max.time = max(value)) %>% 
  ungroup()

ggplot(plot.data, aes(x=value, y=model, fill = model)) + 
  geom_boxplot()+
  geom_text(data = med.df, 
            aes(y = model, x = max.time, label = med.time.label),
            hjust = -0.2) +
  geom_point(data = med.df, 
             aes(y = model, x = max.time*1.25), 
             alpha = 0) +
  facet_nested(.~n+phi, scales = "free_x",  #independent = "x",
               labeller = label_parsed, as.table = F,
               nest_line = element_line(color = "gray", linewidth = 1)) +
    theme_minimal() +
    theme(legend.position = 'none',
          strip.text = element_text(colour = 'black', size = 12),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 8, color = "black"),
          axis.line.x = element_blank(),
          axis.line.y.right = element_blank(),
          panel.border = element_rect(colour = "gray", fill=NA, linewidth=1),
          panel.grid = element_blank(),
          panel.background = element_blank()) +
    scale_fill_viridis_d(name="",option = "G",direction = 1) +
  ylab("") + xlab("Hours")

ggsave(filename = "plots/fit.times.png", width = 8, height = 3)

