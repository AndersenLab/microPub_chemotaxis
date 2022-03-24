library(tidyverse)

# Set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# read data
df <- as.data.frame(read.csv(file = glue::glue("data/Results_Block_3_4_6b_8.csv"), header = T)) %>%
  dplyr::select(-X, -X.1, -X.2) %>%
  dplyr::filter(CI != "#DIV/0!") %>%
  dplyr::filter(strain != "ED3017") %>%
  dplyr::filter(!is.na(Assay),
                Assay != "20180307") %>% #filter block 20180307
  dplyr::mutate(CI = as.numeric(as.character(CI)),
                drug = ifelse(drug == "isa", "isoamly alcohol", "1-octanol"))

glimpse(df)

#===================================================#
# Method correlation
#===================================================#
# Shape data to fit lm with ggplot::geom_smooth to manual counts vs. automated counts
df_corr <- df %>%
  dplyr::select(-Drug_Count, -Control_Count, -count) %>%
  tidyr::spread(counter, CI) %>% 
  dplyr::rename(man_ci = manual, auto_ci = auto)

# Plot the correlation between methods
corr <- ggplot2::ggplot(df_corr) +
  ggplot2::aes(x = man_ci, y = auto_ci) +
  ggplot2::geom_point(aes(x = man_ci, y = auto_ci, color = as.factor(drug))) + 
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
  ggplot2::geom_smooth(method=lm, size = 0.5, se = F) +
  ggpubr::stat_regline_equation(label.y = 0.75, aes(label = ..eq.label..)) +
  ggpubr::stat_regline_equation(label.y = 0.5, aes(label = ..rr.label..)) +
  ggplot2::xlim(-1,1) +
  ggplot2::ylim(-1,1) +
  theme_bw() +
  ggplot2::labs(title = "", color = "compound", x = "Manual (CI)", y = "Automated (CI)") +
  ggplot2::theme(legend.position= c(0.8, 0.2))
corr

# double check ggplot::geom_smooth calculated correct lm and r^2
cor(df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("auto_ci"), df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("man_ci"))^2
lm(df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("auto_ci") ~ df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("man_ci"))

# save the correlation plot
cowplot::ggsave2(corr, filename = 'plots/assay_agreement_CI.png', width = 3.5, height = 3.5)

#========================================================#
# Heritability 
#========================================================#
# Shape data to calculate heritability including three replicates per block
h2_df_assay <- df %>%
  dplyr::rename(trait = drug, value = CI) %>%
  dplyr::select(-Drug_Count, -Control_Count, -Image_ID) %>% 
  dplyr::mutate(strain = as.character(strain)) %>%
  dplyr::arrange(counter, trait, strain) %>%
  dplyr::group_by(trait, counter, strain, Assay) %>%
  dplyr::mutate(value2 = mean(value)) %>%
  dplyr::ungroup()

# Source heritability functions
source("scripts/heritability_functions.R")

# heritability caluculations
oct_H2_auto <- H2.bootstrapping.calc(h2_df_assay %>% dplyr::filter(trait == "1-octanol" & counter == "auto") %>% dplyr::select(strain, value))
oct_H2_man <- H2.bootstrapping.calc(h2_df_assay %>% dplyr::filter(trait == "1-octanol" & counter == "manual") %>% dplyr::select(strain, value))
iso_H2_auto <- H2.bootstrapping.calc(h2_df_assay %>% dplyr::filter(trait == "isoamly alcohol" & counter == "auto") %>% dplyr::select(strain, value))
iso_H2_manual <- H2.bootstrapping.calc(h2_df_assay %>% dplyr::filter(trait == "isoamly alcohol" & counter == "manual") %>% dplyr::select(strain, value))

# Build heritability summary dataframe
h2_dat <- tibble::tibble(type = c("oct_h2_auto","oct_H2_man","iso_H2_auto","iso_H2_man"), dplyr::bind_rows(oct_H2_auto[[1]], oct_H2_man[[1]], iso_H2_auto[[1]], iso_H2_manual[[1]])) %>%
  tidyr::separate(type, into = c("drug", "type", "counter"))

# Export heritability data
rio::export(h2_dat, "data/chemotaxis_heritability_summary.csv")

# plot hetitability
h2_plot <- ggplot(h2_dat) +
  geom_point(aes(x = counter, y = H2.Point.Estimate), size = 3) +
  geom_errorbar(aes(ymin=H2.5.perc, ymax=H2.95.perc, x = counter), colour="black", width=0.05) +
  labs(x = "", y = "broad-sense heritability") +
  ylim(0, 1) +
  theme_bw()+
  facet_grid(~drug) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = .2, color= "dark red", lty = 2)
h2_plot

# save heritability plot
cowplot::ggsave2(h2_plot, filename = "plots/heritability.png", width = 3.5, height = 3.5)

#=================================================================#
# old code and plots
#================================================================#
# # 2
# fig1 <- ggplot(df %>% dplyr::filter(counter == "manual")) +
#   aes(x = strain, y = CI) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = as.factor(Assay)), width = .25, alpha = 0.75, size = 2)+
#   ylim(-1,1) +
#   facet_grid(~drug) +
#   theme_classic() +
#   theme(title = element_text(size = 16, color = "black", face = "bold"),
#         axis.title = element_text(size = 16, color = "black", face = "bold"),
#         axis.text = element_text(size = 16, color = "black"),
#         legend.text = element_text(size = 16, color = "black"),
#         strip.text.x = element_text(size = 16, colour = "black")) + 
#   labs(title = "FIG 1: Manual counting", color = "Block", x = "", y = "Chemotactic Index (ci)") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(legend.position="right")
# fig1
# 
# #ggsave('plots/fig1_manual_CI.pdf', width = 12, height = 8)
# #ggsave('plots/fig1_manual_CI.png', width = 12, height = 8)
# 
# # Plotting automated counts
# 
# fig2 <- ggplot(df %>% dplyr::filter(counter == "auto")) +
#   aes(x = strain, y = CI) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(aes(color = as.factor(Assay)), width = .25, alpha = 0.75, size = 2)+
#   facet_grid(~drug) +
#   ylim(-1,1) +
#   theme_classic() +
#   theme(title = element_text(size = 16, color = "black", face = "bold"),
#         axis.title = element_text(size = 16, color = "black", face = "bold"),
#         axis.text = element_text(size = 16, color = "black"),
#         legend.text = element_text(size = 16, color = "black"),
#         strip.text.x = element_text(size = 16, colour = "black")) + 
#   labs(title = "FIG 2: Automated counting", color = "Block", x = "", y = "Chemotactic Index (ci)") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(legend.position="right")
# fig2
# 
# #ggsave('plots/fig2_auto_CI.pdf', width = 12, height = 8)
# #ggsave('plots/fig2_auto_CI.png', width = 12, height = 8)
