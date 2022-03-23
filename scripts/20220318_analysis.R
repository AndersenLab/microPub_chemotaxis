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

# 2
fig1 <- ggplot(df %>% dplyr::filter(counter == "manual")) +
  aes(x = strain, y = CI) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(Assay)), width = .25, alpha = 0.75, size = 2)+
  ylim(-1,1) +
  facet_grid(~drug) +
  theme_classic() +
  theme(title = element_text(size = 16, color = "black", face = "bold"),
        axis.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        strip.text.x = element_text(size = 16, colour = "black")) + 
  labs(title = "FIG 1: Manual counting", color = "Block", x = "", y = "Chemotactic Index (ci)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position="right")
fig1

#ggsave('plots/fig1_manual_CI.pdf', width = 12, height = 8)
#ggsave('plots/fig1_manual_CI.png', width = 12, height = 8)

# Plotting automated counts

fig2 <- ggplot(df %>% dplyr::filter(counter == "auto")) +
  aes(x = strain, y = CI) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = as.factor(Assay)), width = .25, alpha = 0.75, size = 2)+
  facet_grid(~drug) +
  ylim(-1,1) +
  theme_classic() +
  theme(title = element_text(size = 16, color = "black", face = "bold"),
        axis.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        strip.text.x = element_text(size = 16, colour = "black")) + 
  labs(title = "FIG 2: Automated counting", color = "Block", x = "", y = "Chemotactic Index (ci)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position="right")
fig2

#ggsave('plots/fig2_auto_CI.pdf', width = 12, height = 8)
#ggsave('plots/fig2_auto_CI.png', width = 12, height = 8)

# Fitting correlation with ggplot to two compopunds manual counts vs. automated counts
  df_corr <- df %>%
  dplyr::select(-Drug_Count, -Control_Count, -count) %>%
  tidyr::spread(counter, CI) %>% 
  dplyr::rename(man_ci = manual, auto_ci = auto)

# Add two linear regression plots red = #F8766D, blue = #00BFC4
# show_col(hue_pal()(2))

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

cowplot::ggsave2(corr, filename = 'plots/assay_agreement_CI.png', width = 3.5, height = 3.5)
#ggsave('fig3_assay_agreement_CI.png', width = 12, height = 8)


# Calculate heritability including three replicates per block
h2_df_assay<- df %>%
  dplyr::rename(trait = drug, value = CI) %>%
  dplyr::select(-Drug_Count, -Control_Count, -Image_ID) %>% 
  dplyr::mutate(strain = as.character(strain)) %>%
  dplyr::arrange(counter, trait, strain) %>%
  dplyr::group_by(trait, counter, strain, Assay) %>%
  dplyr::mutate(value2 = mean(value)) %>%
  dplyr::ungroup()

# heritability functions
H2 <- function(d){
  strain.fit <- lme4::lmer(data = d, formula = value ~ 1 + (1|strain))
  variances <- lme4::VarCorr(x = strain.fit)
  A <- as.data.frame(variances)
  Vg <- A$vcov[1]
  Ve <- A$vcov[2]
  H2 <- Vg/(Vg+Ve)
  return(H2)
}

H2.bootstrapping.calc <- function(d, nreps = 100, boot = T){
  
  if(boot == T){
    # Broad-sense Heritability
    H2.point <- H2(d = d)
    #h2.point <- h2(d = d, geno_matrix = genos)
    H2.boots <- list()
    for(i in 1:nreps) {
      if(i %% 10 == 0){
        print(paste0((i/nreps)*100,"%"))
      }
      #################################
      # Bootstrap within strain ##
      #################################
      nested <- d %>%
        dplyr::group_by(strain) %>%
        tidyr::nest()
      boot.strain <- list()
      for(j in 1:length(nested$strain)){
        boot.strain[[j]] <- nested$data[[j]][sample(seq(1:nrow(nested$data[[j]])),replace = T),] %>%
          dplyr::mutate(strain = nested$strain[[j]])
      }
      boot <- boot.strain %>%
        Reduce(rbind,.)
      
      ##################################
      ## Bootstrap agnostic of strain ##
      ##################################
      # boot <- d[sample(seq(1:nrow(d)), replace = T),]
      
      check <- boot %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise(n())
      if(1 %in% check$`n()`){
        print("Only 1 Strain Sampled in Bootstrap - Skipping")
        next
      }
      # Broad-Sense Heritability
      H2.est <- H2(d = boot)
      H2.boots[i] <- H2.est
    }
    
    H2.boots.vec <- unlist(H2.boots)
    H2.quantiles <- c(quantile(H2.boots.vec, probs = seq(0,1,0.05)))
    H2.CI <- data.frame(H2.point, 
                        as.numeric(H2.quantiles[2]), 
                        as.numeric(H2.quantiles[21])) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    
    
    
    return(list(H2.CI,H2.boots.vec))
    
  } else {
    
    H2.point <- H2(d = d)
    # h2.point <- h2(d = d)
    H2.CI <- data.frame(H2.point, 
                        NA, 
                        NA) %>%
      `colnames<-`(c("H2.Point.Estimate","H2.5.perc","H2.95.perc"))
    return(H2.CI)
  }
  
}

# heritability caluculations
oct_H2_auto <- H2.bootstrapping.calc(h2_df_full %>% dplyr::filter(trait == "1-octanol" & counter == "auto") %>% dplyr::select(strain, value))
oct_H2_man <- H2.bootstrapping.calc(h2_df_full %>% dplyr::filter(trait == "1-octanol" & counter == "manual") %>% dplyr::select(strain, value))
iso_H2_auto <- H2.bootstrapping.calc(h2_df_full %>% dplyr::filter(trait == "isoamly alcohol" & counter == "auto") %>% dplyr::select(strain, value))
iso_H2_manual <- H2.bootstrapping.calc(h2_df_full %>% dplyr::filter(trait == "isoamly alcohol" & counter == "manual") %>% dplyr::select(strain, value))

# heritability summary
h2_dat <- tibble::tibble(type = c("oct_h2_auto","oct_H2_man","iso_H2_auto","iso_H2_man"), dplyr::bind_rows(oct_H2_auto[[1]], oct_H2_man[[1]], iso_H2_auto[[1]], iso_H2_manual[[1]])) %>%
  tidyr::separate(type, into = c("drug", "type", "counter"))

rio::export(h2_dat, "data/chemotaxis_heritability_summary.csv")

cor(df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("auto_ci"), df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("man_ci"))^2
lm(df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("auto_ci") ~ df_corr %>% tidyr::drop_na(.) %>% dplyr::pull("man_ci"))


# plot h2
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

cowplot::ggsave2(h2_plot, filename = "plots/heritability.png", width = 3.5, height = 3.5)
