library(tidyverse)

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
