
impute_nas <- function(normData, imputMethod){ 
  
  prots <- normData$Protein.IDs
  
  normData <- select_if(normData, is.numeric)
  
  df <- switch(imputMethod,
               "KNN" = impute.knn(as.matrix(normData), rowmax = 0.9)$data,
               # "QRILC" = imputeLCMD::impute.QRILC(normData)[[1]], # need to normalize beforehand
               "MinDet" = imputeLCMD::impute.MinDet(normData),
               "MinProb" = imputeLCMD::impute.MinProb(normData),
               "Min" = {
                 normData[is.na(normData)] = min(normData, na.rm = T)
                 normData
               }
  )
  
  df <- as.data.frame(df)
  
  df <- df %>% 
    mutate(Protein.IDs = prots) %>% 
    select(Protein.IDs, everything())
  return(df)
}


no_impute <- function(df){
  
  prots <- df$Protein.IDs
  
  df.imp <- df %>%
    select_if(is.numeric) %>% 
    mutate(Protein.IDs = prots) %>% 
    select(Protein.IDs, everything())
  
  return(df.imp)
  
}


compare_den <- function(df1, df2, meta, color = "BioReplicate"){
  
  h1 <- df1 %>%
    select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
    left_join(meta, by = "key") %>%
    ggplot(aes_string(x = "counts", color = color)) +
    geom_density() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_npg() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12),
      legend.position = "bottom"
    )
  
  
  
  h2 <- df2 %>%
    select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
    left_join(meta, by = "key") %>%
    ggplot(aes_string(x = "counts", color = color)) +
    geom_density() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_npg() +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12),
      legend.position = "bottom"
    )
  
  fig <- ggpubr::ggarrange(h1, h2, labels = c("Pre-imputation", "Post-imputation"),
                           common.legend = TRUE)
  
  return(fig)
  
}
