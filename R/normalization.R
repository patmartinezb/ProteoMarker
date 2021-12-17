normalize_prot <- function(df){
  
  prots <- df$Protein.IDs
  
  set.seed(42)
  df_norm <- normalize.quantiles(as.matrix(select_if(df, is.numeric)))
  
  colnames(df_norm) <- colnames(select_if(df, is.numeric))
  
  df_norm <- df_norm %>%
    as.data.frame() %>%
    mutate(Protein.IDs = prots) %>% 
    select(Protein.IDs, everything())
  
  return(df_norm)
  
} # use switch to allow for other normalization methods within shiny


norm_prot <- function(df, normMethod){ 
  
  prots <- df$Protein.IDs
  
  normData <- select_if(df, is.numeric)
  
  df <- switch(normMethod,
               "medianNorm" = medianNorm(normData),
               "meanNorm" = meanNorm(normData),
               "quantNorm" = quantNorm(normData),
               
  )
  
  df <- as.data.frame(df)
  
  df <- df %>% 
    mutate(Protein.IDs = prots) %>% 
    select(Protein.IDs, everything())
  return(df)
}


norm_boxplot <- function(df1, df2, meta){
  
  p1 <- df1 %>%
    select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
    left_join(meta, by = "key") %>%
    ggplot(aes(x = BioReplicate, y = counts, fill = Condition)) +
    geom_boxplot() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_npg() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12))
  
  p2 <- df2 %>%
    select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
    left_join(meta, by = "key") %>%
    ggplot(aes(x = BioReplicate, y = counts, fill = Condition)) +
    geom_boxplot() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_npg() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12))
  
  fig <- ggpubr::ggarrange(p1, p2, labels = c("Pre-normalization", "Post-normalization"),
                           common.legend = TRUE, legend = "bottom")
  
  return(fig)
  
}
