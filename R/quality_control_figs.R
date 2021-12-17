plot_n_prots_sample <- function(df, meta, fill = "Condition"){
  
  reporter_names_clean <- gsub(" ", ".", meta$key)
  
  df_plot <- df %>%
    select(
      Protein.IDs,
      matches(reporter_names_clean)
    ) %>%
    pivot_longer(!Protein.IDs, names_to = "key", values_to = "vals") %>%
    left_join(meta, by = "key") %>%
    select(
      key,
      vals,
      Condition,
      BioReplicate,
      Mixture,
      .data[[fill]]
    ) %>%
    filter(!is.na(vals)) %>%
    group_by(BioReplicate, .data[[fill]]) %>%
    count() %>%
    ggplot(aes(reorder(BioReplicate, n), n)) +
    geom_bar(aes_string(fill = fill), stat = "identity", color = "black") +
    xlab("Sample") +
    ylab("Count") +
    ggtitle("Number of proteins per sample") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) +
    ggsci::scale_fill_npg()
  
  return(df_plot)
  
}


mis_vals <- function(df){
  
  df <- df %>%
    select_if(is.numeric) %>%
    gather(key = "key", value = "val") %>%
    mutate(is.missing = is.na(val)) %>%
    group_by(key, is.missing) %>%
    summarise(num.missing = n()) %>%
    filter(is.missing == T) %>%
    select(-is.missing) %>%
    arrange(desc(num.missing))
  
  return(df)
}


mis_vals_pcts <- function(df, meta){
  
  df <- df %>%
    select_if(is.numeric) %>%
    gather(key = "key", value = "val") %>%
    mutate(isna = is.na(val)) %>%
    group_by(key) %>%
    mutate(total = n()) %>%
    group_by(key, total, isna) %>%
    summarise(num.isna = n()) %>%
    mutate(pct = num.isna / total * 100) %>%
    left_join(meta, by = "key") %>%
    select(key, total, isna, num.isna, pct, Condition, BioReplicate)
  
  levels <- (df %>% filter(isna == T) %>% arrange(desc(pct)))$BioReplicate
  
  pct.plot.NA <- df %>%
    ggplot() +
    geom_bar(aes(x = reorder(BioReplicate, -pct), y = pct, fill = isna),
             stat = "identity", alpha = 0.8
    ) +
    scale_fill_manual(
      name = "",
      values = c("gray43", "mediumpurple1"),
      labels = c("Present", "Missing")
    ) +
    coord_flip() +
    labs(x = "Samples", y = "% of missing values") +
    theme_minimal() +
    scale_x_discrete(labels = levels, limits = levels)
  
  return(pct.plot.NA)
  
}


na_heatmap <- function(df, meta, vars){
  
  plot_na <- df %>%
    select_if(is.numeric) %>%
    mutate_all(list(~ ifelse(is.na(.), 0, 1)))
  
  
  iterations = nrow(meta)
  variables = length(vars)
  
  my_sample_col <- matrix(ncol=variables, nrow=iterations)
  
  for(i in 1:variables){
    my_sample_col[,i] <- meta[[vars[i]]]
  }
  
  colnames(my_sample_col) <- vars
  my_sample_col <- as.data.frame(my_sample_col)
  
  rownames(my_sample_col) <- meta$key
  
  pheatmap(plot_na,
           color = c("skyblue", "grey"),
           show_rownames = FALSE,
           annotation_col = my_sample_col
  )
  
}


reorder_cormat <- function(cormat) {
  # Use correlation between variables as distance
  dd <- as.dist((1 - cormat) / 2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}


cor_plot <- function(df){
  
  cormat <- round(cor(na.omit(select_if(df, is.numeric))), 2)
  
  # Reorder the correlation matrix
  cormat <- reorder_cormat(cormat)
  
  # Melt the correlation matrix
  melted_cormat <- melt(cormat, na.rm = TRUE)
  
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue", high = "red", mid = "white",
      midpoint = 0, limit = c(-1, 1), space = "Lab",
      name = "Pearson\nCorrelation"
    ) +
    theme_minimal() + # minimal theme
    theme(axis.text.x = element_text(
      angle = 45, vjust = 1,
      size = 12, hjust = 1
    )) +
    coord_fixed() +
    xlab("") +
    ylab("")
  
  # Print the heatmap
  print(ggheatmap)
  
}


boxplot_ab <- function(df, meta, fill, reorder = TRUE, var_reorder = "counts"){
  
  if (reorder == TRUE){
    
    df <- df %>%
      select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
      pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
      left_join(meta, by = "key")
    
    df$BioReplicate <- reorder(as.factor(df$BioReplicate), df[[var_reorder]], na.rm = TRUE)
    
    
    plot <- ggplot(df, aes_string(x = "BioReplicate", y = "counts", fill = fill)) +
      geom_boxplot() +
      labs(y = "Log2 of abundance", x = "Samples") +
      scale_fill_npg() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12))
    
    return(plot)
    
  } else {
    
    plot <- df %>%
      select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
      pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
      left_join(meta, by = "key") %>%
      ggplot(aes_string(x = "BioReplicate", y = "counts", fill = fill)) +
      geom_boxplot() +
      labs(y = "Log2 of abundance", x = "Samples") +
      scale_fill_npg() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12))
    
    return(plot)
    
  }
  
} 


den_plot <- function(df, meta, color = "BioReplicate"){
  
  plot <- df %>%
    select(Protein.IDs, which(plyr::colwise(is.numeric)(.) == TRUE)) %>%
    pivot_longer(!Protein.IDs, names_to = "key", values_to = "counts") %>%
    left_join(meta, by = "key") %>%
    ggplot(aes_string(x = "counts", color = color)) +
    geom_density() +
    labs(y = "Log2 of abundance", x = "Samples") +
    scale_fill_npg() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 12)) +
    scale_fill_npg()
  
  return(plot)
  
}


filt_nas <- function(df, threshold = 0.7){
  
  filt <- df %>%
    mutate(filt = case_when(
      rowSums(!is.na(across(where(is.numeric)))) / ncol(across(where(is.numeric))) >= threshold ~ TRUE,
      TRUE ~ FALSE
    )) %>%
    filter(filt == TRUE) %>%
    select(-filt)
  
  return(filt)
  
}
