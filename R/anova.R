anova_prot <- function(df, meta, control.name, covars = NULL){
  
  df.norm <- df
  
  rownames(df) <- df$Protein.IDs
  
  df <- select_if(df, is.numeric)
  
  fac <- as.factor(meta$Condition)
  
  name.treat = levels(fac)[levels(fac) != control.name]
  
  pred <- factor(meta$Condition, levels = c(control.name, name.treat))
  
  
  res <- DA.aov(data = as.matrix(df), predictor = pred, relative = FALSE, allResults = TRUE)
  
  tuk <- DA.TukeyHSD(res)
  
  
  namest <- unique(sub(".*_", "", colnames(tuk)))
  
  listy <- list()
  
  for (i in 1:length(namest)){
    
    listy[[i]] <- tuk %>% 
      select(ends_with(namest[i]))
    
  }
  
  names(listy) <- namest
  
  
  listy2 <- list()
  
  for (i in 1:length(listy)){
    
    p_val <- names(listy[[i]])[grep("pval_", names(listy[[i]]))]
    p.val.adj <- names(listy[[i]])[grep("pval.adj_", names(listy[[i]]))]
    
    dfn <- df.norm[which(df.norm$Protein.IDs %in% rownames(listy[[i]])),]
    
    dat_combine <- bind_cols(dfn, listy[[i]])
    
    name1 <- sub(".*-", "", names(listy[i]))
    name2 <- sub("-.*", "", names(listy[i]))
    
    names1key <- meta$key[which(meta$Condition == name1)]
    names2key <- meta$key[which(meta$Condition == name2)]
    
    dat_fc1 <- dat_combine %>%
      mutate(
        mean_control = rowMeans(df.norm[, names1key]),
        mean_trt = rowMeans(df.norm[, names2key]),
        log_fc = mean_control - mean_trt,
        log_pval = -1 * log10(.data[[p_val]]),
        log_pval_adj = -1 * log10(.data[[p.val.adj]])
      )
    
    dat_fc <- dat_fc1 %>%
      select(
        Protein.IDs,
        .data[[p_val]],
        .data[[p.val.adj]],
        log_fc,
        log_pval,
        log_pval_adj
      ) %>% 
      mutate(res = case_when((log_fc > 0) & (.data[[p.val.adj]] < 0.05) ~ "Up",
                             (log_fc < 0) & (.data[[p.val.adj]] < 0.05) ~ "Down",
                             TRUE ~ "no DE"))
    
    colnames(dat_fc) <- gsub("-", "_", colnames(dat_fc))
    colnames(dat_fc) <- gsub(" ", "_", colnames(dat_fc))
    
    listy2[[i]] <- dat_fc
    
  }
  
  names(listy2) <- namest
  
  
  return(listy2)
  
}


hist_prot <- function(df, var, xlab, bin.width = 0.05){
  
  fig <- df %>%
    ggplot(aes_string(var)) +
    geom_histogram(
      binwidth = bin.width,
      boundary = 0.5,
      fill = "darkblue",
      colour = "white"
    ) +
    xlab(xlab) +
    ylab("Frequency") +
    theme_classic()
  
  return(fig)
  
}




plot_res_aov <- function(df){
  
  var_pval <- colnames(df)[grep("^pval_", colnames(df))]
  
  var_adj <- colnames(df)[grep("^pval.adj_", colnames(df))]
  
  t1 <- hist_prot(df, var = var_pval, xlab = "P value")
  
  t2 <- hist_prot(df, var = var_adj, xlab = "adjusted P value")
  
  t3 <- hist_prot(df, var = "log_fc", xlab = "log2 fold change", bin.width = 0.5)
  
  fig <- ggpubr::ggarrange(t1, t2, t3, ncol = 2, nrow = 2)
  
  return(fig)
}
