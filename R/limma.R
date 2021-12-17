
get.covars <- function(covars, meta){
  
  list.covars <- list()
  
  if (is.null(covars)){
    
    return(NULL)
    
  } else {
    
    for (i in 1:length(covars)){
      
      list.covars[[i]] <- as.factor(meta[[covars[i]]])
      
    }
    
    names(list.covars) <- covars
    
    return(list.covars)
    
  }
}
    
    



de_limma <- function(df, meta, control.name, covars = NULL, coeff = 2){
  
  rownames(df) <- df$Protein.IDs

  df <- select_if(df, is.numeric)

  fac <- as.factor(meta$Condition)

  name.treat = levels(fac)[levels(fac) != control.name]

  pred <- factor(meta$Condition, levels = c(name.treat, control.name))
  
  
  list.covars <- get.covars(covars, meta)
  

  res <- DA.lim(data = as.matrix(df), predictor = pred, relative = FALSE)


  res <- res %>%
    select(- ordering,
           - Method) %>%
    mutate(res = case_when((logFC > 0) & (pval.adj < 0.05) ~ "Up",
                           (logFC < 0) & (pval.adj < 0.05) ~ "Down",
                           TRUE ~ "no DE")) %>%
    rename(Protein.IDs = Feature,
           log_fc = logFC,
           p.val.adj = pval.adj,
           p_val = pval) %>%
    mutate(log_pval_adj = -1 * log10(p.val.adj))
  
  return(res)
  
}