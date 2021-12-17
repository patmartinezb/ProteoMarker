pca_prot <- function(df, meta, var_p, scale = FALSE, labels = "point"){
  
  df <- na.omit(df)
  
  vars <- meta %>% 
    select(-Run,
           -TechRepMixture,
           -Fraction,
           -Channel,
           -BioReplicate,
           -key) %>% 
    select_if(is.character) %>% 
    colnames(.)
  
  out_vars <- colnames(meta)
  
  out_vars <- setdiff(out_vars, vars)
  
  tt <- as.data.frame(t(select_if(df, is.numeric)))
  
  tt <- tt %>%
    mutate(key = rownames(.)) %>%
    inner_join(meta, by = "key") %>%
    select(-all_of(out_vars)) %>%
    select(all_of(vars), everything())
  
  rownames(tt) <- colnames(select_if(df, is.numeric))
  
  res.pca <- prcomp(select_if(tt, is.numeric), scale = FALSE)
  
  fviz_pca_ind(res.pca,
               geom.ind = labels,
               habillage = tt[[var_p]],
               legend.title = "", title = "",
               invisible = "quali",
               pointsize = 4,
               pointshape = 19
  ) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    scale_color_npg()
  
}




cal_z_score <- function(x) {
  (x - mean(x)) / sd(x)
}


heatmap_prot <- function(df, meta, vars, z.score = TRUE, cutrow = NA, cutcol = NA){
  
  df <- as.data.frame(df)
  meta <- as.data.frame(meta)
  
  iterations = nrow(meta)
  variables = length(vars)
  
  my_sample_col <- matrix(ncol=variables, nrow=iterations)
  
  for(i in 1:variables){
    my_sample_col[,i] <- meta[[vars[i]]]
  }
  
  colnames(my_sample_col) <- vars
  my_sample_col <- as.data.frame(my_sample_col)
  
  rownames(my_sample_col) <- meta$BioReplicate
  
  
  
  df_heat <- select_if(df, is.numeric)
  
  colnames(df_heat) <- meta$BioReplicate
  
  if (z.score == TRUE){
    
    df_heat <- t(apply(df_heat, 1, cal_z_score))
    
  }
  
  pheatmap(na.omit(df_heat),
           annotation_col = my_sample_col,
           cutree_rows = cutrow,
           cutree_cols = cutcol
  )
  
}



pvca_plot <- function(df, meta){
  
  df <- na.omit(df)
  df <- df[meta$key]
  
  vars <- meta %>% 
    select(-Run,
           -TechRepMixture,
           -Fraction,
           -Channel,
           -BioReplicate,
           -key) %>% 
    colnames(.)
  
  
  
  v <- vector()
  
  for (i in 1:length(vars)){
    
    if (length(levels(as.factor(meta[[vars[i]]]))) == 1){
      
      v[i] <- vars[i]
    }
  }
  
  
  if (length(v) == 1){
    
    stop("PVCA cannot be computed with only one factor")
    
  }
  
  
  vars <- setdiff(vars, v)
  
  met <- meta %>%
    select(
      key,
      all_of(vars)
    ) %>%
    as.data.frame()
  
  
  met <- as.data.frame(lapply(met, as.factor))
  
  rownames(met) <- met$key
  
  met <- select(met, -key)
  
  phenoData <- new("AnnotatedDataFrame",
                   data = met
  )
  
  exampleSet <- ExpressionSet(
    assayData = as.matrix(select_if(df, is.numeric)),
    phenoData = phenoData
  )
  
  pct_threshold <- 0.6
  batch.factors <- vars
  pvcaObj <- pvcaBatchAssess(exampleSet, batch.factors, pct_threshold)
  
  
  bp <- barplot(pvcaObj$dat,
                xlab = "Effects",
                ylab = "Weighted average proportion variance",
                ylim = c(0, 1.1), col = c("blue"), las = 2,
                main = ""
  )
  
  axis(1, at = bp, labels = pvcaObj$label, xlab = "Effects", cex.axis = 0.7, las = 2)
  title("PVCA estimation bar chart", line = 0)
  values <- pvcaObj$dat
  new_values <- round(values, 3)
  text(bp, pvcaObj$dat, labels = new_values, pos = 3, cex = 0.8)
  
} # Check aesthetics




batch_rm <- function(df, meta, batch1 = NA, batch2 = NA){
  
  prots <- df$Protein.IDs
  
  if (is.na(batch2)){
    
    set.seed(42)
    df_wb <- removeBatchEffect(select_if(df, is.numeric), batch = meta[[batch1]])
    
    df_wb <- as.data.frame(df_wb)
    rownames(df_wb) <- prots
    
  } else {
    
    set.seed(42)
    df_wb <- removeBatchEffect(select_if(df, is.numeric), batch = meta[[batch1]], batch2 = meta[[batch2]])
    
    df_wb <- as.data.frame(df_wb)
    rownames(df_wb) <- prots
    
  }
  
  return(df_wb)
  
}