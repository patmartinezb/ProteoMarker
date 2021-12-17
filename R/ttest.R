t_test <- function(df, grp1, grp2) {
  # Subset control group and convert to numeric
  x <- df[grp1] %>%
    unlist() %>%
    as.numeric()
  # Subset treatment group and convert to numeric
  y <- df[grp2] %>%
    unlist() %>%
    as.numeric()
  # Perform t-test using the mean of x and y
  result <- t.test(x, y)
  # Extract p-values from the results
  p_vals <- tibble(p_val = result$p.value)
  # Return p-values
  return(p_vals)
}


ttest_prot <- function(df, df_norm, meta, name.con){
  
  df <- select_if(as.data.frame(df), is.numeric)
  
  fac <- as.factor(meta$Condition)
  
  name.treat = levels(fac)[levels(fac) != name.con]
  
  dat_pvals <- plyr::adply(as.matrix(df),
                           .margins = 1, .fun = t_test,
                           grp1 = which(fac == name.con),
                           grp2 = which(fac == name.treat)) %>%
    as_tibble() %>%
    mutate(p.val.adj = p.adjust(p_val, method = "BH"))
  
  
  # Select columns and log data
  
  dat_combine <- bind_cols(dat_pvals[, 1], df, dat_pvals[, 2:3])
  
  dat_fc1 <- dat_combine %>%
    mutate(
      mean_control = rowMeans(df[, which(fac == name.con)]),
      mean_trt = rowMeans(df[, which(fac == name.treat)]),
      log_fc = mean_control - mean_trt,
      log_pval = -1 * log10(p_val),
      log_pval_adj = -1 * log10(p.val.adj)
    )
  
  
  dat_fc <- dat_fc1 %>%
    select(
      X1,
      p_val,
      p.val.adj,
      log_fc,
      log_pval,
      log_pval_adj
    ) %>% 
    mutate(res = case_when((log_fc > 0) & (p.val.adj < 0.05) ~ "Up",
                           (log_fc < 0) & (p.val.adj < 0.05) ~ "Down",
                           TRUE ~ "no DE")) %>% 
    rename(Protein.IDs = X1) %>% 
    mutate(Protein.IDs = df_norm$Protein.IDs)
    
  
  return(dat_fc)
  
}


table_des <- function(df){
  
  res <- plyr::count(df, "res")
  
  return(res)
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


plot_res <- function(df){
  
  t1 <- hist_prot(df, var = "p_val", xlab = "P value")
  
  t2 <- hist_prot(df, var = "p.val.adj", xlab = "adjusted P value")
  
  t3 <- hist_prot(df, var = "log_fc", xlab = "log2 fold change", bin.width = 0.5)
  
  fig <- ggpubr::ggarrange(t1, t2, t3, ncol = 2, nrow = 2)
  
  return(fig)
}


volcano_ttest <- function(df){
  
  df <- df %>%
    mutate(threshold = case_when(
      (log_fc >= 1 & log_pval_adj >= 1.3) | (log_fc <= -1 & log_pval_adj >= 1.3) ~ "A",
      (log_fc < 1 & log_pval_adj >= 1.3) | (log_fc > -1 & log_pval_adj >= 1.3) ~ "C",
      TRUE ~ "B"
    ))
  
  ggplot(df, aes(log_fc, log_pval_adj, color = threshold)) +
    geom_point(alpha = 0.5) +
    geom_label_repel(
      data = df %>% filter(threshold == "A"), # Filter data first
      aes(label = Protein.IDs)
    ) + 
    geom_hline(yintercept = 1.3, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = 1, linetype = 2, alpha = 0.5) +
    geom_vline(xintercept = -1, linetype = 2, alpha = 0.5) +
    scale_colour_manual(values = c("A" = "red", "B" = "black", "C" = "blue")) +
    xlab("Log2 FC") +
    ylab("-Log10 adj P value") +
    theme_minimal() +
    theme(legend.position = "none") +
    xlim(-2, 2)
  
} # Maybe allow changing thresholds
