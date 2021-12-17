prep_data_bio <- function(df1, df2, meta, name.control){
  
  df1 <- as.data.frame(df1)
  # df1 <- na.omit(df1)
  df2 <- as.data.frame(df2)
  
  
  if ("Protein.IDs" %in% colnames(df1)){
    
    rownames(df1) <- df1$Protein.IDs
  }
  
  
  de_prot_df <- df1[which(rownames(df1) %in% df2$Protein.IDs[which(df2$res != "no DE")]), ]
  
  de_prot_df <- select_if(de_prot_df, is.numeric)
  
  de_prot_df <- as.data.frame(t(de_prot_df))
  
  vars <- setdiff(colnames(meta), "Condition")
  
  de_prot_df <- de_prot_df %>%
    mutate(key = rownames(.)) %>%
    inner_join(meta, by = "key") %>%
    select(- all_of(vars)) %>%
    select(Condition, everything())
  
  rownames(de_prot_df) <- colnames(select_if(df1, is.numeric))
  
  de_prot_df <- de_prot_df %>% 
    mutate(Condition = case_when(Condition == {{name.control}} ~ 0,
                                 TRUE ~ 1))
  
  de_prot_df <- as.data.frame(de_prot_df)
  
  return(de_prot_df)
}

data_part <- function(df, group, p = 0.7){
  
  set.seed(42)
  train <- createDataPartition(y = df[[group]], p = p, list = FALSE, times = 1)
  
  return(train)
  
}


create_train <- function(df, group, p = 0.7){
  
  datos_train <- df[data_part(df, group, p), ]
  
  return(datos_train)
}


create_test <- function(df, group, p = 0.7){
  
  datos_test <- df[-data_part(df, group, p), ]
  
  return(datos_test)
}

data_part_plot <- function(df_train, df_test, meta, name.control){
  
  fac <- as.factor(meta$Condition)
  
  name.treat = levels(fac)[levels(fac) != name.control] 
  
  v1 <- data.frame(prop.table(table(df_train$Condition)), Origin = "Train")
  v2 <- data.frame(prop.table(table(df_test$Condition)), Origin = "Test")
  
  vf <- rbind(v1, v2) %>% 
    rename(Condition = Var1) %>% 
    mutate(Condition = ifelse(Condition == 0, name.control, name.treat))
  
  figv <- ggplot(vf, aes(Condition, Freq, fill = Origin)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = wesanderson::wes_palette("Darjeeling2")) +
    theme_classic()
  
  
  h1 <- data.frame(table(df_train$Condition), Origin = "Train")
  h2 <- data.frame(table(df_test$Condition), Origin = "Test")
  
  hf <- rbind(h1, h2) %>% 
    rename(Condition = Var1) %>% 
    mutate(Condition = ifelse(Condition == 0, name.control, name.treat))
  
  figh <- ggplot(hf, aes(Condition, Freq, fill = Origin)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = wesanderson::wes_palette("Darjeeling2")) +
    theme_classic()
  
  fig <- ggpubr::ggarrange(figv, figh, common.legend = TRUE, legend = "bottom",
                           labels = c("Proportion", "Raw values"), hjust = -0.6)
  
  return(fig)
  
}


plot_auc <- function(perf, auc){
  
  perf_df <- data.frame(fpr = unlist(perf@x.values),
                        tpr = unlist(perf@y.values))
  
  fig <- ggplot(perf_df, aes(fpr, tpr)) + 
    geom_line(color = "lightsalmon2", size = 2) +
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1) +
    xlab(perf@x.name) + ylab(perf@y.name) +
    annotate("label", x = .95, y = 0.05, label = auc, colour = "lightsalmon4") 
  
  return(fig)
  
}



plot_cv_error <- function(cv_error){
  
  cv_error_df <- data.frame(lambda = cv_error$lambda,
                            cvm = cv_error$cvm,
                            cvsd = cv_error$cvsd,
                            vars = cv_error$nzero)
  
  vertical.lines <- c(log(min(cv_error$lambda)), log(cv_error$lambda.1se))
  # https://stats.stackexchange.com/questions/404795/interpretation-of-cross-validation-plot-for-lasso-regression
  
  fig <- ggplot(cv_error_df, aes(log(lambda), cvm)) +
    geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), color = "gray") +
    geom_point(color = "red") +
    theme_classic() +
    xlab("Log(lambda)") + ylab(cv_error$name) +
    geom_text(
      aes(label = vars, y = max(cv_error$cvm) + sd(cv_error$cvm)),
      size = 1.5
    ) +
    geom_vline(xintercept = vertical.lines, color = "grey", linetype = "dashed")
  
  return(fig)
}


plot_lasso_coef <- function(df_coeficientes){
  
  nColor <- nrow(df_coeficientes[which((df_coeficientes$coeficiente != 0) & (df_coeficientes$predictor != "(Intercept)")),])
  myColor <- randomcoloR::distinctColorPalette(k = nColor)
  
  
  fig <- df_coeficientes %>%
    filter(predictor != "(Intercept)") %>%
    ggplot(aes(x = predictor, y = coeficiente, fill = ifelse(coeficiente != 0, predictor, "no"))) +
    geom_col(color = "black") +
    labs(title = "Lasso's coefficients") +
    theme_classic() +
    xlab("Predictor") + ylab("Coefficient") +
    theme(axis.text.x = element_text(size = 6, angle = 90),
          legend.position = "none") +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    scale_fill_manual(values = c(myColor, "gray"))
  
  return(fig)
  
}


conf_matrix <- function(preds, observed, meta, name.control){
  
  trt = setdiff(levels(as.factor(meta$Condition)), name.control)
  
  preds <- as.data.frame(preds)
  
  preds <- preds %>% 
    mutate(obs = observed) %>% 
    rename(pred = s0) %>% 
    select(obs,
           pred)
  
  preds$obs <- as.factor(ifelse(preds$obs == 0, name.control, trt))
  preds$pred <- as.factor(ifelse(preds$pred == 0, name.control, trt))
  
  cm <- conf_mat(preds, obs, pred)
  
  fig <- autoplot(cm, type = "heatmap") +
    scale_fill_gradient(low = "pink", high = "cyan") +
    ylab("Predicted") + xlab("Observed")
  
  return(fig)
  
}


all_lasso <- function(m_train, labels_train, m_test, labels_test, meta, name.control, standarize = FALSE){
  
  m_train <- as.matrix(m_train)
  m_test <- as.matrix(m_test)
  
  modelo <- glmnet(m_train, labels_train, family = "binomial", alpha = 1, standardize = standarize)
  
  regularizacion <- modelo$beta %>%
    as.matrix() %>%
    t() %>%
    as_tibble() %>%
    mutate(lambda = modelo$lambda)
  
  regularizacion <- regularizacion %>%
    pivot_longer(
      cols = !lambda,
      names_to = "predictor",
      values_to = "coeficientes"
    )
  
  fig1 <- regularizacion %>%
    ggplot(aes(x = lambda, y = coeficientes, color = predictor)) +
    geom_line(size = .8) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(title = "Regularized model coefficients") +
    theme_bw() +
    theme(legend.position = "none") + ylab("Coefficients") +
    xlab("Lambda")
  
  
  
  set.seed(123)
  cv_error <- cv.glmnet(m_train, labels_train,
                        family = "binomial", alpha = 1,
                        nfolds = 10, type.measure = "mse", standardize = standarize)
  
  
  
  fig2 <- plot_cv_error(cv_error)
  
  
  modelo_l <- glmnet(m_train, labels_train, family = "binomial", alpha = 1, standardize = standarize,
                     lambda = cv_error$lambda.1se)
  
  
  df_coeficientes <- coef(modelo_l) %>%
    as.matrix() %>%
    as_tibble(rownames = "predictor") %>%
    rename(coeficiente = s0)
  
  fig3 <- plot_lasso_coef(df_coeficientes)
  
  df_lasso <- df_coeficientes %>%
    rename(Coefficients = coeficiente,
           Predictor = predictor) %>% 
    filter(
      Predictor != "(Intercept)",
      Coefficients != 0)
  
  
  names_lasso <- df_coeficientes %>%
    filter(
      predictor != "(Intercept)",
      coeficiente != 0) %>%
    select(predictor) %>%
    pull()
  
  
  predicciones_train_c <- predict(modelo_l, newx = m_train, type = "class")
  
  m1 <- conf_matrix(predicciones_train_c, labels_train, meta, name.control)
  
  predicciones_train <- predict(modelo_l, newx = m_train)
  
  training_mse <- mean((predicciones_train - labels_train)^2)
  
  
  
  predicciones_test_c <- predict(modelo_l, newx = m_test, type = "class")
  
  m2 <- conf_matrix(predicciones_test_c, labels_test, meta, name.control)
  
  predicciones_test <- predict(modelo_l, newx = m_test, type = "response")
  
  test_mse_lasso <- mean((predicciones_test - labels_test)^2)
  
  
  
  perf <- performance(prediction(predicciones_test, labels_test), "tpr", "fpr")
  
  auc <- paste0("AUC = ", round(unlist(performance(prediction(predicciones_test, labels_test),
                                                   measure = "auc")@y.values), 4) * 100, "%")
  
  fig4 <- plot_auc(perf, auc)
  
  return(list(model = modelo_l,
              plot_first_model = fig1,
              fig_cv_error = fig2,
              fig_coefs = fig3,
              prot_names = names_lasso,
              df_lasso = df_lasso,
              conf_matrix_train = m1,
              conf_matrix_test = m2,
              roc_auc = fig4))
  
}


plot_auc_all <- function(names_lasso, datos_train, x_test, y_test){
  
  allModels <- lapply(names_lasso, function(x) {
    predict.glm(glm(formula = paste0("Condition ~ `", x, "`"), data = datos_train, family = "binomial"),
                newdata = as.data.frame(x_test), type = "response"
    )
  })
  
  names(allModels) <- names_lasso
  
  nColor <- length(names_lasso)
  myColor <- randomcoloR::distinctColorPalette(k = nColor)
  
  names.legend <- vector()
  for (j in 1:length(allModels)) {
    pred <- prediction(
      as.vector(allModels[[j]], mode = "numeric"),
      y_test
    )
    auc <- round(
      unlist(performance(pred, measure = "auc")@y.values),
      4) * 100
    
    names.legend[j] <- paste0(names_lasso[j], " - AUC = ", auc, "%")
  }
  
  # Plot every ROC overlapped
  
  perf <- list()
  for (i in 1:length(names_lasso)) {
    pred <- prediction(
      as.vector(allModels[[i]], mode = "numeric"),
      y_test
    )
    perf[[i]] <- performance(pred, "tpr", "fpr")
    names(perf)[i] <- names_lasso[i]
    
  }
  

  
  perf_df <- m_perf_df(names_lasso, names.legend, perf)
  
  fig <- ggplot(perf_df, aes(fpr, tpr)) + 
    geom_line(aes(color = prots), size = 1.5) +
    theme_classic() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", size = 1) +
    scale_color_manual(values = myColor) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.background = element_rect(size=0.5, linetype="solid", colour = "black"),
          legend.text=element_text(size=10)) +
    xlab("False positive rate (FPR)") +
    ylab("True positive rate (TPR)")
  
  return(fig)
  
}


m_perf_df <- function(names_lasso, names.legend, perf){
  
  # get length of the number of values within perf: fdr
  
  vec_fpr <- vector()
  for (i in 1:length(perf)){
    
    vec_fpr[i] <- length(unlist(perf[[i]]@x.values))
    
  }
  
  
  iterations = length(names_lasso)
  variables = max(vec_fpr)
  
  fpr <- matrix(ncol=iterations, nrow=variables)
  
  for(i in 1:iterations){
    
    if(length(unlist(perf[[i]]@x.values)) != variables){
      
      n <- variables - length(unlist(perf[[i]]@x.values))
      zeroes <- rep(0, n)
      
      y <- c(zeroes, unlist(perf[[i]]@x.values))
      
      fpr[,i] <- y
      
    } else {
      
      fpr[,i] <- unlist(perf[[i]]@x.values)
      
    }
  }
  
  colnames(fpr) <- names.legend
  fpr <- as.data.frame(fpr)
  
  fpr <- fpr %>% 
    mutate(typ = "fpr") %>% 
    select(typ, everything()) %>% 
    pivot_longer(!typ, names_to = "prots", values_to = "fpr") %>% 
    select(-typ)
  
  
  # get length of the number of values within perf: tpr
  
  vec_tpr <- vector()
  for (i in 1:length(perf)){
    
    vec_tpr[i] <- length(unlist(perf[[i]]@y.values))
    
  }
  
  tpr <- matrix(ncol=iterations, nrow=variables)
  
  for(i in 1:iterations){
    
    if(length(unlist(perf[[i]]@y.values)) != variables){
      
      n <- variables - length(unlist(perf[[i]]@y.values))
      zeroes <- rep(0, n)
      
      y <- c(zeroes, unlist(perf[[i]]@y.values))
      
      tpr[,i] <- y
      
    } else {
      
      tpr[,i] <- unlist(perf[[i]]@y.values)
      
    }
  }
  
  colnames(tpr) <- names.legend
  tpr <- as.data.frame(tpr)
  
  tpr <- tpr %>% 
    mutate(typ = "tpr") %>% 
    select(typ, everything()) %>% 
    pivot_longer(!typ, names_to = "prots", values_to = "tpr") %>% 
    select(-typ)
  
  
  perf_df <- cbind(fpr, tpr)[,-3]
  
  return(perf_df)
  
}


pca_lasso <- function(df, meta, names_lasso, scale = FALSE, labels = "point"){
  
  df <- df[,which(colnames(df) %in% names_lasso)]
  
  rownames(df) <- meta$key
  
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
  
  out_vars <- setdiff(out_vars, "Condition")
  
  tt <- as.data.frame(select_if(df, is.numeric))
  
  rownames(tt) <- meta$key
  
  tt <- tt %>%
    mutate(key = rownames(.)) %>%
    inner_join(meta, by = "key") %>%
    select(-all_of(out_vars)) %>%
    select(Condition, everything())
  
  rownames(tt) <- rownames(df, is.numeric)
  
  res.pca <- prcomp(select_if(tt, is.numeric), scale = FALSE)
  
  fviz_pca_ind(res.pca,
               geom.ind = labels,
               habillage = tt$Condition,
               legend.title = "", title = "",
               invisible = "quali",
               pointsize = 4,
               pointshape = 19
  ) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    scale_color_npg()
  
}


heatmap_lasso <- function(df, meta, names_lasso, z.score = TRUE, cutrow = NA, cutcol = NA){
  
  df <- df[,which(colnames(df) %in% names_lasso)]
  
  df <- df %>%
    t() %>% 
    as.data.frame(.)
  
  df_heat <- select_if(df, is.numeric)
  
  colnames(df_heat) <- meta$BioReplicate
  
  
  data_subset_norm <- t(apply(df_heat, 1, cal_z_score))
  
  my_sample_col <- data.frame(Condition = meta$Condition)
  rownames(my_sample_col) <- meta$BioReplicate
  
  pheatmap(data_subset_norm, annotation_col = my_sample_col, cutree_rows = cutrow, cutree_cols = cutcol)
  
}
