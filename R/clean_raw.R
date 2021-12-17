filter_reference <- function(df, ref_channel = "Norm"){
  df <- df %>%
    filter(Condition != ref_channel)
  
  df$key <- gsub(" ", ".", df$key)
  
  return(df)
}

clean_key <- function(df){
  
  df$key <- make.names(df$key)
  
  return(df)
  
}

select_mq_vars <- function(df, meta){
  
  reporter_names_clean <- make.names(meta$key)
  
  df <- df %>%
    rename_all(make.names) %>%
    select(
      Protein.IDs,
      Q.value,
      matches(reporter_names_clean),
      Only.identified.by.site,
      Reverse,
      Potential.contaminant,
      Protein.names,
      Gene.names,
      Fasta.headers)
  
  return(df)
}

clean_mq_data <- function(df){
  
  df <- df %>% 
    mutate(Protein.IDs = sub(";.*", "", Protein.IDs))
  
  names_q <- df$Protein.IDs[which(df$Q.value > 0.01)]
  
  df <- df %>%
    select(-Q.value) %>%
    filter(is.na(Potential.contaminant)) %>%
    filter(is.na(Reverse)) %>%
    filter(is.na(Only.identified.by.site != "+")) %>%
    select(
      -Reverse,
      -Only.identified.by.site,
      -Potential.contaminant
    ) %>%
    mutate(
      across(where(is.numeric), na_if, 0),
      across(where(is.numeric), log2),
      across(where(is.numeric), na_if, Inf),
      across(where(is.numeric), na_if, -Inf),
    )
  
  names_q_after <- df$Protein.IDs[df$Protein.IDs %in% names_q]
  
  if (length(names_q_after) != 0){
    
    df <- df[-which(df$Protein.IDs %in% names_q_after),]
  }
  
  
  if (any(duplicated(na.omit(df$Gene.names)))) {
    
    df <- df %>% distinct(Gene.names, .keep_all = TRUE)
    
  }
  
  
  if (any(duplicated(na.omit(df$Protein.IDs)))) {
    
    df <- df %>% distinct(Protein.IDs, .keep_all = TRUE)
    
  }
  
  return(df)
  
}


select_no_mq <- function(df, meta){
  
  reporter_names_clean <- make.names(meta$key)
  
  df <- df %>% 
    rename_all(make.names) %>% 
    select(Accession,
           matches(reporter_names_clean)) %>% 
    rename(Protein.IDs = Accession)
  
  return(df)
  
}


clean_no_mq <- function(df){
  
  df <- df %>% 
    mutate(Protein.IDs = sub(";.*", "", Protein.IDs))
  
  df <- df %>%
    mutate(
      across(where(is.numeric), na_if, 0),
      across(where(is.numeric), log2)
    )
  
  
  if (any(duplicated(na.omit(df$Protein.IDs)))) {
    
    df <- df %>% distinct(Protein.IDs, .keep_all = TRUE)
    
  }
  
  return(df)
  
}