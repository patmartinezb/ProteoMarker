
medianNorm <- function(logDat) {
  rowmed <- apply(logDat, 1, function(x){median(x, na.rm = TRUE)})
  
  out <- sweep(logDat, 1, rowmed, "-")
  
  return(as.matrix(out))
}

meanNorm <- function(logDat) {
  rowmean <- apply(logDat, 1, function(x){mean(x, na.rm = TRUE)})
  
  out <- sweep(logDat, 1, rowmean, "-")
  
  return(as.matrix(out))
}

# vsnNorm <- function(dat) {
#   vsnNormed <- suppressMessages(vsn::justvsn(as.matrix(dat)))
#   colnames(vsnNormed) <- colnames(dat)
#   row.names(vsnNormed) <- rownames(dat)
#   return(as.matrix(vsnNormed))
# }

quantNorm <- function(logDat) {
  quantNormed <- preprocessCore::normalize.quantiles(as.matrix(logDat), copy=FALSE)
  colnames(quantNormed) <- colnames(logDat)
  row.names(quantNormed) <- rownames(logDat)
  return(as.matrix(quantNormed))
}

# cycLoessNorm <- function(logDat) {
#   cycLoessNormed <- limma::normalizeCyclicLoess(as.matrix(logDat), method="fast")
#   colnames(cycLoessNormed) <- colnames(logDat)
#   row.names(cycLoessNormed) <- rownames(logDat)
#   return(as.matrix(cycLoessNormed))
# }
# 
# rlrNorm <- function(logDat) {
#   rlrNormed <- NormalyzerDE::performGlobalRLRNormalization(as.matrix(logDat), noLogTransform=TRUE)
#   colnames(rlrNormed) <- colnames(logDat)
#   row.names(rlrNormed) <- rownames(logDat)
#   return(as.matrix(rlrNormed))
# }
# 
# giNorm <- function(logDat) {
#   giNormed <- NormalyzerDE::globalIntensityNormalization(as.matrix(logDat), noLogTransform=TRUE)
#   colnames(giNormed) <- colnames(logDat)
#   row.names(giNormed) <- rownames(logDat)
#   return(as.matrix(giNormed))
# }