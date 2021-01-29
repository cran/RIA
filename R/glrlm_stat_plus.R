# statistics for glrlms
# (c): Márton Kolossváry, 2018

sre <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w <- col(data); w <- w^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(data/w)/sum_data, 0)
  return(output)
}

lre <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w <- col(data); w <- w^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(data*w)/sum_data, 0)
  return(output)
}

gln <- function(data) {
  if(length(data) == 0) {return (0)}
  
  row_s <- rowSums(data); row_s <- row_s^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(row_s)/sum_data, 0)
  return(output)
}

rln <- function(data) {
  if(length(data) == 0) {return (0)}
  
  row_s <- colSums(data); row_s <- row_s^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(row_s)/sum_data, 0)
  return(output)
}

rp <- function(data) {
  if(length(data) == 0) {return (0)}
  
  sum_data <- sum(data, na.rm = T)
  all_data <- sum(col(data)*data)
  output <- ifelse( all_data != 0, sum_data/all_data, 0)
  return(output)
}

lglre <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w <- row(data); w <- w^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(data/w)/sum_data, 0)
  return(output)
}

hglre <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w <- row(data); w <- w^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(data*w)/sum_data, 0)
  return(output)
}

srlgle <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w1 <- row(data); w1 <- w1^2
  w2 <- col(data); w2 <- w2^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(data/(w1*w2))/sum_data, 0)
  return(output)
}

lrhgle <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w1 <- row(data); w1 <- w1^2
  w2 <- col(data); w2 <- w2^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum(data*(w1*w2))/sum_data, 0)
  return(output)
}

srhgle <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w1 <- row(data); w1 <- w1^2
  w2 <- col(data); w2 <- w2^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum((data*w1)/(w2))/sum_data, 0)
  return(output)
}

lrlgle <- function(data) {
  if(length(data) == 0) {return (0)}
  
  w1 <- row(data); w1 <- w1^2
  w2 <- col(data); w2 <- w2^2
  sum_data <- sum(data, na.rm = T)
  output <- ifelse( sum_data != 0, sum((data*w2)/(w1))/sum_data, 0)
  return(output)
}


