# Additional first-order statistics
# (c): Márton Kolossváry, 2018

mode <- function(data) {
  uniq <- unique(data)
  tabeb <- tabulate(match(data, uniq))
  uniq[tabeb == max(tabeb)]
}

geo_mean <- function(data) {
  data_boo <- data<0
  data_abs <- abs(data)
  data_abs[data_abs==0] <- 1
  data_log <- log2(data_abs)
  
  output <- exp(mean(data_log))
  return(output)
}

geo_mean2 <- function(data) {
  data_boo <- data<0
  data_abs <- abs(data)
  data_abs[data_abs==0] <- 1
  data_log <- log2(data_abs)
  data_log[data_boo] <- data_log[data_boo] * -1
  
  output <- exp(mean(data_log))
  return(output)
}

geo_mean3 <- function(data) {
  min_data <- min(data)
  dev <- abs(min_data)
  if(min_data <= 0) data <- data + dev + 1
  output <- exp(mean(log(data)))
  return(output)
}

har_mean <- function(data) {
  data[data==0] <- 1
  output <- ifelse(mean(1/data) != 0, 1/mean(1/data), 0)
  return(output)
}

mn_AD_mn<- function(data){
  output <- mean(abs(data-mean(data)))
  return(output)
}

mn_AD_md <- function(data){
  output <- mean(abs(data-stats::median(data)))
  return(output)
}

md_AD_mn<- function(data){
  output <- stats::median(abs(data-mean(data)))
  return(output)
}

md_AD_md <- function(data){
  output <- stats::median(abs(data-stats::median(data)))
  return(output)
}

max_AD_mn<- function(data){
  output <- max(abs(data-mean(data)))
  return(output)
}
max_AD_md <- function(data){
  output <- max(abs(data-stats::median(data)))
  return(output)
}

skew <- function(data) {
  avg <- mean(data)
  SD  <- ifelse(length(data)>1, stats::sd(data), 0)
  output <- ifelse((SD != 0 | is.na(SD)), mean(((data-avg)^3))/(SD)^3, 0)
  return(output)
}

kurtosis <- function(data) {
  avg <- mean(data)
  SD  <- ifelse(length(data)>1, stats::sd(data), 0)
  output <- ifelse((SD != 0 | is.na(SD)), mean(((data-avg)^4))/(SD)^4, 0)
  return(output-3)
}

energy <- function(data) {
  output <- sum(data^2)
  return(output)
}

rms <- function(data) {
  output <- ifelse(length(data) != 0, sqrt(sum(data^2)/length(data)), 0)
  return(output)
}

uniformity <- function(data) {
  output <- ifelse(length(data) != 0, sum((table(data)/length(data))^2), 0)
  return(output)
}

entropy <- function (data, base = 2)
{
  p <- ifelse(length(data) != 0, table(data)/length(data), 0)
  l <- ifelse(p > 0, logb(p, base), 0)
  H <- sum(p * l)*-1
  return(H)
}
