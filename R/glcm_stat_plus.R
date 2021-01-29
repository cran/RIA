# Additional statistics on NxN symmetric matrixes
# (c): Márton Kolossváry, 2018

contrast <- function(data, type_in = "single", base = 2) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- (abs(row(ind_m)-col(ind_m)))^2
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  return(sum((ind_m)*data))
}

homogeneity2 <- function(data, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- 1/(((abs(row(ind_m)-col(ind_m)))^2)+1)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
}

dissimilarity <- function(data, type_in = "single", base = 2) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- (abs(row(ind_m)-col(ind_m)))
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  
  return(sum((ind_m)*data))
}


homogeneity1 <- function(data, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- 1/((abs(row(ind_m)-col(ind_m)))+1)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
}

dmn <- function(data, type_in = "single", base = 2) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- ((abs(row(ind_m)-col(ind_m)))^2)/dim_m^2
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  return(sum((ind_m)*data))
}

dn <- function(data, type_in = "single", base = 2) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- (abs(row(ind_m)-col(ind_m)))/dim_m
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  return(sum((ind_m)*data))
}

idmn <- function(data, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- 1/((((abs(row(ind_m)-col(ind_m)))^2)/dim_m^2)+1)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
}

idn <- function(data, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  
  ind_m <- ifelse(((((abs(row(ind_m)-col(ind_m))))/dim_m)+1) != 0, 1/((((abs(row(ind_m)-col(ind_m))))/dim_m)+1), 0)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
}

autocorrelation <- function(data, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- row(ind_m)*col(ind_m)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
}

inv_autocorrelation <- function(data, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  
  ind_m <- ifelse((row(ind_m)*col(ind_m)) != 0, 1/(row(ind_m)*col(ind_m)), 0)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
}

gauss <- function(data, inverse = FALSE, type_in = "single", base = 2, diag = TRUE, loc = 3) {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  
  if(loc == 1) {
    mu <- 1 #left polar
  } else if(loc == 2) {
    til  <- ceiling(dim_m/2)
    mu  <- mean(1:til, na.rm = TRUE) #left focus
  } else if(loc == 4) {
    from <- floor(dim_m/2)+1
    mu  <- mean(from:length(row(data)[,1]), na.rm = TRUE) #right focus
  } else if(loc == 5) {
    mu  <- dim_m #right polar
  } else {
    mu  <- mean(row(data)[,1], na.rm = TRUE) #center
  }
  
  sig <- stats::sd(row(data)[,1], na.rm = TRUE)
  
  if(is.na(sig) | sig == 0) {return(0)}
  
  const <- 1#/(sqrt(2*pi)*sig)
  expon <- exp(-1*((1/(2*sig^2))*(row(data)[,1]- mu)^2))
  
  if(inverse) {x<- matrix(1/(const*expon), nrow = dim_m, ncol = dim_m, byrow = TRUE)
  } else x<- matrix(const*expon, nrow = dim_m, ncol = dim_m, byrow = TRUE)
  if(inverse) {y<- matrix(1/(const*expon), nrow = dim_m, ncol = dim_m, byrow = FALSE)
  } else y<- matrix(const*expon, nrow = dim_m, ncol = dim_m, byrow = FALSE)
  
  ind_m <- x*y
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
  
}


gauss2f <- function(data, inverse = FALSE, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  til  <- ceiling(dim_m/2)
  from <- floor(dim_m/2)+1
  mu1  <- base::mean(1:til, na.rm = TRUE)
  mu2  <- base::mean(from:length(row(data)[,1]), na.rm = TRUE)
  sig <- stats::sd(row(data)[,1], na.rm = TRUE)
  
  if(is.na(sig) | sig == 0) {return(0)}
  
  const <- 1#/(sqrt(2*pi)*sig)
  expon1 <- exp(-1*((1/(2*sig^2))*(row(data)[,1]- mu1)^2))
  expon2 <- exp(-1*((1/(2*sig^2))*(row(data)[,1]- mu2)^2))
  
  if(inverse) {x1<- matrix(1/(const*expon1), nrow = dim_m, ncol = dim_m, byrow = TRUE)
  x2<- matrix(1/(const*expon2), nrow = dim_m, ncol = dim_m, byrow = TRUE)
  } else {x1<- matrix(const*expon1, nrow = dim_m, ncol = dim_m, byrow = TRUE)
  x2<- matrix(const*expon2, nrow = dim_m, ncol = dim_m, byrow = TRUE)}
  if(inverse) {y1<- matrix(1/(const*expon1), nrow = dim_m, ncol = dim_m, byrow = FALSE)
  y2<- matrix(1/(const*expon2), nrow = dim_m, ncol = dim_m, byrow = FALSE)
  } else {y1<- matrix(const*expon1, nrow = dim_m, ncol = dim_m, byrow = FALSE)
  y2<- matrix(const*expon2, nrow = dim_m, ncol = dim_m, byrow = FALSE)}
  
  ind_m <- x1*y1 + x2*y2
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
  
}

gauss2p <- function(data, inverse = FALSE, type_in = "single", base = 2, diag = TRUE) {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  sig   <- stats::sd(row(data)[,1], na.rm = TRUE)
  
  if(is.na(sig) | sig == 0) {return(0)}
  
  const <- 1#/(sqrt(2*pi)*sig)
  expon1 <- base::exp(-1*((1/(2*sig^2))*(row(data)[,1]- 1)^2))
  expon2 <- base::exp(-1*((1/(2*sig^2))*(row(data)[,1]- dim_m)^2))
  
  if(inverse) {x1<- matrix(1/(const*expon1), nrow = dim_m, ncol = dim_m, byrow = TRUE)
  x2<- matrix(1/(const*expon2), nrow = dim_m, ncol = dim_m, byrow = TRUE)
  } else {x1<- matrix(const*expon1, nrow = dim_m, ncol = dim_m, byrow = TRUE)
  x2<- matrix(const*expon2, nrow = dim_m, ncol = dim_m, byrow = TRUE)}
  if(inverse) {y1<- matrix(1/(const*expon1), nrow = dim_m, ncol = dim_m, byrow = FALSE)
  y2<- matrix(1/(const*expon2), nrow = dim_m, ncol = dim_m, byrow = FALSE)
  } else {y1<- matrix(const*expon1, nrow = dim_m, ncol = dim_m, byrow = FALSE)
  y2<- matrix(const*expon2, nrow = dim_m, ncol = dim_m, byrow = FALSE)}
  
  ind_m <- x1*y1 + x2*y2
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(ind_m) <- 0}
  
  return(sum((ind_m)*data))
  
}

cluster <- function(data, pow = 4, base = 2, inverse = FALSE, type_in = "single", diag = TRUE) {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  
  ind_m <- row(ind_m)+col(ind_m)
  ind_mc <- row(ind_m)
  ind_mr <- col(ind_m)
  
  ind_mcp <- ind_mc * data
  ind_mrp <- ind_mr * data
  
  col_avg   <- base::colMeans(ind_mcp, na.rm = TRUE, dims = 1)
  col_avg_m <- matrix(col_avg, dim_m, dim_m, byrow = TRUE)
  row_avg   <- base::rowMeans(ind_mrp, na.rm = TRUE, dims = 1)
  row_avg_m <- matrix(row_avg, dim_m, dim_m, byrow = FALSE)
  
  w <- (abs(ind_m - col_avg_m - row_avg_m))^pow
  
  if(inverse) w <- ifelse(w != 0, 1/w, 0)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  if (diag == FALSE) {diag(w) <- 0}
  
  
  
  return(sum( w*data ))
}

avg <- function(data, base = 2, type_in = "single") {
  dim_m <- dim(data)[1]
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- row(ind_m)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  return(sum(ind_m * data))
}

variance <- function(data, base = 2, type_in = "single") {
  dim_m <- dim(data)[1]
  
  if(dim_m == 1 & data[1,1] == TRUE) return (0)
  
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- row(ind_m)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  mu <- sum(ind_m * data)
  
  func <- (ind_m - mu)^2
  
  return(sum(data*func))
}

correlation <- function(data, base = 2, type_in = "single") {
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_mr <- row(ind_m)
  ind_mc <- col(ind_m)
  
  if (type_in == "squared") {data <- data^2}
  if (type_in == "entropy") {data <- ifelse(data > 0, -1*data*logb(data, base), 0)}
  
  mu <- sum(ind_mr * data)
  func <- (ind_mr - mu)^2
  sig <- sum(data*func)
  
  if(is.na(sig) | sig == 0) {return(1)}
  
  return(sum(data*( ((ind_mr - mu) * (ind_mc - mu)) / sig)))
}

sum_f <- function(data, base = 2, inverse = FALSE, type_in = "single")
{
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- row(ind_m)+col(ind_m)
  
  sum_all <- 0
  
  if (type_in == "variance") sum_func_var <- sum_entropy(data)
  
  for(i in 2: (2*dim_m)) {
    sum_p <- sum(data[ind_m==i], na.rm = TRUE)
    
    if (type_in == "squared") sum_p <- sum_p^2
    
    if (type_in == "single" | type_in == "squared") sum_func <- i
    if (type_in == "entropy") sum_func <- ifelse(sum_p > 0, -1*logb(sum_p, base), 0)
    if (type_in == "variance") sum_func <- (i-sum_func_var)^2
    
    if(inverse == TRUE) sum_func <- 1/sum_func
    
    
    sum_all <- sum_all + (sum_p*sum_func)
  }
  
  return(sum_all)
  
}

dif_f <- function(data, base = 2, inverse = FALSE, type_in = "single")
{
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- base::abs(row(ind_m)-col(ind_m))
  
  sum_all <- 0
  
  if (type_in == "variance") sum_func_var <- dif_entropy(data)
  if (type_in == "squared") data <- data^2
  
  for(i in 0: (dim_m-1)) {
    sum_p <- sum(data[ind_m==i], na.rm = TRUE)
    if (type_in == "squared") sum_p <- sum_p^2
    
    if (type_in == "single" | type_in == "squared") sum_func <- i
    if (type_in == "entropy") sum_func <- ifelse(sum_p > 0, -1*logb(sum_p, base), 0)
    if (type_in == "variance") sum_func <- (i-sum_func_var)^2
    
    if(inverse == TRUE) sum_func <- ifelse(sum_func != 0, 1/sum_func, 0)
    
    
    sum_all <- sum_all + (sum_p*sum_func)
  }
  
  return(sum_all)
}

sum_entropy <- function(data, base = 2) {
  #helper function for sum_f
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- row(ind_m)+col(ind_m)
  
  sum_all <- 0
  for(i in 2: (2*dim_m)) {
    sum_p <- sum(data[ind_m==i], na.rm = TRUE)
    sum_func <- ifelse(sum_p > 0, logb(sum_p, base), 0)
    
    sum_all <- sum_all + (sum_p*sum_func*-1)
  }
  
  return(sum_all)
}

dif_entropy <- function(data, base = 2) {
  #helper function for sum_f
  dim_m <- dim(data)[1]
  ind_m <- matrix(NA, dim_m, dim_m)
  ind_m <- abs(row(ind_m)-col(ind_m))
  
  sum_all <- 0
  for(i in 0: (dim_m-1)) {
    sum_p <- sum(data[ind_m==i], na.rm = TRUE)
    sum_func <- ifelse(sum_p > 0, logb(sum_p, base), 0)
    
    
    sum_all <- sum_all + (sum_p*sum_func*-1)
  }
  
  return(sum_all)
}


energy_m <- function(data) {
  
  data <- data^2
  
  return(sum(data))
}

entropy_m <- function(data, base = 2) {
  
  data <- ifelse(data > 0, -1*data*logb(data, base), 0)
  
  return(sum(data))
}


imc1 <-  function(data, base = 2) {
  H <- entropy_m(data)
  
  col_avg <- colMeans(data, na.rm = TRUE, dims = 1)
  row_avg <- rowMeans(data, na.rm = TRUE, dims = 1)
  marg_avg <- col_avg%*%t(row_avg)
  l_marg_avg <- ifelse(marg_avg > 0, logb(marg_avg, base), 0)
  
  HXY1 <- sum(data*l_marg_avg)*-1
  
  lx <- ifelse(row_avg > 0, logb(row_avg, base), 0)
  HX <- sum(row_avg*lx)*-1
  
  ly <- ifelse(col_avg > 0, logb(col_avg, base), 0)
  HY <- sum(col_avg*ly)*-1
  
  if(max(HX, HY) == 0) {return(0)}
  
  return( (H-HXY1)/max(HX, HY))
  
}

imc2 <-  function(data, base = 2) {
  H <- entropy(data)
  
  col_avg <- colMeans(data, na.rm = TRUE, dims = 1)
  row_avg <- rowMeans(data, na.rm = TRUE, dims = 1)
  marg_avg <- col_avg%*%t(row_avg)
  l_marg_avg <- ifelse(marg_avg > 0, logb(marg_avg, base), 0)
  
  HXY2 <- sum(marg_avg*l_marg_avg)*-1
  
  return( sqrt(abs(1 - exp(-2*(HXY2-H)))) )
  
}

