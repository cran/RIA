# Additional functions for geometry parameters
# (c): Márton Kolossváry, 2018

volume <- function(data, xy_dim = 0.5, z_dim = 0.4) {
  
  v_voxel <- xy_dim*xy_dim*z_dim
  
  data <- as.vector(data)
  data <- data[!is.na(data)]
  data[!is.na(data)] <- 1
  voxel_number <- sum(data, na.rm = T)
  
  return(voxel_number*v_voxel)
}

surface <- function(data, xy_dim, z_dim) {
  
  base_m <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  base_m[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- data
  
  xl <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  xl[1:(dim(data)[1]), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- data
  
  #xr <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  #xr[3:(dim(data)[1]+2), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- data
  
  yd <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  yd[2:(dim(data)[1]+1), 1:(dim(data)[2]), 2:(dim(data)[3]+1)] <- data
  
  #yu <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  #yu[2:(dim(data)[1]+1), 3:(dim(data)[2]+2), 2:(dim(data)[3]+1)] <- data
  
  zb <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  zb[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 1:(dim(data)[3])] <- data
  
  #zf <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  #zf[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 3:(dim(data)[3]+2)] <- data
  
  xl_m <- xl-base_m #; xr_m <- xr-base_m
  yd_m <- yd-base_m #; yu_m <- yu-base_m
  zb_m <- zb-base_m #; zf_m <- zf-base_m
  
  base_NA <- length(base_m[is.na(base_m)])
  xl_NA <- length(xl_m[is.na(xl_m)]); xl_NA_diff  <- xl_NA-base_NA
  #xr_NA <- length(xr_m[is.na(xr_m)]); xr_NA_diff  <- xr_NA-base_NA
  yd_NA <- length(yd_m[is.na(yd_m)]); yd_NA_diff  <- yd_NA-base_NA
  #yu_NA <- length(yu_m[is.na(yu_m)]); yu_NA_diff  <- yu_NA-base_NA
  zb_NA <- length(zb_m[is.na(zb_m)]); zb_NA_diff  <- zb_NA-base_NA
  #zf_NA <- length(zf_m[is.na(zf_m)]); zf_NA_diff  <- zf_NA-base_NA
  
  
  x_surf <- 2*xl_NA_diff*xy_dim*z_dim
  y_surf <- 2*yd_NA_diff*xy_dim*z_dim
  z_surf <- 2*zb_NA_diff*z_dim*z_dim
  
  surf <- x_surf+y_surf+z_surf
}


compactness1 <- function(vol, surf) {
  output <- ifelse( ( sqrt(pi) * (surf^(2/3)) )  > 0, vol/( sqrt(pi) * (surf^(2/3)) ), 0)
  return(output)
}

compactness2 <- function(vol, surf) {
  output <- ifelse( (surf^3)  > 0, 36*pi*((vol^2) / (surf^3)), 0)
  return(output)
}

spherical_dis <- function(vol, surf) {
  rs <- ((3*vol)/(4*pi))^(1/3)
  output <- ifelse( (4*pi*(rs^2)) > 0, surf/(4*pi*(rs^2)), 0)
  return(output)
}

sphericity <- function(vol, surf) {
  output <- ifelse( surf > 0, (((6*vol)^(2/3))*(pi^(1/3)))/surf, 0)
  return(output)
}


surface_dis <- function(data, xy_dim, z_dim) {
  
  base_m <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  base_m[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- data
  
  xl <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  xl[1:(dim(data)[1]), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- data
  
  xr <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  xr[3:(dim(data)[1]+2), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- data
  
  yd <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  yd[2:(dim(data)[1]+1), 1:(dim(data)[2]), 2:(dim(data)[3]+1)] <- data
  
  yu <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  yu[2:(dim(data)[1]+1), 3:(dim(data)[2]+2), 2:(dim(data)[3]+1)] <- data
  
  zb <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  zb[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 1:(dim(data)[3])] <- data
  
  zf <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  zf[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 3:(dim(data)[3]+2)] <- data
  
  xl_m <- xl-base_m ; xr_m <- xr-base_m
  yd_m <- yd-base_m ; yu_m <- yu-base_m
  zb_m <- zb-base_m ; zf_m <- zf-base_m
  
  null_m <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  
  null_m_xl <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  null_m_xl[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- xl_m[1:(dim(data)[1]), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)]
  null_m[(!is.na(base_m) & is.na(null_m_xl))] <- 1
  
  null_m_xr <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  null_m_xr[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- xr_m[3:(dim(data)[1]+2), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)]
  null_m[(!is.na(base_m) & is.na(null_m_xr))] <- 1
  
  null_m_yd <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  null_m_yd[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- yd_m[2:(dim(data)[1]+1), 1:(dim(data)[2]), 2:(dim(data)[3]+1)]
  null_m[(!is.na(base_m) & is.na(null_m_yd))] <- 1
  
  null_m_yu <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  null_m_yu[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- yu_m[2:(dim(data)[1]+1), 3:(dim(data)[2]+2), 2:(dim(data)[3]+1)]
  null_m[(!is.na(base_m) & is.na(null_m_yu))] <- 1
  
  null_m_zb <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  null_m_zb[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- zb_m[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 1:(dim(data)[3])]
  null_m[(!is.na(base_m) & is.na(null_m_zb))] <- 1
  
  null_m_zf <- array(NA, dim = c(dim(data)[1]+2, dim(data)[2]+2, dim(data)[3]+2))
  null_m_zf[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 2:(dim(data)[3]+1)] <- zf_m[2:(dim(data)[1]+1), 2:(dim(data)[2]+1), 3:(dim(data)[3]+2)]
  null_m[(!is.na(base_m) & is.na(null_m_zf))] <- 1
  
  null_m_df <- volume_to_dataframe(null_m, xy_dim = xy_dim, z_dim = z_dim)
  
  max_val <- 0
  for(i in 1: dim(null_m_df)[1]){
    for(j in i: dim(null_m_df)[1]) {
      
      dif <- null_m_df[i,1:3]-null_m_df[j,1:3]
      
      if(!is.na(null_m_df[i,4])&!is.na(null_m_df[j,4])) {
        dif <- dif^2
        dif <- sum(dif)
        dif <- sqrt(dif)
        
        if(dif > max_val) max_val <- dif
      }
    }
  }
  
  return(max_val)
  
}




# Claculates fractal dimenstion of dataset
#
# Creates gray-level co-occurrence matrix from \emph{RIA_image}, matrix or array.
#
# RIA_data_in \emph{RIA_image}, matrix or array
#
# return fractal dimensions


fractal <- function(RIA_data_in)
{
  data_in <-  RIA_data_in; if(all(is.na(data_in))) {
    output <- as.data.frame(matrix(c(0, 0, 0), nrow = 1, ncol = 3, byrow = T))
    colnames(output) <- c("Box_cunting_D", "Information_D", "Correlation_D")
    return(output)
  }
  
  if(length(dim(data_in)) != 3) {stop(paste0("Data is not 3D, but ", length(dim(data_in))," dimentional! Only 3D data is supported!"))}
  
  dim_x <- dim(data_in)[1]
  dim_y <- dim(data_in)[2]
  dim_z <- ifelse(!is.na(dim(data_in)[3]), dim(data_in)[3], 1)
  
  max_dim <- max(dim_x, dim_y, dim_z)
  size <- ceiling(logb(max_dim, 2))
  
  base_m <- array(NA, dim = c(2^size, 2^size, 2^size))
  base_m[1:dim_x, 1:dim_y, 1:dim_z] <- data_in
  base_m[!is.na(base_m)] <- 1
  r_dim <- array(NA, c(4, size+1))
  
  for(i in 0:size) {
    v <- 2^i
    iv <- 2^(size-i)
    scale_m <- array(NA, dim = c(v,v,v))
    
    for(x in 0:(v-1)) {
      for(y in 0:(v-1)) {
        for(z in 0:(v-1)) {
          
          scale_m[x+1,y+1,z+1] <- sum(base_m[(1+x*iv):((x+1)*iv), (1+y*iv):((y+1)*iv), (1+z*iv):((z+1)*iv)], na.rm = T)/(iv^3)
        }
      }
    }
    
    scale_m <- scale_m/sum(scale_m, na.rm = T)
    
    scale_m[scale_m == 0] <- NA
    r_dim[1,i+1] <- 2^i
    r_dim[2,i+1] <- sum(apply(scale_m, 2, function(x) sum(x!=0, na.rm = T)))
    r_dim[3,i+1] <- sum(ifelse(scale_m != 0, -1*scale_m*logb(scale_m, 2), 0), na.rm = T)
    r_dim[4,i+1] <- sum(scale_m^2, na.rm = T)
    
  }
  
  r_dim[1,] <- logb(r_dim[1,], 2); r_dim[2,] <- logb(r_dim[2,], 2); r_dim[4,] <- logb(r_dim[4,], 2)
  fit1 <- stats::lm(r_dim[2,] ~ r_dim[1,])
  fit2 <- stats::lm(r_dim[3,] ~ r_dim[1,])
  fit3 <- stats::lm(r_dim[4,] ~ r_dim[1,])
  
  output <- as.data.frame(matrix(c(fit1$coefficients[2], fit2$coefficients[2], -1*fit3$coefficients[2]), nrow = 1, ncol = 3, byrow = T))
  colnames(output) <- c("Box_cunting_D", "Information_D", "Correlation_D")
  
  return(output)
  
}

