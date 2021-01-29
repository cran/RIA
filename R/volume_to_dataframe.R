# Changes 3D array information to coordinate space matrix
#
# Changes 3D array to 2D matrix where first 2 or three columns are the coordinaltes
# and the last one is the voxel value.
#
# data 2-3D array
# xy_dim inplane resolution
# y_dim crossplane resolution
# return 3 or 4 column matrix with coordinates and voxel values
# (c): Márton Kolossváry, 2018


volume_to_dataframe <- function(data, xy_dim = 0.5, z_dim = 0.5) {
  dims <- dim(data)
  x_length <- 1:dims[1]
  y_length <- 1:dims[2]
  z_length <- 1:dims[3]
  
  null_mx <- array(NA, dims); null_my <- array(NA, dims); null_mz <- array(NA, dims)
  
  for(i in 1:dims[2]) {
    for(j in 1:dims[3]) {
      null_mx[,i,j] <- x_length
      
    }
  }
  
  for(i in 1:dims[1]) {
    for(j in 1:dims[3]) {
      null_my[i,,j] <- y_length
      
    }
  }
  
  for(i in 1:dims[1]) {
    for(j in 1:dims[2]) {
      null_mz[i,j,] <- z_length
      
    }
  }
  
  null_mx <-  null_mx*xy_dim;  null_my <-  null_my*xy_dim; null_mz<- null_mz*z_dim
  
  length_m <- dims[1]*dims[2]*dims[3]
  df <- matrix(NA, nrow = length_m, ncol = 4)
  
  p <- 1
  for(i in 1:dims[1]) {
    for(j in 1:dims[2]) {
      for(k in 1:dims[3]) {
        
        df[p, 1] <- null_mx[i,j,k]
        df[p, 2] <- null_my[i,j,k]
        df[p, 3] <- null_mz[i,j,k]
        df[p, 4] <- data[i,j,k]
        p <- p+1
      }
      
    }
    
  }
  
  return(df)
  
}
