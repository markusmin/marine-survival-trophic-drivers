H_matrix_prior_predictive_check_for_ggplot <- function(ln_H_1_mean, ln_H_1_sd,
                                                       ln_H_2_mean, ln_H_2_sd,
                                                       range){
  # density for parameter H_1
  ln_H_1 <- rnorm(1000, ln_H_1_mean, ln_H_1_sd)
  ln_H_2 <- rnorm(1000, ln_H_2_mean, ln_H_2_sd)
  H <- array(NA, dim = c(2, 2, 1000))
  
  Eigen <- list()
  
  # calculate H
  for (i in 1:1000){
    # H[1,1,i] = exp(ln_H_1[i])
    # H[2,1,i] = exp(ln_H_2[i])
    # H[1,2,i] = exp(ln_H_2[i])
    # H[2,2,i] = (1+exp(ln_H_2[i])*exp(ln_H_2[i])) / exp(ln_H_1[i])
    
    # don't exponentiate H_2, otherwise we are constraining it to be positive
    H[1,1,i] = exp(ln_H_1[i])
    H[2,1,i] = ln_H_2[i]
    H[1,2,i] = ln_H_2[i]
    H[2,2,i] = (1+ln_H_2[i]*ln_H_2[i]) / exp(ln_H_1[i])
    
    # don't exponentiate either?
    # this leads to some crazy anisotropy that doesn't make any sense
    # H[1,1,i] = ln_H_1[i]
    # H[2,1,i] = ln_H_2[i]
    # H[1,2,i] = ln_H_2[i]
    # H[2,2,i] = (1+ln_H_2[i]*ln_H_2[i]) / ln_H_1[i]
    
    # get eigen
    # get eigen decomposition
    Eigen[[i]] = eigen(H[,,i])
  }
  # convert to df
  eigen_values <- lapply(Eigen, function(x) x[[1]])
  eigen_vectors <- lapply(Eigen, function(x) as.vector(x[[2]]))
  
  eigen_values <- do.call(rbind, eigen_values)
  eigen_vectors <- do.call(rbind, eigen_vectors)
  
  eigen_df <- as.data.frame(cbind(eigen_values, eigen_vectors))
  colnames(eigen_df) <- c("value1", "value2", "vector1", "vector2", "vector3", "vector4")
  rownames(eigen_df) <- NULL
  
  # make plots
  # get max to set the range
  Major_1 = Minor_1 = Major_2 = Minor_2 = NA
  # use estimated range to get these values
  Major_1 = c(eigen_df$vector1, eigen_df$vector3)*eigen_df$value1 * range
  Minor_1 = c(eigen_df$vector2, eigen_df$vector4)*eigen_df$value2 * range
  Range = 1.1 * c(-1,1) * max(abs( cbind(Major_1,Minor_1) ),na.rm=TRUE)
  
  ellipse_plotting_df <- data.frame(semi_major = rep(NA, 1000), 
                                    semi_minor = rep(NA, 1000),
                                    angle = rep(NA, 1000))
  
  for (i in 1:1000){
    
    rss = function(V) sqrt(sum(V[1]^2+V[2]^2))
    Major_1 = Minor_1 = Major_2 = Minor_2 = NA
    # use estimated range to get these values
    Major_1 = Eigen[[i]]$vectors[,1]*Eigen[[i]]$values[1] * range
    Minor_1 = Eigen[[i]]$vectors[,2]*Eigen[[i]]$values[2] * range
    
    # plot the ellipse
    # shape::plotellipse( rx=rss(Major_1), ry=rss(Minor_1), angle=-1*(atan(Major_1[1]/Major_1[2])/(2*pi)*360-90), lcol="green", lty="solid")
    
    # store this ellipse in the df
    ellipse_plotting_df$semi_major[i] <- rss(Major_1)
    ellipse_plotting_df$semi_minor[i] <- rss(Minor_1)
    ellipse_plotting_df$angle[i] <- -1*(atan(Major_1[1]/Major_1[2])/(2*pi)*360-90)
    
    
  }
  return(ellipse_plotting_df)
}