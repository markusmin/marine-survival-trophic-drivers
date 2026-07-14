# Create the objects for geometric anisotropy
# Original author: Jim Thorson

make_anisotropy_spde <-
  function( inla_mesh ){
    
    spde = fm_fem( inla_mesh )
    Dset <- 1:2
    TV <- inla_mesh$graph$tv
    V0 <- inla_mesh$loc[TV[, 1], Dset]
    V1 <- inla_mesh$loc[TV[, 2], Dset]
    V2 <- inla_mesh$loc[TV[, 3], Dset]
    E0 <- V2 - V1
    E1 <- V0 - V2
    E2 <- V1 - V0
    TmpFn <- function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area <- rep(NA, nrow(E0))
    for (i in seq_len(length(Tri_Area))){
      Tri_Area[i] <- TmpFn(E0[i,], E1[i,])/2
    }
    ret <- list( n_s = inla_mesh$n,
                 n_tri = nrow(TV),
                 Tri_Area = Tri_Area,
                 E0 = E0,
                 E1 = E1,
                 E2 = E2,
                 TV = TV - 1,
                 G0 = spde$c0,
                 G0_inv = as(Matrix::diag(1/Matrix::diag(spde$c0)),"TsparseMatrix") )
    return(ret)
  }