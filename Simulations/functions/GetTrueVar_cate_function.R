#' Function that estimates the true variance.
#' 
#' Using draws of potential outcomes under the intervention and corresponding
#' estimates, we approximate the true variability of our estimator.
#' 
#' @param estimates List with elements corresponding to time points from the
#' start of the intervention until the end. Each element is a list in itself
#' including estimates according to our method for different draws according
#' to the intervention.
#' @param point.basis The new points projected on the selected spline basis

GetTrueVar <- function(estimates,weights, points.basis,E.basis,t.basis, separate_model=TRUE,estimand_beta,npixel = 128) {
  
  dims <- dim(estimates)
  num_interv <- dims[1]
  time_points <- dims[2]
  p <- dims[4]
  N <- nrow(points.basis)
  
  if(separate_model == FALSE){
    keep_index <- 1:((time_points-M+1)*npixel^2)
    Z <- na.omit((cbind(1,E.basis)%.%cbind(1,t.basis))[keep_index,])
    Q <- solve(t(Z)%*%Z)
  }
  
  mom_tau <- array(NA, dim = c(num_interv, num_interv,  N, 2, time_points))
  dimnames(mom_tau) <- list(interv1 = 1 : num_interv, interv2 = 1 : num_interv,
                            point = 1 : N, stat = c('var', 'mom2'),
                            time = 1 : time_points)
  
  mom_tau_haj <- array(NA, dim = c(num_interv, num_interv,  N, 2, time_points))
  dimnames(mom_tau_haj) <- list(interv1 = 1 : num_interv, interv2 = 1 : num_interv,
                            point = 1 : N, stat = c('var', 'mom2'),
                            time = 1 : time_points)
  
  tmp_mean <- array(NA,dim = c(num_interv, num_interv,  N, time_points))
  tmp_V <- array(NA, c(2*p+2,2*p+2, time_points))
  for (ii in 1 : (num_interv - 1)) {
    for (jj in (ii + 1) : num_interv) {

      for (tt in 1 : time_points) {
        
          # -------------------For IPW estimator ----------------------------
          if(separate_model==FALSE){
            var_beta <- var((estimates[jj, tt, , ]-estimates[ii, tt, , ])%*%t(Q))
            smom_beta <- array(rowMeans(apply((estimates[jj, tt, , ]-estimates[ii, tt, , ])%*%t(Q), 1, function(x) x%*%t(x))),c(p,p))
            mom_tau[ii,jj, , 1, tt] <- diag(points.basis%*%var_beta%*%t(points.basis))
            mom_tau[ii,jj, , 2, tt] <- diag(points.basis%*%smom_beta%*%t(points.basis))
          }else{
            tmp_mean[ii,jj,,tt] <- points.basis%*%colMeans(estimates[jj,tt, , ]-estimates[ii,tt, , ])
            var_beta <- var(estimates[jj, tt, , ]-estimates[ii, tt, , ])
            smom_beta <- array(rowMeans(apply(estimates[jj, tt, , ]-estimates[ii, tt, , ], 1, function(x) x%*%t(x))),c(p,p))
            mom_tau[ii,jj, , 1, tt] <- diag(points.basis%*%var_beta%*%t(points.basis))
            mom_tau[ii,jj, , 2, tt] <- diag(points.basis%*%smom_beta%*%t(points.basis))
          }

        
          
          #-----------------------For Hajek estimator--------------------------
          if(separate_model==FALSE){
            J <- cbind(-diag(p),
                       diag(p),
                       apply(po_for_var$estimates[ii, , , ], 3, function(x) mean(x,na.rm = TRUE)),
                       -apply(po_for_var$estimates[jj, , , ], 3, function(x) mean(x,na.rm = TRUE)))
            tmp_t <- cbind(estimates[ii, tt, , ],estimates[jj, tt, , ],weights[ii,tt,],weights[jj,tt,])

            var_beta <- Q%*%J%*%var(tmp_t)%*%t(J)%*%t(Q)
            smom_beta <- Q%*%J%*%array(rowMeans(apply(tmp_t, 1, function(x)x%*%t(x))),c(2*p+2,2*p+2))%*%t(J)%*%t(Q)
            mom_tau_haj[ii,jj, , 1, tt] <- diag(points.basis%*%var_beta%*%t(points.basis))
            mom_tau_haj[ii,jj, , 2, tt] <- diag(points.basis%*%smom_beta%*%t(points.basis))
            
          }else{
            
            J <- cbind(-diag(p),
                       diag(p),
                       estimand_beta[ii,],
                       -estimand_beta[jj,])
            
            # J <- cbind(-diag(p),
            #            diag(p),
            #            apply(po_for_var$estimates[ii, , , ], 3, function(x) mean(x,na.rm = TRUE)),
            #            -apply(po_for_var$estimates[jj, , , ], 3, function(x) mean(x,na.rm = TRUE)))
            #apply(po_for_var$estimates[1, , , ], 3, function(x) mean(x,na.rm = TRUE))
            tmp_t <- cbind(estimates[ii, tt, , ],estimates[jj, tt, , ],weights[ii,tt,],weights[jj,tt,])
            if(ii==1 & jj==3){
              tmp_V[,,tt] <- array(rowMeans(apply(tmp_t, 1, function(x)x%*%t(x))),c(2*p+2,2*p+2))
            }
            
            var_beta <- J%*%var(tmp_t)%*%t(J)
            smom_beta <- J%*%array(rowMeans(apply(tmp_t, 1, function(x)x%*%t(x))),c(2*p+2,2*p+2))%*%t(J)
            mom_tau_haj[ii,jj, , 1, tt] <- diag(points.basis%*%var_beta%*%t(points.basis))
            mom_tau_haj[ii,jj, , 2, tt] <- diag(points.basis%*%smom_beta%*%t(points.basis))
          }

          
        }

      }
    }
  
 
  return(list(moments_tau = mom_tau,moments_tau_haj = mom_tau_haj,tmp_mean=tmp_mean,tmp_V = rowMeans(tmp_V,dims = 2,na.rm = TRUE)))
  
}
