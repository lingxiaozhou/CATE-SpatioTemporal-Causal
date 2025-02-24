# Date: 2/23/2025
# General description: Simulations over the Iraq map. For covariates I use
# distance from the road network, distance from cities, as well as point
# patterns that are generated using the observed treatment and outcome
# intensities.
# We perform the following:
# (1) generate data (for spatio effect modifier setting),
# (2) calculate estimand (since the estimand is data specific), and the projected estimand is considered
# (3) for the specific intervenion, we can calculate the estimator's
#     asymptotic variance, and the bound of the asymptotic variance we are
#     using.
# (4) Calculate the estimate based on IPW and Hajek
# (5) Calculate the estimate of the variance bound
# (6) Perform hypothesis test for no heterogeneity effect



library(spatstat)
library(spatstat.geom)
library(splines)
library(mgcv)
library(abind)  # For combine 3D arrays


# -----------------------Step 1: Specifications----------------------------------

wh_sims <- 14       # name of the simulation

index <- 2          # index of the simulation
time_points <- 200 # number of total time periods
type_interv <- 1   # length of the intervention: 1 for single time period, 2 for 3 periods and 5 for 7 time periods
type_em <- 1       # type of the effect modifier: 1 for spatial and 3 for spatio-temporal 
separate_model <- 1   # Fit separate regression for each time period (always set to 1)


interaction <- type_em


# for spatial effect modifier (em_type=1) use 8,otherwise use 11
if(type_em==1){
  wh_data_specs <- 1
}

if(type_em==3){
  wh_data_specs <- 2
}


   
wh_interventions <- 1


nbaseE <- 6   # number of basis for the splines
nbaset <- 4

points <- 1:100/100   # on which we evaluate the cate
npixel <- 128         # Number of npixel cells per dimension (space is divided into npixel x npixel)



#--------------------------------specify the path and which task will be done------------------

separate_path <- ifelse(separate_model,"/separate","/unified")
out_path <- paste0('~/Github/CATE-SpatioTemporal-Causal/Simulations/Output/', wh_sims,
                    '_sims/Results',separate_path,'/T', time_points, '_type', type_interv,'_em',type_em, '/')

load_data_path <- '~/Github/CATE-SpatioTemporal-Causal/Simulations/Iraq_info/'
source_path <- '~/Github/CATE-SpatioTemporal-Causal/Simulations/functions/'


perform_estimand <- TRUE
perform_estimation <- TRUE


#--------------------------load functions and data--------------------------------------

file.sources <- list.files(source_path, pattern="*.R$")
sapply(paste0(source_path, file.sources), source, .GlobalEnv)



load(paste0('~/Github/CATE-SpatioTemporal-Causal/Simulations/Data_specs/data_specs_', wh_data_specs, '.dat'))
load(paste0('~/Github/CATE-SpatioTemporal-Causal/Simulations/Data_specs/interventions', wh_interventions, '.dat'))
load(paste0(load_data_path, 'cities_dist.dat'))
load(paste0(load_data_path, 'routes_dist.dat'))


Xt_int <- data_specs$Xt_int
Xt_int2 <- data_specs$Xt_int2
trt_coefs <- data_specs$trt_coefs
trt_lags <- data_specs$trt_lags
out_coefs <- data_specs$out_coefs
out_lags <- data_specs$out_lags
confounder_priority <- ifelse('confounder_priority' %in% names(data_specs),
                              data_specs$confounder_priority, 2)



# ---------- LOADING THE DATA WE NEED FOR IRAQ SPECIFIC SIMULATIONS -------- 

load(paste0(load_data_path, 'baghdad_window.dat'))
load(paste0(load_data_path, 'mosul_window.dat'))
load(paste0(load_data_path, 'mosul_window.dat'))
load(paste0(load_data_path, 'obs_dens.dat'))


# Creating the priority column that will be used:
obs_dens$Priority <- with(obs_dens, Density)
if ('interv_cov' %in% names(data_specs)) {
  interv_cov <- data_specs$interv_cov  # this is phi(w)
} else {
  interv_cov <- log(obs_dens$Priority[[1]] + 1e-5)
}


# -------- CREATING IRAQ SPECIFIC DATA -------- 

# The windows:
# Using the Iraq window from the roads because it's rougher hence smaller:

window <- Window(routes_dist)
B <- window

# Making the time invariant covariates
# One of them is the Iraq road network and the other is distance from Baghdad

if(npixel!=128){
  routes_dist <- as.im(routes_dist,dimyx = c(npixel,npixel))
  cities_dist$distance[1][[1]] <- as.im(cities_dist$distance[1][[1]],dimyx = c(npixel,npixel))
  interv_cov <- as.im(interv_cov,dimyx = c(npixel,npixel))
}



roads <- hyperframe(Dist = list(routes_dist, cities_dist$distance[1][[1]]),
                    row.names = c('routes', 'baghdad'))
roads$Priority <- with(roads, exp(- Dist))
roads$Priority[1][[1]] <- roads$Priority[1][[1]] ^ 3

# roads$Priority[1][[1]] is the X^1 in the paper
roads$Priority[1][[1]] <- roads$Priority[1][[1]] + log(bdist.pixels(window,dimyx = c(npixel,npixel))) # bdist.pixel compute the distance to the border





ps_coefs <- trt_coefs
ps_lags <- trt_lags
intervention <- interventions[[type_interv]]
if (length(dim(intervention)) == 2) {
  intervention <- array(intervention, dim = c(dim(intervention), 1))
}

M <- dim(intervention)[3]   # length of the intervention
reps <- 40
spat_kernel <- 'gaussian'
spat_adjust <- 10 / (time_points ^ (2 / 3)) * sqrt(max(diff(window$xrange), diff(window$yrange)))
edge_cor <- TRUE
truncation_level <- 0.99

if (index == 1) {
  specs <- list(wh_data_specs = wh_data_specs, trt_coefs = trt_coefs,
                out_coefs = out_coefs, Xt_int = Xt_int, Xt_int2 = Xt_int2,
                trt_lags = trt_lags, out_lags = out_lags, ps_coefs = ps_coefs,
                ps_lags = ps_lags, intervention = intervention,
                window = window, B = B, reps = reps, spat_kernel = spat_kernel,
                spat_adjust = spat_adjust, edge_cor = edge_cor,
                truncation_level = truncation_level)
  save(specs, file = paste0(out_path, 'specs.dat'))
}


# ----------- STEP 2: Generating the data ----------- 

set.seed(index)

# The first covariate will be most predictive of treatments, while the second
# most predictive of outcomes.
# Xt is X^3 in the paper
Xt <- GenConf(time_points = time_points, roads = obs_dens[1, ], Xt_int = Xt_int,
              priority_coef = confounder_priority,dimyx = c(npixel = npixel))
# Xt2 is X^4 in the paper
Xt2 <- GenConf(time_points = time_points, roads = obs_dens[2, ], Xt_int = Xt_int2,
               priority_coef = confounder_priority,dimyx = c(npixel = npixel))
# plot(Xt2$Priority[[2]])

dta <- GenTrtOut(roads = roads, Xt = Xt, Xt2 = Xt2, trt_coefs = trt_coefs,
                 trt_lags = trt_lags, out_coefs = out_coefs,
                 out_lags = out_lags,interaction = interaction,dimyx=c(npixel,npixel),linear=linear)
Wt <- dta$Wt
Yt <- dta$Yt
mean(sapply(1:time_points, function(x) Yt$Points[[x]]$n))
mean(sapply(1:time_points, function(x) Wt$Points[[x]]$n))
# dta$YtPriority <- GetPriority(Yt)$Priority
# get the previous lag 1 outcome priority
if(type_em == 2){
  dta$prev_out <- lapply(1:time_points, function(x) GetHistory(tt = x, Wt = Wt$Points, Yt = Yt$Points,
                                                               center_variable = 'Y', lags = c(1,1),
                                                               window = window,priority_coefs=c(2,2))$Yt_hist$Priority[[1]])
}



na.matrix <- ifelse(is.na(as.matrix(interv_cov)),NA,1)


# Get the effect modifiers in the matrix form
if(type_em==1){
  E.mat <- imls_to_mat(roads$Priority[2][[1]],1,1,na.matrix)
  E.mat <- array(rep(c(E.mat),time_points), c(dim(E.mat)[1:2],time_points))
}


if(type_em==2){
  E.mat <- imls_to_mat(dta$prev_out,1,time_points,na.matrix) 

}


if(type_em==3){
  E.mat <- imls_to_mat(Xt2$Priority,1,time_points,na.matrix) 
  
}




rm(dta)
na.index <- is.na(E.mat[,,1])
if(type_em!=4){
  # E.basis <- ns(na.omit(c(E.mat)),df = nbaseE-1, Boundary.knots=c(0,1),knots = specify_knots(nbaseE-2))
  E.basis <- ns(na.omit(c(E.mat)),df = nbaseE-1, Boundary.knots=c(0,1))
  if(separate_model==TRUE){
    points.basis <- cbind(1,predict(E.basis,newx = points))
  }else{
    t.basis <- as.matrix(ns(sort(rep(1:time_points,prod(dim(E.mat)[1:2]))),df = nbaset-1))
    points.basis <- cbind(1,predict(E.basis,newx = points))%.%matrix(rep(c(1,colMeans(t.basis)),nrow(points.basis)),
                                                                     ncol=nbaset,
                                                                     byrow=T)
  }
}





# ---------- STEP 3: Calculating the estimand ----------- 
# Generating potential outcomes under the intervention:
if (perform_estimand) {
  
  total_pixel <- 8451   # total number of pixels in iraq when use 128 X 128 grid
  counter <- GenPotOut(intervention = intervention,
                       interv_cov = interv_cov,
                       roads = roads, Xt = Xt, Xt2 = Xt2, Wt = Wt, Yt = Yt,
                       out_coefs = out_coefs, out_lags = out_lags,
                       reps = reps, short = FALSE,separate_prev_trt = (interaction!=1),interaction = interaction,linear = FALSE)
  
  
  Ap <- integral(counter$av_surface[[1]])/sum(counter$av_surface[[1]])  # area of a pixel
  
  
  estimand <- array(NA,c(length(points),3,dim(intervention)[1]))
  estimand_beta <- array(NA,c(dim(intervention)[1], ncol(points.basis)))
  
  dimnames(estimand)[[2]] <- c("proj","true","meaneffect")
  for(this_interv in 1: dim(intervention)[1]){
    tau_t_proj <- array(NA,c(nrow(points.basis),time_points-M+1))
    tau_t <- array(NA,c(nrow(points.basis),time_points-M+1))
    meaneffect <- c()
    for (tt in M:time_points-M+1) {
      index_tmp <- (1:(total_pixel))+(tt-1)*total_pixel
      
      if(type_em!=4){
        df_tmp <- cbind(1,E.basis[index_tmp,], Y = na.omit(Ap*c(as.matrix(counter$surfaces[[this_interv]][[tt]]))), E = na.omit(c(E.mat[,,tt])))
      }else{
        df_tmp <- cbind(1,E.basis[index_tmp,], Y = na.omit(Ap*c(as.matrix(counter$surfaces[[this_interv]][[tt]]))), E = na.omit(E.basis[index_tmp,]))
      }
      
      df_tmp <- na.omit(df_tmp)
      
      
      
      tau_t[,tt] <- ksmooth(df_tmp[,nbaseE+2],df_tmp[,nbaseE+1],bandwidth=0.15,x.points = points)$y
      meaneffect <- c(meaneffect, mean(df_tmp[,nbaseE+1]))
      
      
    }
    surface <- imls_to_mat(counter$surfaces[[this_interv]],start = 1,end = time_points-M+1,na.matrix = na.matrix)
    res_estimand <- get_cate1(0, surface,E.basis,t.basis,point.basis,
                              used_time=(M):time_points-M+1, C=1/Ap,
                              Hajek = FALSE,weights = NULL, separate_model = separate_model)
    estimand[,1,this_interv] <- res_estimand$est
    estimand[,2,this_interv] <- rowMeans(tau_t,na.rm = TRUE)
    estimand[,3,this_interv] <- mean(meaneffect)
    estimand_beta[this_interv,] <- res_estimand$beta
    rownames(estimand_beta) <- c("c=3","c=5","c=7")
  }
  
  save(estimand,estimand_beta, file = paste0(out_path, 'estimand_', index, '.dat'))
  rm(counter)
  rm(surface)
  
}

# ---------- STEP 4: Getting the estimate (CATE, Variance bound, Hypothesis test) ----------- 

# ----- PART A: Without truncation, based on the true PS.

if (perform_estimation) {
  
  

  ps_coefs <- trt_coefs
  estimates <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2, Wt = Wt, Yt = Yt,
                             ps_coefs = ps_coefs, ps_lags = trt_lags,
                             intervention = intervention,
                             interv_cov = interv_cov,
                             truncation_level = NULL,
                             spat_kernel = spat_kernel,
                             spat_adjust = spat_adjust, edge_cor = edge_cor,
                             return_smoothed_outcome = TRUE)
  if(M>1){
    intervention_M1 <- interventions[[1]]
    
    intervention_M1 <- array(intervention_M1, dim = c(dim(intervention_M1), 1))
  
    weights_M1 <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2, Wt = Wt, Yt = Yt,
                               ps_coefs = ps_coefs, ps_lags = trt_lags,
                               intervention = intervention_M1,
                               interv_cov = interv_cov,
                               truncation_level = NULL,
                               spat_kernel = spat_kernel,
                               spat_adjust = spat_adjust, edge_cor = edge_cor,
                               return_smoothed_outcome = TRUE)$weights
  }else{
    weights_M1 <- estimates$weights
  }

  
  smoothed_outcome <- estimates$smoothed_outcome
  
  used_time <- as.numeric(names(estimates$weights[1,]))-M+1
  
  for (interv1 in 1:dim(intervention)[1]) {
    for (interv2 in 1:dim(intervention)[1]){
      if(interv1<interv2){
        cat(interv1,interv2)
        
        cate <- get_cate1(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                         E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = FALSE,
                         weights = estimates$weights[c(interv1,interv2),], separate_model=separate_model)
        cate_haj <- get_cate2(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                              E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = TRUE,
                              weights = estimates$weights[c(interv1,interv2),],M=M, separate_model=separate_model,weights_M1=weights_M1[c(interv1,interv2),])
        
        
        save(cate, cate_haj, file = paste0(out_path, 'cate_',interv1,"_",interv2,"_", index, '.dat'))
      }
      
    }
  }
  
  
  rm(estimates)
  
}



# ----- PART B: With truncation, based on the true PS.


if (perform_estimation) {
  
 
  estimates <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2,
                             Wt = Wt, Yt = Yt,
                             ps_coefs = trt_coefs, ps_lags = trt_lags,
                             intervention = intervention,
                             interv_cov = interv_cov,
                             truncation_level = truncation_level,
                             spat_kernel = spat_kernel,
                             spat_adjust = spat_adjust, edge_cor = edge_cor,
                             smoothed_outcome = smoothed_outcome)
  

  
  save(estimates, file = paste0(out_path, 'estimate_trun_', index, '.dat'))
  used_time <- as.numeric(names(estimates$weights[1,]))-M+1
  
  for (interv1 in 1:dim(intervention)[1]) {
    for (interv2 in 1:dim(intervention)[1]){
      if(interv1<interv2){
        cat(interv1,interv2)
        cate <- get_cate1(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                         E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = FALSE,
                         weights = estimates$weights[c(interv1,interv2),], separate_model=separate_model)
        cate_haj <- get_cate2(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                              E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = TRUE,
                              weights = estimates$weights[c(interv1,interv2),],separate_model=separate_model,M=M,weights_M1 = weights_M1[c(interv1,interv2),])
        

        
        save(cate, cate_haj, file = paste0(out_path, 'cate_trun_',interv1,"_",interv2,"_", index, '.dat'))
      }
      
    }
  }
  
  rm(estimates)
  
}





# -------- PART D: Using the estimated PS ----- #


# Estimating the propensity score.
if (perform_estimation) {
  
  start_time <- max(2, dim(intervention)[1])
  cat(start_time,"and",time_points,"\n")
  cat(length(Wt$Points[start_time : time_points]),",",length(Xt$Priority[start_time : time_points]),",",length(Xt2$Priority[start_time : time_points]),"\n")
  ps_dta <- hyperframe(Points_Wt = Wt$Points[start_time : time_points],
                       road_priority = roads$Priority[[1]],
                       road2_priority = roads$Priority[[2]],
                       Xt_priority = Xt$Priority[start_time : time_points],
                       Xt2_priority = Xt2$Priority[start_time : time_points],
                       Wt_prev = NULL, Yt_prev = NULL)
  cat("finish ps_dta")
  for (tt in start_time : time_points) {
    hist <- GetHistory(tt = tt, Wt = Wt$Points, Yt = Yt$Points,
                       center_variable = 'W', lags = ps_lags,
                       window = window)
    ps_dta$Wt_prev[[tt - start_time + 1]] <- hist$Wt_hist$Priority[[1]]
    ps_dta$Yt_prev[[tt - start_time + 1]] <- hist$Yt_hist$Priority[[1]]
  }
  
  
  est_ps_mod <- mppm(Points_Wt ~ road_priority + Xt_priority + Wt_prev +
                       Yt_prev + road2_priority + Xt2_priority,
                     data = ps_dta)
  est_ps_coef <- as.numeric(summary(est_ps_mod)$coef)
  save(est_ps_coef, file = paste0(out_path, 'est_ps_coef_', index,  '.dat'))
  
}



if (perform_estimation) {
  
  # Without truncation:
  estimates <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2, Wt = Wt, Yt = Yt,
                             ps_coefs = est_ps_coef, ps_lags = trt_lags,
                             intervention = intervention,
                             interv_cov = interv_cov,
                             truncation_level = NULL,
                             spat_kernel = spat_kernel,
                             spat_adjust = spat_adjust, edge_cor = edge_cor,
                             smoothed_outcome = smoothed_outcome)
  
  used_time <- as.numeric(names(estimates$weights[1,]))-M+1
  
  for (interv1 in 1:dim(intervention)[1]) {
    for (interv2 in 1:dim(intervention)[1]){
      if(interv1<interv2){
        cat(interv1,interv2)
        cate <- get_cate1(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                         E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = FALSE,
                         weights = estimates$weights[c(interv1,interv2),], separate_model=separate_model)
        cate_haj <- get_cate2(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                              E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = TRUE,
                              weights = estimates$weights[c(interv1,interv2),],separate_model=separate_model,M=M,weights_M1 = weights_M1[c(interv1,interv2),])
        
        save(cate, cate_haj, file = paste0(out_path, 'cate_estPS_',interv1,"_",interv2,"_", index, '.dat'))
      }
      
    }
  }
  

  
  
  
  # With truncation:
  estimates <- CalcEstimates(roads = roads, Xt = Xt, Xt2 = Xt2, Wt = Wt, Yt = Yt,
                             ps_coefs = est_ps_coef, ps_lags = trt_lags,
                             intervention = intervention,
                             interv_cov = interv_cov,
                             truncation_level = truncation_level,
                             spat_kernel = spat_kernel,
                             spat_adjust = spat_adjust, edge_cor = edge_cor,
                             smoothed_outcome = smoothed_outcome)
  

  used_time <- as.numeric(names(estimates$weights[1,]))-M+1
  for (interv1 in 1:dim(intervention)[1]) {
    for (interv2 in 1:dim(intervention)[1]){
      if(interv1<interv2){
        cat(interv1,interv2)
        cate <- get_cate1(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                         E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = FALSE,
                         weights = estimates$weights[c(interv1,interv2),], separate_model=separate_model)
        cate_haj <- get_cate2(estimates$surfaces[[interv1]], estimates$surfaces[[interv2]],
                              E.basis,t.basis,point.basis,used_time=used_time, C=1,Hajek = TRUE,
                              weights = estimates$weights[c(interv1,interv2),], separate_model=separate_model,M=M,weights_M1 = weights_M1[c(interv1,interv2),])
        

        save(cate, cate_haj, file = paste0(out_path, 'cate_estPS_trun_',interv1,"_",interv2,"_", index, '.dat'))
      }
      
    }
  }
  
  
}





