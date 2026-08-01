# devtools::install_github("lingxiaozhou/geocausal_fork",ref = "adaptiveIntervention")

library(geocausal)

library(splines)
library(spatstat)
library(tidyverse)
library(dplyr)

outpath <- "~/GitHub/CATE-SpatioTemporal-Causal/Applications/Results/"
save_res <- TRUE

airstr <- airstrikes
insurg <- insurgencies
airstr_2006 <- airstrikes_2006 %>%
  filter(type == "Airstrike") %>%   # keep only rows where type == "Airstrike"
  select(-type)                     # then drop the column 'type'

airstr$time <- as.numeric(airstr$date - min(airstr$date) + 1)
insurg$time <- as.numeric(insurg$date - min(insurg$date) + 1)


treatment_hfr <- get_hfr(data = airstr, col = "type", window = iraq_window,
                         time_col = "time", time_range = c(1, 499),
                         coordinates = c("longitude", "latitude"),
                         combine = TRUE)
outcome_hfr <- get_hfr(data = insurg, col = "type", window = iraq_window,
                       time_col = "time", time_range = c(1, 499),
                       coordinates = c("longitude", "latitude"),
                       combine = TRUE)



dat_hfr <- spatstat.geom::cbind.hyperframe(treatment_hfr, outcome_hfr[, -1])
names(dat_hfr)[names(dat_hfr) == "all_combined"] <- "all_treatment"
names(dat_hfr)[names(dat_hfr) == "all_combined.1"] <- "all_outcome"





#-------------------------transform covariates-----------------
# dataset version 4.0
library(httr)

# helper function for loading data
read_dataverse_rds <- function(fileid, server = "https://dataverse.harvard.edu") {
  tf <- tempfile(fileext = ".rds")
  httr::GET(
    paste0(server, "/api/access/datafile/", fileid),
    query = list(format = "original"),
    httr::write_disk(tf, overwrite = TRUE)
  ) |> httr::stop_for_status()
  readRDS(tf)
}


# cities
cities_dist <- read_dataverse_rds(8079278)
cities_dist <- lapply(1:5, function(x) exp(-(2 * x) * cities_dist$distance[[x]]))
dat_hfr$all_cities <- cities_dist[[1]] + cities_dist[[2]] + 
  cities_dist[[3]] + cities_dist[[4]] + cities_dist[[5]]

# rivers
rivers_dist   <- read_dataverse_rds(8079275)
dat_hfr$rivers_dist <- exp(-3 * rivers_dist)

# routes
routes_dist <- read_dataverse_rds(8079274)
dat_hfr$routes_dist <- exp(-3 * routes_dist)

# settles
settle_dist <- read_dataverse_rds(8079279)
for (jj in 1 : 18) {
  dat_hfr$settles <- exp(-12 * settle_dist$distance[[jj]])
  col_name <- paste0("settles_dist_", jj)
  names(dat_hfr)[length(names(dat_hfr))]<- col_name
}


#---------------------get treatment and outcome histories--------

get_out_hx <- function(timeint, out, window, lag) {
  temp <- lapply(timeint, get_hist, Xt = out, lag = lag, window = window)
  purrr::map(temp, 1)
}


dat_hfr$hist_allout_1 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$all_outcome,
                                    window = iraq_window, lag = 1)
dat_hfr$hist_allout_7 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$all_outcome,
                                    window = iraq_window, lag = 7)
dat_hfr$hist_allout_30 <- get_out_hx(timeint = dat_hfr$time, 
                                     out = dat_hfr$all_outcome,
                                     window = iraq_window, lag = 30)



dat_hfr$hist_alltre_1 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$Airstrike,
                                    window = iraq_window, lag = 1)
dat_hfr$hist_alltre_7 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$Airstrike,
                                    window = iraq_window, lag = 7)
dat_hfr$hist_alltre_30 <- get_out_hx(timeint = dat_hfr$time, 
                                     out = dat_hfr$Airstrike,
                                     window = iraq_window, lag = 30)



dat_hfr$hist_sof_1 <- get_out_hx(timeint = dat_hfr$time, 
                                 out = dat_hfr$SOF,
                                 window = iraq_window, lag = 1)
dat_hfr$hist_sof_7 <- get_out_hx(timeint = dat_hfr$time, 
                                 out = dat_hfr$SOF,
                                 window = iraq_window, lag = 7)
dat_hfr$hist_sof_30 <- get_out_hx(timeint = dat_hfr$time, 
                                  out = dat_hfr$SOF,
                                  window = iraq_window, lag = 30)


# help function for transforming the point pattern to distance to the point pattern
get_prior_hist <- function(x, prior_coef, window) {
  
  dist_mp <- spatstat.geom::distmap(x)
  dims <- dim(dist_mp$v)
  if (length(unique(as.numeric(dist_mp$v))) == 1) {
    dist_mp$v <- matrix(Inf, nrow = dims[1], ncol = dims[2])
  }
  prior_hx <- exp(- prior_coef * dist_mp)
  spatstat.geom::Window(prior_hx) <- window
  
  return(list(prior_hx = prior_hx))

}

dat_hfr$hist_allout_1 <- purrr::map(lapply(dat_hfr$hist_allout_1, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_allout_7 <- purrr::map(lapply(dat_hfr$hist_allout_7, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_allout_30 <- purrr::map(lapply(dat_hfr$hist_allout_30, get_prior_hist, 
                                            prior_coef = 6, window = iraq_window), 1)



dat_hfr$hist_alltre_1 <- purrr::map(lapply(dat_hfr$hist_alltre_1, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_alltre_7 <- purrr::map(lapply(dat_hfr$hist_alltre_7, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_alltre_30 <- purrr::map(lapply(dat_hfr$hist_alltre_30, get_prior_hist, 
                                            prior_coef = 6, window = iraq_window), 1)



dat_hfr$hist_sof_1 <- purrr::map(lapply(dat_hfr$hist_sof_1, get_prior_hist, 
                                        prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_sof_7 <- purrr::map(lapply(dat_hfr$hist_sof_7, get_prior_hist, 
                                        prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_sof_30 <- purrr::map(lapply(dat_hfr$hist_sof_30, get_prior_hist, 
                                         prior_coef = 6, window = iraq_window), 1)


#-------------------other covariate--------------------------

iraq_district <- read_dataverse_rds(8138903)
iraq_district <- sf::st_as_sf(iraq_district)

aid_district <- read_dataverse_rds(8079277)





# Create image data for aid
dim <- 256
date_range <- seq(range(airstr$date)[1], range(airstr$date)[2], by = "1 day")
dat_hfr$aid <- NA
get_aid <- function (i) {
  aid_oneday <- aid_district[which(aid_district$USE_DATE == date_range[i]), 
                             c("District", "log_construct")]
  aid_oneday <- dplyr::left_join(iraq_district, aid_oneday, 
                                 by = c("ADM3NAME" = "District"))
  aid_oneday$log_construct <- ifelse(is.na(aid_oneday$log_construct), 0, 
                                     aid_oneday$log_construct)
  
  r <- terra::rast(terra::ext(aid_oneday), nrow = dim, ncol = dim)
  r <- terra::rasterize(terra::vect(aid_oneday), r, field = "log_construct")
  
  raster_matrix <- matrix(as.matrix(r, layer = 1), nrow = dim, ncol = dim)
  raster_matrix <- t(raster_matrix)
  raster_matrix <- raster_matrix[nrow(raster_matrix):1, ]
  
  return(spatstat.geom::as.im(raster_matrix, iraq_window))
}  

dat_hfr$aid <- lapply(dat_hfr$time, get_aid)


# get ethnicity
ethnicity <- read_dataverse_rds(8079276)



# get log population
dat_hfr$pop <- log(ethnicity$total_pop)

#-------------------------non-spatial covariates-------------------------

# splines
time_splines <- splines::ns(1 : nrow(dat_hfr), df = 3)
dat_hfr$time_1 <- time_splines[, 1]
dat_hfr$time_2 <- time_splines[, 2]
dat_hfr$time_3 <- time_splines[, 3]


# surge
min_s <- as.numeric(as.Date("2007-03-25") - min(airstr$date) + 1)
max_s <- as.numeric(as.Date("2008-01-01") - min(airstr$date) + 1)
dat_hfr$surge <- ifelse(dat_hfr$time >= min_s & dat_hfr$time <= max_s, 1, 0)



#----------------------estimate propensity scores------------------------

rhs_var <- c("pop", "rivers_dist", "routes_dist", "all_cities", "aid",
             "hist_alltre_1", "hist_alltre_7", "hist_alltre_30", 
             "hist_allout_1", "hist_allout_7", "hist_allout_30",
             "hist_sof_1", "hist_sof_7", "hist_sof_30",
             "time_1", "time_2", "time_3",
             "settles_dist_1", "settles_dist_2", "settles_dist_3", 
             "settles_dist_4", "settles_dist_5", "settles_dist_6", 
             "settles_dist_7", "settles_dist_8", "settles_dist_9", 
             "settles_dist_10", "settles_dist_11", "settles_dist_12", 
             "settles_dist_13", "settles_dist_14", "settles_dist_15",
             "settles_dist_16", "settles_dist_17", "settles_dist_18",
             "surge")
obs_density <- get_obs_dens(dat_hfr[c(31:nrow(dat_hfr)), ], "Airstrike", 
                            rhs_var, ngrid = 100, window = iraq_window)



#----------------------Performing out-of-sample prediction---------------------
obs_density_8 <- predict_obs_dens(dat_hfr[c(31:nrow(dat_hfr)), ], 
                                  ratio = 0.8, "Airstrike", 
                                  rhs_var, ngrid = 100, iraq_window)


#-----------get baseline distribution for intervention-------------------------
# Use airstrike locations in 2006 to obtain baseline distribution for airstrikes
baseline <- get_base_dens(option = "out", out_data = airstr_2006, 
                          out_coordinates = c("longitude", "latitude"),
                          window = iraq_window, ndim = 100)
plot(baseline, main = "Baseline density", window = iraq_window)

cf_1 <- get_cf_dens(base_dens = baseline, expected_number = 1,
                    power_dens = NA, window = iraq_window) 
cf_6 <- get_cf_dens(base_dens = baseline, expected_number = 6,
                    power_dens = NA, window = iraq_window)

dat_hfr$pixel_count_IED <- pixel_count_ppp(data = dat_hfr$IED)
dat_hfr$pixel_count_SAF <- pixel_count_ppp(data = dat_hfr$SAF)

aid_lag <- 1
aid_lag_name <- "onemonth"

#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.9,0.95,0.98)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        
        
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        em <- em>0
        E.mat <- matrix(em,ncol = 1)
        
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        points <- matrix(c(0,1),ncol=1)
        cat("Start M =",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points)

        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}











#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
    }
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        E.mat <- array(em,c(128,128,length(31:(498-mm+1))))
        
        myindex <- E.mat>0 & !is.na(E.mat)
        myindex2 <- E.mat==0 & !is.na(E.mat)
        
        # choose the knots
        e_min <- min(E.mat[myindex],na.rm = TRUE)
        
        e_max <- max(E.mat[myindex],na.rm = TRUE)
        
        
        my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
        
        E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
        E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
        
        E.basis[myindex] <- E.basis.tmp
        E.basis[myindex2] <- 0
        E.basis <- cbind(E.basis,c(E.mat==0))
        
        tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
        points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
        points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
        points <- c(0,points)
        points.basis <- rbind(0,points.basis)
        points.basis <- cbind(points.basis,c(1,rep(0,20)))
        points.basis[1,1] <- 1
        points.basis <- points.basis[,-1]
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis,bound = 3,M=mm)
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}









#-------------------Binary aid analysis+intensity------------------------------
for (al in 1:length(aid_lag)) {
  cat("Start al: ",al, "\n")
  
  save_path <- paste0(outpath,aid_lag_name[al],"/")
  
  
  load(paste0(outpath,"phi_im.dat"))
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      if(file.exists(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat"))){
        cat("Exist!\n")
      }else{

        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,dimyx = c(128,128))
        em <- em>0
        E.basis <- matrix(em,ncol = 1)

        
        
        
        inten <- rep(c(as.matrix(phi_im6-phi_im1)),length(31:(498-mm+1)))
        inten <- sqrt(inten)
        imedian <- median(unique(c(as.matrix(phi_im6-phi_im1))),na.rm=T)
        imedian <- sqrt(imedian)
        E.basis <- cbind(E.basis,inten,E.basis*inten)
        
        points <- c(1)
        #points.basis <- cbind(1,points)
        points.basis <- cbind(points,imedian,imedian)
        #points.basis <- cbind(points,0.01,0.01)
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=0.95, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis)
        plot(res)
        points.basis <- cbind(0,points.basis)
        
        est_diff <- points.basis%*%res$est_beta
        est_diff_var <- points.basis%*%res$V_beta%*%t(points.basis)
        
        (est_diff+qnorm(c(0.025,0.975))*sqrt(est_diff_var))*81
        
        save(res,est_diff,est_diff_var,file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat"))
      }
      
    }
    
  }
  
  
}


#------------------------Adaptive intervention-----------------------------------

outpath <- "~/GitHub/CATE-SpatioTemporal-Causal/Applications/Results/adaptive/"


airstr2006 <- airstrikes_2006
airstr2006$time <- as.numeric(airstr2006$date - min(airstr2006$date) + 1)
insurg2006 <- insurgencies_2006
insurg2006$time <- as.numeric(insurg2006$date - min(insurg2006$date) + 1)





treatment_hfr_baseline <- get_hfr(data = airstr2006, col = "type", window = iraq_window,
                                  time_col = "time", time_range = c(1, 276),
                                  coordinates = c("longitude", "latitude"),
                                  combine = TRUE)
outcome_hfr_baseline <- get_hfr(data = insurg2006, col = "type", window = iraq_window,
                                time_col = "time", time_range = c(1, 276),
                                coordinates = c("longitude", "latitude"),
                                combine = TRUE)


dat_hfr_baseline <- spatstat.geom::cbind.hyperframe(treatment_hfr_baseline, outcome_hfr_baseline[, -1])
names(dat_hfr_baseline)[names(dat_hfr_baseline) == "all_combined"] <- "all_treatment"
names(dat_hfr_baseline)[names(dat_hfr_baseline) == "all_combined.1"] <- "all_outcome"


#-------------------------transform covariates-----------------

# cities
cities_dist   <- read_dataverse_rds(8079278)
cities_dist <- lapply(1:5, function(x) exp(-(2 * x) * cities_dist$distance[[x]]))
dat_hfr_baseline$all_cities <- cities_dist[[1]] + cities_dist[[2]] + 
  cities_dist[[3]] + cities_dist[[4]] + cities_dist[[5]]

# rivers
rivers_dist   <- read_dataverse_rds(8079275)
dat_hfr_baseline$rivers_dist <- exp(-3 * rivers_dist)

# routes
routes_dist   <- read_dataverse_rds(8079274)
dat_hfr_baseline$routes_dist <- exp(-3 * routes_dist)

# settles
settle_dist   <- read_dataverse_rds(8079279)
for (jj in 1 : 18) {
  dat_hfr_baseline$settles <- exp(-12 * settle_dist$distance[[jj]])
  col_name <- paste0("settles_dist_", jj)
  names(dat_hfr_baseline)[length(names(dat_hfr_baseline))]<- col_name
}


#---------------------get treatment and outcome histories--------

get_out_hx <- function(timeint, out, window, lag) {
  temp <- lapply(timeint, get_hist, Xt = out, lag = lag, window = window)
  purrr::map(temp, 1)
}


dat_hfr_baseline$hist_allout_1 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$all_outcome,
                                             window = iraq_window, lag = 1)
dat_hfr_baseline$hist_allout_7 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$all_outcome,
                                             window = iraq_window, lag = 7)
dat_hfr_baseline$hist_allout_30 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                              out = dat_hfr_baseline$all_outcome,
                                              window = iraq_window, lag = 30)



dat_hfr_baseline$hist_alltre_1 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$Airstrike,
                                             window = iraq_window, lag = 1)
dat_hfr_baseline$hist_alltre_7 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$Airstrike,
                                             window = iraq_window, lag = 7)
dat_hfr_baseline$hist_alltre_30 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                              out = dat_hfr_baseline$Airstrike,
                                              window = iraq_window, lag = 30)



dat_hfr_baseline$hist_sof_1 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                          out = dat_hfr_baseline$SOF,
                                          window = iraq_window, lag = 1)
dat_hfr_baseline$hist_sof_7 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                          out = dat_hfr_baseline$SOF,
                                          window = iraq_window, lag = 7)
dat_hfr_baseline$hist_sof_30 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                           out = dat_hfr_baseline$SOF,
                                           window = iraq_window, lag = 30)


# help function for transforming the point pattern to distance to the point pattern
get_prior_hist <- function(x, prior_coef, window) {
  
  dist_mp <- spatstat.geom::distmap(x)
  dims <- dim(dist_mp$v)
  if (length(unique(as.numeric(dist_mp$v))) == 1) {
    dist_mp$v <- matrix(Inf, nrow = dims[1], ncol = dims[2])
  }
  prior_hx <- exp(- prior_coef * dist_mp)
  spatstat.geom::Window(prior_hx) <- window
  
  return(list(prior_hx = prior_hx))
  
}

dat_hfr_baseline$hist_allout_1 <- purrr::map(lapply(dat_hfr_baseline$hist_allout_1, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_allout_7 <- purrr::map(lapply(dat_hfr_baseline$hist_allout_7, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_allout_30 <- purrr::map(lapply(dat_hfr_baseline$hist_allout_30, get_prior_hist, 
                                                     prior_coef = 6, window = iraq_window), 1)



dat_hfr_baseline$hist_alltre_1 <- purrr::map(lapply(dat_hfr_baseline$hist_alltre_1, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_alltre_7 <- purrr::map(lapply(dat_hfr_baseline$hist_alltre_7, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_alltre_30 <- purrr::map(lapply(dat_hfr_baseline$hist_alltre_30, get_prior_hist, 
                                                     prior_coef = 6, window = iraq_window), 1)



dat_hfr_baseline$hist_sof_1 <- purrr::map(lapply(dat_hfr_baseline$hist_sof_1, get_prior_hist, 
                                                 prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_sof_7 <- purrr::map(lapply(dat_hfr_baseline$hist_sof_7, get_prior_hist, 
                                                 prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_sof_30 <- purrr::map(lapply(dat_hfr_baseline$hist_sof_30, get_prior_hist, 
                                                  prior_coef = 6, window = iraq_window), 1)
# get ethnicity
ethnicity <- read_dataverse_rds(8079276)


# get log population
dat_hfr_baseline$pop <- log(ethnicity$total_pop)



rhs_var_baseline <- c("pop", "rivers_dist", "routes_dist", "all_cities",
                      "hist_allout_1", "hist_allout_7", "hist_allout_30",
                      "hist_sof_1", "hist_sof_7", "hist_sof_30",
                      "settles_dist_1", "settles_dist_2", "settles_dist_3", 
                      "settles_dist_4", "settles_dist_5", "settles_dist_6", 
                      "settles_dist_7", "settles_dist_8", "settles_dist_9", 
                      "settles_dist_10", "settles_dist_11", "settles_dist_12", 
                      "settles_dist_13", "settles_dist_14", "settles_dist_15",
                      "settles_dist_16", "settles_dist_17", "settles_dist_18")


baseline_den <- get_adaptive_baseline_dens(dat_hfr_baseline[c(31:nrow(dat_hfr_baseline)), ],dat_hfr[c(31:nrow(dat_hfr)), ],"Airstrike",
                                           rhs_var_baseline, ngrid = 100, window = iraq_window)




cf_1 <- get_cf_dens_adaptive(baseline_den, scale_factor=1/2)
cf_2 <- get_cf_dens_adaptive(baseline_den, scale_factor=2)

aid_lag <- 1
aid_lag_name <- "onemonth"

#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        
        
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        em <- em>0
        E.mat <- matrix(em,ncol = 1)
        
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        points <- matrix(c(0,1),ncol=1)
        cat("Start M =",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_2, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points)
        
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}





#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
    }
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        E.mat <- array(em,c(128,128,length(31:(498-mm+1))))
        
        myindex <- E.mat>0 & !is.na(E.mat)
        myindex2 <- E.mat==0 & !is.na(E.mat)
        
        # choose the knots
        e_min <- min(E.mat[myindex],na.rm = TRUE)
        
        e_max <- max(E.mat[myindex],na.rm = TRUE)
        
        
        my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
        
        E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
        E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
        
        E.basis[myindex] <- E.basis.tmp
        E.basis[myindex2] <- 0
        E.basis <- cbind(E.basis,c(E.mat==0))
        
        tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
        points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
        points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
        points <- c(0,points)
        points.basis <- rbind(0,points.basis)
        points.basis <- cbind(points.basis,c(1,rep(0,20)))
        points.basis[1,1] <- 1
        points.basis <- points.basis[,-1]
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_2, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis,bound = 3,M=mm)
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}




#-------------------------sensitivity:pixel---------------------
pixel_x_vec <- c(256,64,32)
pixel_y_vec <- c(256,64,32)


aid_lag <- 1
aid_lag_name <- "onemonth"

#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    for (pixel_index in 3) {
      pixel_x <- pixel_x_vec[pixel_index]
      pixel_y <- pixel_y_vec[pixel_index]
      
      dat_hfr$pixel_count_IED <- pixel_count_ppp(data = dat_hfr$IED, ngrid = c(pixel_y,pixel_x))
      dat_hfr$pixel_count_SAF <- pixel_count_ppp(data = dat_hfr$SAF, ngrid = c(pixel_y,pixel_x))
    
      cat("Start al: ",al, "\n")
      
      
      save_path <- paste0(outpath,"pixel/",pixel_y,pixel_x,"/")
      
      
      for (oo in c("SAF","IED")) {
        for (mm in 1:10) {
          
          
          em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(pixel_y,pixel_x))
          em <- em>0
          E.mat <- matrix(em,ncol = 1)
          
          
          if(oo=="IED"){
            pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
          }else{
            pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
          }
          
          points <- matrix(c(0,1),ncol=1)
          cat("Start M =",mm,"Outcome: ",oo,"\n")
          res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                          E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                          nbase = 6, spline_type = "ns",intercept = TRUE,
                          eval_values = points, eval_mat = points)
          
          if(save_res){
            save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
          }
        }
          
        }
    }
    
    
  }
}











#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- FALSE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    for (pixel_index in 1:3) {
      pixel_x <- pixel_x_vec[pixel_index]
      pixel_y <- pixel_y_vec[pixel_index]
      dat_hfr$pixel_count_IED <- pixel_count_ppp(data = dat_hfr$IED, ngrid = c(pixel_y,pixel_x))
      dat_hfr$pixel_count_SAF <- pixel_count_ppp(data = dat_hfr$SAF, ngrid = c(pixel_y,pixel_x))
      
      cat("Start al: ",al, "\n")
      
      if(firsttime){
        get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
          res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
          return(res)
        }
        dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
      }
      
      
      save_path <- paste0(outpath,"pixel/",pixel_y,pixel_x,"/")
      
      
      for (oo in c("SAF","IED")) {
        for (mm in c(3,7,10)) {
          
          
          em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(pixel_y,pixel_x))
          E.mat <- array(em,c(pixel_y,pixel_x,length(31:(498-mm+1))))
          
          myindex <- E.mat>0 & !is.na(E.mat)
          myindex2 <- E.mat==0 & !is.na(E.mat)
          
          # choose the knots
          e_min <- min(E.mat[myindex],na.rm = TRUE)
          
          e_max <- max(E.mat[myindex],na.rm = TRUE)
          
          
          my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
          
          E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
          E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
          
          E.basis[myindex] <- E.basis.tmp
          E.basis[myindex2] <- 0
          E.basis <- cbind(E.basis,c(E.mat==0))
          
          tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
          points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
          points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
          points <- c(0,points)
          points.basis <- rbind(0,points.basis)
          points.basis <- cbind(points.basis,c(1,rep(0,20)))
          points.basis[1,1] <- 1
          points.basis <- points.basis[,-1]
          
          if(oo=="IED"){
            pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
          }else{
            pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
          }
          
          
          cat("Start M=",mm,"Outcome: ",oo,"\n")
          res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                          E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                          nbase = 6, spline_type = "ns",intercept = TRUE,
                          eval_values = points, eval_mat = points.basis,bound = 3,M=mm,dimyx = c(pixel_y,pixel_x))
          if(save_res){
            save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
          }
          
        }
      }
    }
    
  }
}









#--------------------------district---------------------------


## Step 1: add numeric district ID
iraq_district$district_id <- seq_len(nrow(iraq_district))

## Step 2: convert polygon object to terra SpatVector
iraq_district_vect <- terra::vect(iraq_district)

## Step 3: get extent from spatstat window
xr <- iraq_window$xrange
yr <- iraq_window$yrange

## Step 4: create  raster on the same window
r_id <- terra::rast(
  xmin = xr[1], xmax = xr[2],
  ymin = yr[1], ymax = yr[2],
  nrow = 256, ncol = 256
)

## Step 5: rasterize district IDs
r_id <- terra::rasterize(iraq_district_vect, r_id, field = "district_id")

## Step 6: convert raster to matrix with the same orientation as your previous code
id_matrix <- matrix(as.matrix(r_id, layer = 1), nrow = 256, ncol = 256)
id_matrix <- t(id_matrix)
id_matrix <- id_matrix[nrow(id_matrix):1, ]

## Step 7: convert to spatstat image
district_im <- spatstat.geom::as.im(id_matrix, W = iraq_window)


aggregate_to_district <- function(im_obj, district_im, n_district = 104,
                                  district_names = NULL,
                                  fun = c("mean", "sum")) {
  
  ## match argument
  fun <- match.arg(fun)
  
  ## Step 1: force resolution to be the same as district_im
  if (!all(im_obj$dim == district_im$dim)) {
    im_obj <- spatstat.geom::as.im(
      im_obj,
      W = spatstat.geom::Window(district_im),
      dimyx = district_im$dim
    )
  }
  
  ## Step 2: extract matrices
  val_mat <- spatstat.geom::as.matrix.im(im_obj)
  id_mat  <- spatstat.geom::as.matrix.im(district_im)
  
  ## choose aggregation function
  agg_fun <- switch(fun,
                    mean = function(x) mean(x, na.rm = TRUE),
                    sum  = function(x) sum(x, na.rm = TRUE))
  
  ## Step 3: compute district aggregates
  out <- sapply(seq_len(n_district), function(k) {
    v <- val_mat[id_mat == k]
    if (length(v) == 0 || all(is.na(v))) {
      NA_real_
    } else {
      agg_fun(v)
    }
  })
  
  ## Step 4: return 2 x n_district/2 matrix
  out_mat <- matrix(out, nrow = 1)
  
  if (!is.null(district_names)) {
    colnames(out_mat) <- district_names
  }
  
  out_im <- spatstat.geom::as.im(array(out_mat, c(2, n_district / 2)))
  
  return(out_im)
}


dat_hfr$aid_dist <- lapply(dat_hfr$aid, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME
  )
})

dat_hfr$pixel_count_IED_dist <- lapply(dat_hfr$pixel_count_IED, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME,
    fun = "sum"
  )
})

dat_hfr$pixel_count_SAF_dist <- lapply(dat_hfr$pixel_count_SAF, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME,
    fun = "sum"
  )
})

dat_hfr$pop_dist <- lapply(dat_hfr$pop, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME
  )
})



#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    
      
      cat("Start al: ",al, "\n")
      
      
      save_path <- paste0(outpath,"district/")
      
      
      for (oo in c("SAF","IED")) {
        for (mm in 1:10) {
          
          
          em <- get_em_vec(dat_hfr$aid_dist[31:499],time_after = TRUE,lag=mm,entire_window = NULL)
          em <- em>0
          E.mat <- matrix(em,ncol = 1)
          
          
          if(oo=="IED"){
            pixel_count_out <- dat_hfr$pixel_count_IED_dist[31:499]
          }else{
            pixel_count_out <- dat_hfr$pixel_count_SAF_dist[31:499]
          }
          
          points <- matrix(c(0,1),ncol=1)
          cat("Start M =",mm,"Outcome: ",oo,"\n")
          res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                          E_list = dat_hfr$aid_dist[31:499],E_mat = E.mat,
                          nbase = 6, spline_type = "ns",intercept = TRUE,
                          eval_values = points, eval_mat = points)
          
          if(save_res){
            save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
          }
        }
        
      }
    
    
    
  }
}











#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {

      
      cat("Start al: ",al, "\n")
      
      
    save_path <- paste0(outpath,"district/")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap_dist <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid_dist, pop=dat_hfr$pop_dist,dimyx = c(2,52)))
    }
      
      
      for (oo in c("SAF","IED")) {
        for (mm in c(3,7,10)) {
          
          
          em <- get_em_vec(dat_hfr$aid_per_cap_dist[31:499],time_after = TRUE,lag=mm,entire_window = NULL)
          E.mat <- array(em,c(2,52,length(31:(498-mm+1))))
          
          myindex <- E.mat>0 & !is.na(E.mat)
          myindex2 <- E.mat==0 & !is.na(E.mat)
          
          # choose the knots
          e_min <- min(E.mat[myindex],na.rm = TRUE)
          
          e_max <- max(E.mat[myindex],na.rm = TRUE)
          
          
          my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
          
          E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
          E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
          
          E.basis[myindex] <- E.basis.tmp
          E.basis[myindex2] <- 0
          E.basis <- cbind(E.basis,c(E.mat==0))
          
          tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
          points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
          points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
          points <- c(0,points)
          points.basis <- rbind(0,points.basis)
          points.basis <- cbind(points.basis,c(1,rep(0,20)))
          points.basis[1,1] <- 1
          points.basis <- points.basis[,-1]
          
          if(oo=="IED"){
            pixel_count_out <- dat_hfr$pixel_count_IED_dist[31:499]
          }else{
            pixel_count_out <- dat_hfr$pixel_count_SAF_dist[31:499]
          }
          
          
          cat("Start M=",mm,"Outcome: ",oo,"\n")
          res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                          E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                          nbase = 6, spline_type = "ns",intercept = TRUE,
                          eval_values = points, eval_mat = points.basis,bound = 3,M=mm,dimyx = c(2,52))
          if(save_res){
            save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
          }
          
        }
      }
    
    
  }
}


#--------------------------------------sensitivity analysis----------------

truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    
    
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,"sensitivity/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        
        
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid =c(128,128))
        em <- em>0
        E.mat <- matrix(em,ncol = 1)
        
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        points <- matrix(c(0,1),ncol=1)
        cat("Start M =",mm,"Outcome: ",oo,"\n")
        res <- get_cate_sensitivity(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points,gamma_vals = seq(1,1.5,0.01))
        
        print(res$gamma_summary)
        
        print(res$gamma_summary_beta)
        
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        }
      }
      
    }
    
    
    
  }
}


#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
    }
    
    
    save_path <- paste0(outpath,"sensitivity/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        E.mat <- array(em,c(128,128,length(31:(498-mm+1))))
        
        myindex <- E.mat>0 & !is.na(E.mat)
        myindex2 <- E.mat==0 & !is.na(E.mat)
        
        # choose the knots
        e_min <- min(E.mat[myindex],na.rm = TRUE)
        
        e_max <- max(E.mat[myindex],na.rm = TRUE)
        
        
        my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
        
        E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
        E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
        
        E.basis[myindex] <- E.basis.tmp
        E.basis[myindex2] <- 0
        E.basis <- cbind(E.basis,c(E.mat==0))
        
        tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
        points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
        points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
        points <- c(0,points)
        points.basis <- rbind(0,points.basis)
        points.basis <- cbind(points.basis,c(1,rep(0,20)))
        points.basis[1,1] <- 1
        points.basis <- points.basis[,-1]
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate_sensitivity(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis,bound = 3,M=mm,gamma_vals = seq(1,1.5,0.01))
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}



devtools::install_github("lingxiaozhou/geocausal_fork",ref = "adaptiveIntervention")



library(geocausal) 
library(splines)
library(spatstat)
library(tidyverse)
library(dplyr)

outpath <- "~/GitHub/CATE-SpatioTemporal-Causal/Applications/Results/"
save_res <- TRUE

airstr <- airstrikes
insurg <- insurgencies
airstr_2006 <- airstrikes_2006 %>%
  filter(type == "Airstrike") %>%   # keep only rows where type == "Airstrike"
  select(-type)                     # then drop the column 'type'

airstr$time <- as.numeric(airstr$date - min(airstr$date) + 1)
insurg$time <- as.numeric(insurg$date - min(insurg$date) + 1)


treatment_hfr <- get_hfr(data = airstr, col = "type", window = iraq_window,
                         time_col = "time", time_range = c(1, 499),
                         coordinates = c("longitude", "latitude"),
                         combine = TRUE)
outcome_hfr <- get_hfr(data = insurg, col = "type", window = iraq_window,
                       time_col = "time", time_range = c(1, 499),
                       coordinates = c("longitude", "latitude"),
                       combine = TRUE)



dat_hfr <- spatstat.geom::cbind.hyperframe(treatment_hfr, outcome_hfr[, -1])
names(dat_hfr)[names(dat_hfr) == "all_combined"] <- "all_treatment"
names(dat_hfr)[names(dat_hfr) == "all_combined.1"] <- "all_outcome"





#-------------------------transform covariates-----------------
# dataset version 4.0
library(httr)

# helper function for loading data
read_dataverse_rds <- function(fileid, server = "https://dataverse.harvard.edu") {
  tf <- tempfile(fileext = ".rds")
  httr::GET(
    paste0(server, "/api/access/datafile/", fileid),
    query = list(format = "original"),
    httr::write_disk(tf, overwrite = TRUE)
  ) |> httr::stop_for_status()
  readRDS(tf)
}


# cities
cities_dist <- read_dataverse_rds(8079278)
cities_dist <- lapply(1:5, function(x) exp(-(2 * x) * cities_dist$distance[[x]]))
dat_hfr$all_cities <- cities_dist[[1]] + cities_dist[[2]] + 
  cities_dist[[3]] + cities_dist[[4]] + cities_dist[[5]]

# rivers
rivers_dist   <- read_dataverse_rds(8079275)
dat_hfr$rivers_dist <- exp(-3 * rivers_dist)

# routes
routes_dist <- read_dataverse_rds(8079274)
dat_hfr$routes_dist <- exp(-3 * routes_dist)

# settles
settle_dist <- read_dataverse_rds(8079279)
for (jj in 1 : 18) {
  dat_hfr$settles <- exp(-12 * settle_dist$distance[[jj]])
  col_name <- paste0("settles_dist_", jj)
  names(dat_hfr)[length(names(dat_hfr))]<- col_name
}


#---------------------get treatment and outcome histories--------

get_out_hx <- function(timeint, out, window, lag) {
  temp <- lapply(timeint, get_hist, Xt = out, lag = lag, window = window)
  purrr::map(temp, 1)
}


dat_hfr$hist_allout_1 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$all_outcome,
                                    window = iraq_window, lag = 1)
dat_hfr$hist_allout_7 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$all_outcome,
                                    window = iraq_window, lag = 7)
dat_hfr$hist_allout_30 <- get_out_hx(timeint = dat_hfr$time, 
                                     out = dat_hfr$all_outcome,
                                     window = iraq_window, lag = 30)



dat_hfr$hist_alltre_1 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$Airstrike,
                                    window = iraq_window, lag = 1)
dat_hfr$hist_alltre_7 <- get_out_hx(timeint = dat_hfr$time, 
                                    out = dat_hfr$Airstrike,
                                    window = iraq_window, lag = 7)
dat_hfr$hist_alltre_30 <- get_out_hx(timeint = dat_hfr$time, 
                                     out = dat_hfr$Airstrike,
                                     window = iraq_window, lag = 30)



dat_hfr$hist_sof_1 <- get_out_hx(timeint = dat_hfr$time, 
                                 out = dat_hfr$SOF,
                                 window = iraq_window, lag = 1)
dat_hfr$hist_sof_7 <- get_out_hx(timeint = dat_hfr$time, 
                                 out = dat_hfr$SOF,
                                 window = iraq_window, lag = 7)
dat_hfr$hist_sof_30 <- get_out_hx(timeint = dat_hfr$time, 
                                  out = dat_hfr$SOF,
                                  window = iraq_window, lag = 30)


# help function for transforming the point pattern to distance to the point pattern
get_prior_hist <- function(x, prior_coef, window) {
  
  dist_mp <- spatstat.geom::distmap(x)
  dims <- dim(dist_mp$v)
  if (length(unique(as.numeric(dist_mp$v))) == 1) {
    dist_mp$v <- matrix(Inf, nrow = dims[1], ncol = dims[2])
  }
  prior_hx <- exp(- prior_coef * dist_mp)
  spatstat.geom::Window(prior_hx) <- window
  
  return(list(prior_hx = prior_hx))
  
}

dat_hfr$hist_allout_1 <- purrr::map(lapply(dat_hfr$hist_allout_1, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_allout_7 <- purrr::map(lapply(dat_hfr$hist_allout_7, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_allout_30 <- purrr::map(lapply(dat_hfr$hist_allout_30, get_prior_hist, 
                                            prior_coef = 6, window = iraq_window), 1)



dat_hfr$hist_alltre_1 <- purrr::map(lapply(dat_hfr$hist_alltre_1, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_alltre_7 <- purrr::map(lapply(dat_hfr$hist_alltre_7, get_prior_hist, 
                                           prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_alltre_30 <- purrr::map(lapply(dat_hfr$hist_alltre_30, get_prior_hist, 
                                            prior_coef = 6, window = iraq_window), 1)



dat_hfr$hist_sof_1 <- purrr::map(lapply(dat_hfr$hist_sof_1, get_prior_hist, 
                                        prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_sof_7 <- purrr::map(lapply(dat_hfr$hist_sof_7, get_prior_hist, 
                                        prior_coef = 6, window = iraq_window), 1)
dat_hfr$hist_sof_30 <- purrr::map(lapply(dat_hfr$hist_sof_30, get_prior_hist, 
                                         prior_coef = 6, window = iraq_window), 1)


#-------------------other covariate--------------------------

iraq_district <- read_dataverse_rds(8138903)
iraq_district <- sf::st_as_sf(iraq_district)

aid_district <- read_dataverse_rds(8079277)





# Create image data for aid
dim <- 256
date_range <- seq(range(airstr$date)[1], range(airstr$date)[2], by = "1 day")
dat_hfr$aid <- NA
get_aid <- function (i) {
  aid_oneday <- aid_district[which(aid_district$USE_DATE == date_range[i]), 
                             c("District", "log_construct")]
  aid_oneday <- dplyr::left_join(iraq_district, aid_oneday, 
                                 by = c("ADM3NAME" = "District"))
  aid_oneday$log_construct <- ifelse(is.na(aid_oneday$log_construct), 0, 
                                     aid_oneday$log_construct)
  
  r <- terra::rast(terra::ext(aid_oneday), nrow = dim, ncol = dim)
  r <- terra::rasterize(terra::vect(aid_oneday), r, field = "log_construct")
  
  raster_matrix <- matrix(as.matrix(r, layer = 1), nrow = dim, ncol = dim)
  raster_matrix <- t(raster_matrix)
  raster_matrix <- raster_matrix[nrow(raster_matrix):1, ]
  
  return(spatstat.geom::as.im(raster_matrix, iraq_window))
}  

dat_hfr$aid <- lapply(dat_hfr$time, get_aid)


# get ethnicity
ethnicity <- read_dataverse_rds(8079276)



# get log population
dat_hfr$pop <- log(ethnicity$total_pop)

#-------------------------non-spatial covariates-------------------------

# splines
time_splines <- splines::ns(1 : nrow(dat_hfr), df = 3)
dat_hfr$time_1 <- time_splines[, 1]
dat_hfr$time_2 <- time_splines[, 2]
dat_hfr$time_3 <- time_splines[, 3]


# surge
min_s <- as.numeric(as.Date("2007-03-25") - min(airstr$date) + 1)
max_s <- as.numeric(as.Date("2008-01-01") - min(airstr$date) + 1)
dat_hfr$surge <- ifelse(dat_hfr$time >= min_s & dat_hfr$time <= max_s, 1, 0)



#----------------------estimate propensity scores------------------------

rhs_var <- c("pop", "rivers_dist", "routes_dist", "all_cities", "aid",
             "hist_alltre_1", "hist_alltre_7", "hist_alltre_30", 
             "hist_allout_1", "hist_allout_7", "hist_allout_30",
             "hist_sof_1", "hist_sof_7", "hist_sof_30",
             "time_1", "time_2", "time_3",
             "settles_dist_1", "settles_dist_2", "settles_dist_3", 
             "settles_dist_4", "settles_dist_5", "settles_dist_6", 
             "settles_dist_7", "settles_dist_8", "settles_dist_9", 
             "settles_dist_10", "settles_dist_11", "settles_dist_12", 
             "settles_dist_13", "settles_dist_14", "settles_dist_15",
             "settles_dist_16", "settles_dist_17", "settles_dist_18",
             "surge")
obs_density <- get_obs_dens(dat_hfr[c(31:nrow(dat_hfr)), ], "Airstrike", 
                            rhs_var, ngrid = 100, window = iraq_window)



#----------------------Performing out-of-sample prediction---------------------
obs_density_8 <- predict_obs_dens(dat_hfr[c(31:nrow(dat_hfr)), ], 
                                  ratio = 0.8, "Airstrike", 
                                  rhs_var, ngrid = 100, iraq_window)


#-----------get baseline distribution for intervention-------------------------
# Use airstrike locations in 2006 to obtain baseline distribution for airstrikes
baseline <- get_base_dens(option = "out", out_data = airstr_2006, 
                          out_coordinates = c("longitude", "latitude"),
                          window = iraq_window, ndim = 100)
plot(baseline, main = "Baseline density", window = iraq_window)

cf_1 <- get_cf_dens(base_dens = baseline, expected_number = 1,
                    power_dens = NA, window = iraq_window) 
cf_6 <- get_cf_dens(base_dens = baseline, expected_number = 6,
                    power_dens = NA, window = iraq_window)

dat_hfr$pixel_count_IED <- pixel_count_ppp(data = dat_hfr$IED)
dat_hfr$pixel_count_SAF <- pixel_count_ppp(data = dat_hfr$SAF)

aid_lag <- 1
aid_lag_name <- "onemonth"

#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.9,0.95,0.98)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        
        
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        em <- em>0
        E.mat <- matrix(em,ncol = 1)
        
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        points <- matrix(c(0,1),ncol=1)
        cat("Start M =",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points)
        
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}











#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
    }
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        E.mat <- array(em,c(128,128,length(31:(498-mm+1))))
        
        myindex <- E.mat>0 & !is.na(E.mat)
        myindex2 <- E.mat==0 & !is.na(E.mat)
        
        # choose the knots
        e_min <- min(E.mat[myindex],na.rm = TRUE)
        
        e_max <- max(E.mat[myindex],na.rm = TRUE)
        
        
        my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
        
        E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
        E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
        
        E.basis[myindex] <- E.basis.tmp
        E.basis[myindex2] <- 0
        E.basis <- cbind(E.basis,c(E.mat==0))
        
        tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
        points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
        points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
        points <- c(0,points)
        points.basis <- rbind(0,points.basis)
        points.basis <- cbind(points.basis,c(1,rep(0,20)))
        points.basis[1,1] <- 1
        points.basis <- points.basis[,-1]
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis,bound = 3,M=mm)
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}









#-------------------Binary aid analysis+intensity------------------------------
for (al in 1:length(aid_lag)) {
  cat("Start al: ",al, "\n")
  
  save_path <- paste0(outpath,aid_lag_name[al],"/")
  
  
  load(paste0(outpath,"phi_im.dat"))
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      if(file.exists(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat"))){
        cat("Exist!\n")
      }else{
        
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,dimyx = c(128,128))
        em <- em>0
        E.basis <- matrix(em,ncol = 1)
        
        
        
        
        inten <- rep(c(as.matrix(phi_im6-phi_im1)),length(31:(498-mm+1)))
        inten <- sqrt(inten)
        imedian <- median(unique(c(as.matrix(phi_im6-phi_im1))),na.rm=T)
        imedian <- sqrt(imedian)
        E.basis <- cbind(E.basis,inten,E.basis*inten)
        
        points <- c(1)
        #points.basis <- cbind(1,points)
        points.basis <- cbind(points,imedian,imedian)
        #points.basis <- cbind(points,0.01,0.01)
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=0.95, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis)
        plot(res)
        points.basis <- cbind(0,points.basis)
        
        est_diff <- points.basis%*%res$est_beta
        est_diff_var <- points.basis%*%res$V_beta%*%t(points.basis)
        
        (est_diff+qnorm(c(0.025,0.975))*sqrt(est_diff_var))*81
        
        save(res,est_diff,est_diff_var,file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat"))
      }
      
    }
    
  }
  
  
}


#------------------------Adaptive intervention-----------------------------------

outpath <- "~/GitHub/CATE-SpatioTemporal-Causal/Applications/Results/adaptive/"


airstr2006 <- airstrikes_2006
airstr2006$time <- as.numeric(airstr2006$date - min(airstr2006$date) + 1)
insurg2006 <- insurgencies_2006
insurg2006$time <- as.numeric(insurg2006$date - min(insurg2006$date) + 1)





treatment_hfr_baseline <- get_hfr(data = airstr2006, col = "type", window = iraq_window,
                                  time_col = "time", time_range = c(1, 276),
                                  coordinates = c("longitude", "latitude"),
                                  combine = TRUE)
outcome_hfr_baseline <- get_hfr(data = insurg2006, col = "type", window = iraq_window,
                                time_col = "time", time_range = c(1, 276),
                                coordinates = c("longitude", "latitude"),
                                combine = TRUE)


dat_hfr_baseline <- spatstat.geom::cbind.hyperframe(treatment_hfr_baseline, outcome_hfr_baseline[, -1])
names(dat_hfr_baseline)[names(dat_hfr_baseline) == "all_combined"] <- "all_treatment"
names(dat_hfr_baseline)[names(dat_hfr_baseline) == "all_combined.1"] <- "all_outcome"


#-------------------------transform covariates-----------------

# cities
cities_dist   <- read_dataverse_rds(8079278)
cities_dist <- lapply(1:5, function(x) exp(-(2 * x) * cities_dist$distance[[x]]))
dat_hfr_baseline$all_cities <- cities_dist[[1]] + cities_dist[[2]] + 
  cities_dist[[3]] + cities_dist[[4]] + cities_dist[[5]]

# rivers
rivers_dist   <- read_dataverse_rds(8079275)
dat_hfr_baseline$rivers_dist <- exp(-3 * rivers_dist)

# routes
routes_dist   <- read_dataverse_rds(8079274)
dat_hfr_baseline$routes_dist <- exp(-3 * routes_dist)

# settles
settle_dist   <- read_dataverse_rds(8079279)
for (jj in 1 : 18) {
  dat_hfr_baseline$settles <- exp(-12 * settle_dist$distance[[jj]])
  col_name <- paste0("settles_dist_", jj)
  names(dat_hfr_baseline)[length(names(dat_hfr_baseline))]<- col_name
}


#---------------------get treatment and outcome histories--------

get_out_hx <- function(timeint, out, window, lag) {
  temp <- lapply(timeint, get_hist, Xt = out, lag = lag, window = window)
  purrr::map(temp, 1)
}


dat_hfr_baseline$hist_allout_1 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$all_outcome,
                                             window = iraq_window, lag = 1)
dat_hfr_baseline$hist_allout_7 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$all_outcome,
                                             window = iraq_window, lag = 7)
dat_hfr_baseline$hist_allout_30 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                              out = dat_hfr_baseline$all_outcome,
                                              window = iraq_window, lag = 30)



dat_hfr_baseline$hist_alltre_1 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$Airstrike,
                                             window = iraq_window, lag = 1)
dat_hfr_baseline$hist_alltre_7 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                             out = dat_hfr_baseline$Airstrike,
                                             window = iraq_window, lag = 7)
dat_hfr_baseline$hist_alltre_30 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                              out = dat_hfr_baseline$Airstrike,
                                              window = iraq_window, lag = 30)



dat_hfr_baseline$hist_sof_1 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                          out = dat_hfr_baseline$SOF,
                                          window = iraq_window, lag = 1)
dat_hfr_baseline$hist_sof_7 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                          out = dat_hfr_baseline$SOF,
                                          window = iraq_window, lag = 7)
dat_hfr_baseline$hist_sof_30 <- get_out_hx(timeint = dat_hfr_baseline$time, 
                                           out = dat_hfr_baseline$SOF,
                                           window = iraq_window, lag = 30)


# help function for transforming the point pattern to distance to the point pattern
get_prior_hist <- function(x, prior_coef, window) {
  
  dist_mp <- spatstat.geom::distmap(x)
  dims <- dim(dist_mp$v)
  if (length(unique(as.numeric(dist_mp$v))) == 1) {
    dist_mp$v <- matrix(Inf, nrow = dims[1], ncol = dims[2])
  }
  prior_hx <- exp(- prior_coef * dist_mp)
  spatstat.geom::Window(prior_hx) <- window
  
  return(list(prior_hx = prior_hx))
  
}

dat_hfr_baseline$hist_allout_1 <- purrr::map(lapply(dat_hfr_baseline$hist_allout_1, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_allout_7 <- purrr::map(lapply(dat_hfr_baseline$hist_allout_7, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_allout_30 <- purrr::map(lapply(dat_hfr_baseline$hist_allout_30, get_prior_hist, 
                                                     prior_coef = 6, window = iraq_window), 1)



dat_hfr_baseline$hist_alltre_1 <- purrr::map(lapply(dat_hfr_baseline$hist_alltre_1, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_alltre_7 <- purrr::map(lapply(dat_hfr_baseline$hist_alltre_7, get_prior_hist, 
                                                    prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_alltre_30 <- purrr::map(lapply(dat_hfr_baseline$hist_alltre_30, get_prior_hist, 
                                                     prior_coef = 6, window = iraq_window), 1)



dat_hfr_baseline$hist_sof_1 <- purrr::map(lapply(dat_hfr_baseline$hist_sof_1, get_prior_hist, 
                                                 prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_sof_7 <- purrr::map(lapply(dat_hfr_baseline$hist_sof_7, get_prior_hist, 
                                                 prior_coef = 6, window = iraq_window), 1)
dat_hfr_baseline$hist_sof_30 <- purrr::map(lapply(dat_hfr_baseline$hist_sof_30, get_prior_hist, 
                                                  prior_coef = 6, window = iraq_window), 1)
# get ethnicity
ethnicity <- read_dataverse_rds(8079276)


# get log population
dat_hfr_baseline$pop <- log(ethnicity$total_pop)



rhs_var_baseline <- c("pop", "rivers_dist", "routes_dist", "all_cities",
                      "hist_allout_1", "hist_allout_7", "hist_allout_30",
                      "hist_sof_1", "hist_sof_7", "hist_sof_30",
                      "settles_dist_1", "settles_dist_2", "settles_dist_3", 
                      "settles_dist_4", "settles_dist_5", "settles_dist_6", 
                      "settles_dist_7", "settles_dist_8", "settles_dist_9", 
                      "settles_dist_10", "settles_dist_11", "settles_dist_12", 
                      "settles_dist_13", "settles_dist_14", "settles_dist_15",
                      "settles_dist_16", "settles_dist_17", "settles_dist_18")


baseline_den <- get_adaptive_baseline_dens(dat_hfr_baseline[c(31:nrow(dat_hfr_baseline)), ],dat_hfr[c(31:nrow(dat_hfr)), ],"Airstrike",
                                           rhs_var_baseline, ngrid = 100, window = iraq_window)




cf_1 <- get_cf_dens_adaptive(baseline_den, scale_factor=1/2)
cf_2 <- get_cf_dens_adaptive(baseline_den, scale_factor=2)

aid_lag <- 1
aid_lag_name <- "onemonth"

#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        
        
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        em <- em>0
        E.mat <- matrix(em,ncol = 1)
        
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        points <- matrix(c(0,1),ncol=1)
        cat("Start M =",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_2, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points)
        
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}





#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    cat("Start al: ",al, "\n")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
    }
    
    
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(128,128))
        E.mat <- array(em,c(128,128,length(31:(498-mm+1))))
        
        myindex <- E.mat>0 & !is.na(E.mat)
        myindex2 <- E.mat==0 & !is.na(E.mat)
        
        # choose the knots
        e_min <- min(E.mat[myindex],na.rm = TRUE)
        
        e_max <- max(E.mat[myindex],na.rm = TRUE)
        
        
        my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
        
        E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
        E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
        
        E.basis[myindex] <- E.basis.tmp
        E.basis[myindex2] <- 0
        E.basis <- cbind(E.basis,c(E.mat==0))
        
        tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
        points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
        points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
        points <- c(0,points)
        points.basis <- rbind(0,points.basis)
        points.basis <- cbind(points.basis,c(1,rep(0,20)))
        points.basis[1,1] <- 1
        points.basis <- points.basis[,-1]
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_2, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis,bound = 3,M=mm)
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}




#-------------------------sensitivity:pixel---------------------
pixel_x_vec <- c(256,64,32)
pixel_y_vec <- c(256,64,32)


aid_lag <- 1
aid_lag_name <- "onemonth"

#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    for (pixel_index in 3) {
      pixel_x <- pixel_x_vec[pixel_index]
      pixel_y <- pixel_y_vec[pixel_index]
      
      dat_hfr$pixel_count_IED <- pixel_count_ppp(data = dat_hfr$IED, ngrid = c(pixel_y,pixel_x))
      dat_hfr$pixel_count_SAF <- pixel_count_ppp(data = dat_hfr$SAF, ngrid = c(pixel_y,pixel_x))
      
      cat("Start al: ",al, "\n")
      
      
      save_path <- paste0(outpath,"pixel/",pixel_y,pixel_x,"/")
      
      
      for (oo in c("SAF","IED")) {
        for (mm in 1:10) {
          
          
          em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(pixel_y,pixel_x))
          em <- em>0
          E.mat <- matrix(em,ncol = 1)
          
          
          if(oo=="IED"){
            pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
          }else{
            pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
          }
          
          points <- matrix(c(0,1),ncol=1)
          cat("Start M =",mm,"Outcome: ",oo,"\n")
          res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                          E_list = dat_hfr$aid[31:499],E_mat = E.mat,
                          nbase = 6, spline_type = "ns",intercept = TRUE,
                          eval_values = points, eval_mat = points)
          
          if(save_res){
            save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
          }
        }
        
      }
    }
    
    
  }
}











#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- FALSE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    for (pixel_index in 1:3) {
      pixel_x <- pixel_x_vec[pixel_index]
      pixel_y <- pixel_y_vec[pixel_index]
      dat_hfr$pixel_count_IED <- pixel_count_ppp(data = dat_hfr$IED, ngrid = c(pixel_y,pixel_x))
      dat_hfr$pixel_count_SAF <- pixel_count_ppp(data = dat_hfr$SAF, ngrid = c(pixel_y,pixel_x))
      
      cat("Start al: ",al, "\n")
      
      if(firsttime){
        get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
          res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
          return(res)
        }
        dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
      }
      
      
      save_path <- paste0(outpath,"pixel/",pixel_y,pixel_x,"/")
      
      
      for (oo in c("SAF","IED")) {
        for (mm in c(3,7,10)) {
          
          
          em <- get_em_vec(dat_hfr$aid_per_cap[31:499],time_after = TRUE,lag=mm,entire_window = NULL,ngrid = c(pixel_y,pixel_x))
          E.mat <- array(em,c(pixel_y,pixel_x,length(31:(498-mm+1))))
          
          myindex <- E.mat>0 & !is.na(E.mat)
          myindex2 <- E.mat==0 & !is.na(E.mat)
          
          # choose the knots
          e_min <- min(E.mat[myindex],na.rm = TRUE)
          
          e_max <- max(E.mat[myindex],na.rm = TRUE)
          
          
          my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
          
          E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
          E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
          
          E.basis[myindex] <- E.basis.tmp
          E.basis[myindex2] <- 0
          E.basis <- cbind(E.basis,c(E.mat==0))
          
          tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
          points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
          points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
          points <- c(0,points)
          points.basis <- rbind(0,points.basis)
          points.basis <- cbind(points.basis,c(1,rep(0,20)))
          points.basis[1,1] <- 1
          points.basis <- points.basis[,-1]
          
          if(oo=="IED"){
            pixel_count_out <- dat_hfr$pixel_count_IED[31:499]
          }else{
            pixel_count_out <- dat_hfr$pixel_count_SAF[31:499]
          }
          
          
          cat("Start M=",mm,"Outcome: ",oo,"\n")
          res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                          E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                          nbase = 6, spline_type = "ns",intercept = TRUE,
                          eval_values = points, eval_mat = points.basis,bound = 3,M=mm,dimyx = c(pixel_y,pixel_x))
          if(save_res){
            save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
          }
          
        }
      }
    }
    
  }
}









#--------------------------district---------------------------


## Step 1: add numeric district ID
iraq_district$district_id <- seq_len(nrow(iraq_district))

## Step 2: convert polygon object to terra SpatVector
iraq_district_vect <- terra::vect(iraq_district)

## Step 3: get extent from spatstat window
xr <- iraq_window$xrange
yr <- iraq_window$yrange

## Step 4: create  raster on the same window
r_id <- terra::rast(
  xmin = xr[1], xmax = xr[2],
  ymin = yr[1], ymax = yr[2],
  nrow = 256, ncol = 256
)

## Step 5: rasterize district IDs
r_id <- terra::rasterize(iraq_district_vect, r_id, field = "district_id")

## Step 6: convert raster to matrix with the same orientation as your previous code
id_matrix <- matrix(as.matrix(r_id, layer = 1), nrow = 256, ncol = 256)
id_matrix <- t(id_matrix)
id_matrix <- id_matrix[nrow(id_matrix):1, ]

## Step 7: convert to spatstat image
district_im <- spatstat.geom::as.im(id_matrix, W = iraq_window)


aggregate_to_district <- function(im_obj, district_im, n_district = 104,
                                  district_names = NULL,
                                  fun = c("mean", "sum")) {
  
  ## match argument
  fun <- match.arg(fun)
  
  ## Step 1: force resolution to be the same as district_im
  if (!all(im_obj$dim == district_im$dim)) {
    im_obj <- spatstat.geom::as.im(
      im_obj,
      W = spatstat.geom::Window(district_im),
      dimyx = district_im$dim
    )
  }
  
  ## Step 2: extract matrices
  val_mat <- spatstat.geom::as.matrix.im(im_obj)
  id_mat  <- spatstat.geom::as.matrix.im(district_im)
  
  ## choose aggregation function
  agg_fun <- switch(fun,
                    mean = function(x) mean(x, na.rm = TRUE),
                    sum  = function(x) sum(x, na.rm = TRUE))
  
  ## Step 3: compute district aggregates
  out <- sapply(seq_len(n_district), function(k) {
    v <- val_mat[id_mat == k]
    if (length(v) == 0 || all(is.na(v))) {
      NA_real_
    } else {
      agg_fun(v)
    }
  })
  
  ## Step 4: return 2 x n_district/2 matrix
  out_mat <- matrix(out, nrow = 1)
  
  if (!is.null(district_names)) {
    colnames(out_mat) <- district_names
  }
  
  out_im <- spatstat.geom::as.im(array(out_mat, c(2, n_district / 2)))
  
  return(out_im)
}


dat_hfr$aid_dist <- lapply(dat_hfr$aid, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME
  )
})

dat_hfr$pixel_count_IED_dist <- lapply(dat_hfr$pixel_count_IED, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME,
    fun = "sum"
  )
})

dat_hfr$pixel_count_SAF_dist <- lapply(dat_hfr$pixel_count_SAF, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME,
    fun = "sum"
  )
})

dat_hfr$pop_dist <- lapply(dat_hfr$pop, function(z) {
  aggregate_to_district(
    im_obj = z,
    district_im = district_im,
    n_district = 104,
    district_names = iraq_district$ADM3NAME
  )
})



#-------------------Binary aid analysis------------------------------
truncation_vec <- c(0.95)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    
    
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,"district/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        
        
        em <- get_em_vec(dat_hfr$aid_dist[31:499],time_after = TRUE,lag=mm,entire_window = NULL)
        em <- em>0
        E.mat <- matrix(em,ncol = 1)
        
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED_dist[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF_dist[31:499]
        }
        
        points <- matrix(c(0,1),ncol=1)
        cat("Start M =",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid_dist[31:499],E_mat = E.mat,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points)
        
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        }
      }
      
    }
    
    
    
  }
}











#-------------------------Continuous aid--------------------------
nbaseE <-5
firsttime <- TRUE

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    
    
    cat("Start al: ",al, "\n")
    
    
    save_path <- paste0(outpath,"district/")
    
    if(firsttime){
      get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
        res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
        return(res)
      }
      dat_hfr$aid_per_cap_dist <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid_dist, pop=dat_hfr$pop_dist,dimyx = c(2,52)))
    }
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        em <- get_em_vec(dat_hfr$aid_per_cap_dist[31:499],time_after = TRUE,lag=mm,entire_window = NULL)
        E.mat <- array(em,c(2,52,length(31:(498-mm+1))))
        
        myindex <- E.mat>0 & !is.na(E.mat)
        myindex2 <- E.mat==0 & !is.na(E.mat)
        
        # choose the knots
        e_min <- min(E.mat[myindex],na.rm = TRUE)
        
        e_max <- max(E.mat[myindex],na.rm = TRUE)
        
        
        my_knots <- e_min+1:(nbaseE-2)*(e_max-e_min)/(nbaseE-1)
        
        E.basis <- array(NA, c(prod(dim(E.mat)),nbaseE-1))
        E.basis.tmp <- ns(c(E.mat[myindex]),knots = my_knots)
        
        E.basis[myindex] <- E.basis.tmp
        E.basis[myindex2] <- 0
        E.basis <- cbind(E.basis,c(E.mat==0))
        
        tmp_points <- quantile(E.mat[myindex],c(0.1,0.9),na.rm = TRUE)
        points <- seq(tmp_points[1],tmp_points[2],length.out = 20)
        points.basis <- cbind(1,predict(E.basis.tmp,newx = points))
        points <- c(0,points)
        points.basis <- rbind(0,points.basis)
        points.basis <- cbind(points.basis,c(1,rep(0,20)))
        points.basis[1,1] <- 1
        points.basis <- points.basis[,-1]
        
        if(oo=="IED"){
          pixel_count_out <- dat_hfr$pixel_count_IED_dist[31:499]
        }else{
          pixel_count_out <- dat_hfr$pixel_count_SAF_dist[31:499]
        }
        
        
        cat("Start M=",mm,"Outcome: ",oo,"\n")
        res <- get_cate(obs=obs_density, cf1 = cf_1, cf2 = cf_6, treat = dat_hfr$Airstrike[31:499], pixel_count_out=pixel_count_out,lag=mm, trunc_level=truncation, time_after=TRUE,entire_window = iraq_window,
                        E_list = dat_hfr$aid[31:499],E_mat = E.basis,
                        nbase = 6, spline_type = "ns",intercept = TRUE,
                        eval_values = points, eval_mat = points.basis,bound = 3,M=mm,dimyx = c(2,52))
        if(save_res){
          save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        }
        
      }
    }
    
    
  }
}







#--------------------------------------sensitivity analysis----------------

truncation_vec <- c(0.95)

gamma_cate <- array(NA,c(2,10,2))
dimnames(gamma_cate) <- list(c("SAF","IED"),1:10,1:2)

gamma_beta <- array(NA,c(2,10,1))
dimnames(gamma_beta) <- list(c("SAF","IED"),1:10,1)


for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {

    save_path <- paste0(outpath,"sensitivity/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
  
        load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        
        gamma_cate[oo,paste0(mm),] <- res$gamma_summary$pointwise_first_gamma
        gamma_beta[oo,paste0(mm),] <- res$gamma_summary_beta$pointwise_first_gamma
          
      }
      
    }
    
  }
}


gamma_cate <- array(NA,c(2,3,21))
dimnames(gamma_cate) <- list(c("SAF","IED"),c(3,7,10),1:21)

gamma_beta <- array(NA,c(2,3,5))
dimnames(gamma_beta) <- list(c("SAF","IED"),c(3,7,10),1:5)

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    save_path <- paste0(outpath,"sensitivity/")
    
    
    for (oo in c("SAF","IED")) {
      for (mm in c(3,7,10)) {
        
        
        load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
        gamma_cate[oo,paste0(mm),] <- res$gamma_summary$pointwise_first_gamma
        gamma_beta[oo,paste0(mm),] <- res$gamma_summary_beta$pointwise_first_gamma
        
      }
    }
    
    
  }
}

apply(gamma_beta, c(1,2), max, na.rm = TRUE)
