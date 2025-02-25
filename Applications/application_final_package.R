obs_density <- readRDS("~/Iraq_causal/application_code/obs_density.rds")
dat_hfr <- readRDS("~/Iraq_causal/application_code/dat_hfr.rds")
library(geocausal)
library(splines)
library(spatstat)
library(ggplot2)

airstr <- airstrikes
insurg <- insurgencies
airstr_2006 <- airstrikes_base


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




# dat_hfr$smooth_allout_mcl <- smooth_ppp(data = dat_hfr$all_outcome,
#                                         method = "mclust", sampling = 0.05)
# 
# 
# 
# dat_hfr$smooth_allout <- smooth_ppp(data = dat_hfr$all_outcome, 
#                                     method = "abramson")
# dat_hfr$smooth_IED <- smooth_ppp(data = dat_hfr$IED, 
#                                  method = "abramson")
# dat_hfr$smooth_SAF <- smooth_ppp(data = dat_hfr$SAF, 
#                                  method = "abramson")


#-------------------------transform covariates-----------------

# cities
cities_dist <- dataverse::get_dataframe_by_name(
  filename = "cities_dist.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)
cities_dist <- lapply(1:5, function(x) exp(-(2 * x) * cities_dist$distance[[x]]))
dat_hfr$all_cities <- cities_dist[[1]] + cities_dist[[2]] + 
  cities_dist[[3]] + cities_dist[[4]] + cities_dist[[5]]

# rivers
rivers_dist <- dataverse::get_dataframe_by_name(
  filename = "rivers_dist.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)
dat_hfr$rivers_dist <- exp(-3 * rivers_dist)

# routes
routes_dist <- dataverse::get_dataframe_by_name(
  filename = "routes_dist.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)
dat_hfr$routes_dist <- exp(-3 * routes_dist)

# settles
settle_dist <- dataverse::get_dataframe_by_name(
  filename = "settle_dist.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)
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


iraq_district <- dataverse::get_dataframe_by_name(
  filename = "iraq_dist_shapefile.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)
iraq_district <- sf::st_as_sf(iraq_district)



aid_district <- dataverse::get_dataframe_by_name(
  filename = "aid_district.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)




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
ethnicity <- dataverse::get_dataframe_by_name(
  filename = "ethnicity.rds", dataset = "doi:10.7910/DVN/IU8RQK",
  server = "dataverse.harvard.edu", .f = readRDS, original = TRUE)


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
truncation <- 0.95
for (al in 1:length(aid_lag)) {
  cat("Start al: ",al, "\n")

  
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")


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
      plot(res,result = "beta",scale = 81)
     
      # save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
    }
  }
  

}










#-------------------------Continuous aid--------------------------
nbaseE <-5
truncation <- 0.98
firsttime <- TRUE
for (al in 1:length(aid_lag)) {
  cat("Start al: ",al, "\n")

  if(firsttime){
    get_aid_percapita <- function(t,aid,pop,dimyx=c(128,128)){
      res <- log(as.im(exp(aid[[t]])-1,dimyx = dimyx)/as.im(exp(pop[[t]])-1,dimyx=dimyx)+1)
      return(res)
    }
    dat_hfr$aid_per_cap <- lapply(dat_hfr$time, function(x) get_aid_percapita(x,aid = dat_hfr$aid, pop=dat_hfr$pop))
  }
  
  
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  
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
      plot.cate(res,categorical = c(1,rep(0,20)),maineffect = TRUE,scale = 81)
      save(res, file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
    }
  }
  
  
}

library(tidyverse)

# Calculate density of x values
density_x <- density(E.mat,na.rm=TRUE)
data <- data.frame(x=c(E.mat))
data <- data %>%filter(!is.na(x))
data <- data %>% filter(x>tmp_points[1] & x < tmp_points[2] | x==0)
# Scale the density to fit under the x-axis
density_x$y <- density_x$y / max(density_x$y) * 0.05 # Scale to 5% of y-axis range

# Add polygon to plot
p+ geom_hdr_rug(aes( x = x),data = data,show.legend = FALSE)






truncation <- 0.98
for (al in 1:length(aid_lag)) {
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  # make df for plot
  df_plot <- data.frame(array(0,c(0,8)))
  colnames(df_plot) <- c("x","m","ub","lb","M","Outcome","Meaneff","pvalue")
  for (oo in c("SAF","IED")) {
    for (mm in c(3,7,10)) {
      scale <- 81
      load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_percapita_",truncation*100,".dat"))
      m <- res$est_eval *scale
      V <- diag(res$V_eval)*scale^2
      x <- res$specification$eval_values
      ub <- m+qnorm(c(0.975))*sqrt(V)
      lb <- m+qnorm(c(0.025))*sqrt(V)
      meaneffect <- res$mean_effect*scale
      df_tmp <- cbind.data.frame(x = x,m,ub,lb,M = mm,Outcome = oo,Meaneff = meaneffect,pvalue = res$p.value)
      df_plot <- rbind.data.frame(df_plot,df_tmp)
      
      
    }
    
  }
  
  Mlist <- c(3,7,10)

  print(max(df_plot$x))
  print(min(df_plot$x))
  df_plot <- df_plot[df_plot$M%in%Mlist,]
  df_plot$M <- as.factor(df_plot$M)
  df_plot$x <- exp(df_plot$x)-1
  levels(df_plot$M) <-  paste0("M=",Mlist)
  
  a <- df_plot[,c(5,6,8)]
  a <- a %>% distinct()
  
  p <- ggplot() +
    geom_point(data = subset(df_plot, df_plot$x==0), aes(x = x, y = m), color = "black") +
    geom_errorbar(data = subset(df_plot, df_plot$x==0), aes(x = x, ymin = lb, ymax = ub), width = 0.1) +
    geom_line(data = subset(df_plot, df_plot$x!=0), aes(x = x, y = m), color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_hline(data = df_plot, aes(yintercept = Meaneff), linetype = "dashed", color = "#D55E00")+
    geom_ribbon(data = subset(df_plot, df_plot$x!=0), aes(x = x, ymin = lb, ymax = ub), fill = "grey", alpha = 0.5) +
    facet_grid(Outcome~M, scales = "free_y") +theme_bw()+
    labs(x = "Aid per capita", y = "CATE (District Effect)")+
    theme(
      axis.text.x = element_text(size = 12),  # Increase size of x-axis tick labels
      axis.text.y = element_text(size = 12)   # Increase size of y-axis tick labels
    ) + geom_text(
      size    = 4.5,
      data    = a,
      mapping = aes(x = Inf, y = Inf, label = sprintf("p-value=%.3f", pvalue)),
      hjust   = 1.05,
      vjust   = 10.8
    )
  
  
  

  
  print(p)
  ggsave(file = paste0(save_path,"plot_continuous_",aid_lag_name[al],"_",truncation*100,".pdf"), plot = p, width = 7, height = 3.8, units = "in",device = "pdf")
  # ggsave(file = paste0(save_path,"plot_continuous_poster",aid_lag_name[al],".png"), plot = p, width = 6, height = 5, units = "in", dpi = 300)
}



#-------------------Binary aid analysis+intensity+voilence------------------------------
for (al in 1:length(aid_lag)) {
  cat("Start al: ",al, "\n")

  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  
  load("~/phi_im.dat")
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      if(al==2){
        em <- get_em_vec(dat_hfr$aid[31:499],time_after = TRUE,lag=mm,entire_window = NULL,dimyx = c(128,128))
        em <- em>0
        E.basis <- matrix(em,ncol = 1)
      }else{
        E.mat <- aid_st[,,31:(498-mm+1)]
        E.mat[E.mat>0] <- 1
        E.mat[(E.mat==0)] <- 0
        E.basis <- matrix(c(E.mat),ncol = 1)
      }
      
      violence <- get_em_vec(dat_hfr$hist_allout_30[31:499],time_after = TRUE,lag=mm,entire_window = NULL,dimyx = c(128,128))
      vmedian <- median(violence,na.rm=TRUE)
      inten <- rep(c(as.matrix(phi_im6-phi_im1)),length(31:(498-mm+1)))
      inten <- sqrt(inten)
      imedian <- median(unique(c(as.matrix(phi_im6-phi_im1))),na.rm=T)
      imedian <- sqrt(imedian)
      E.basis <- cbind(E.basis,inten,E.basis*inten,violence,E.basis*violence)
      
      points <- c(1)
      
      #points.basis <- cbind(1,points)
      points.basis <- cbind(points,imedian,imedian,vmedian,vmedian)
      
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
      save(res,est_diff,est_diff_var,file = paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten_violence.dat"))

        
      
    }
    
  }
  
  
}





# plots
for (al in 1:3) {
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  # make df for plot
  df_plot <- data.frame(array(0,c(0,5)))
  colnames(df_plot) <- c("m","ub","lb","M","Outcome")
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      scale <- 81
      load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten_violence.dat"))
      
      est_diff <- points.basis%*%res$est_beta
      est_diff_var <- points.basis%*%res$V_beta%*%t(points.basis)
      
      m <- est_diff*scale
      V <- est_diff_var*scale^2
      
      ub <- m+qnorm(c(0.975))*sqrt(V)
      lb <- m+qnorm(c(0.025))*sqrt(V)
      df_tmp <- cbind.data.frame(m,ub,lb,M = mm,Outcome = oo)
      
      
      
      df_plot <- rbind.data.frame(df_plot,df_tmp)
      
      
    }
    
  }
  
  
  
  df_plot$M <- as.factor(df_plot$M)
  levels(df_plot$M) <-  paste0("M=",1:10)
  
  p <- ggplot(df_plot, aes(x = M, y = m)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red")+
    facet_grid(Outcome ~ ., scales = "free_y") +
    labs(x = "", y = "Expected difference") +
    theme_bw()+
    theme(legend.position = "bottom")
  p
  
  ggsave(file = paste0(save_path,"plot_continuous_inten_voilence",aid_lag_name[al],".png"), plot = p, width = 6, height = 3.5, units = "in", dpi = 300)
}





#-------------------Binary aid analysis+intensity------------------------------
checkfile <- TRUE
for (al in 1:length(aid_lag)) {
  cat("Start al: ",al, "\n")
  
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  
  load("~/phi_im.dat")
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      if(file.exists(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat")) & checkfile){
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


for (al in 1:3) {
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  # make df for plot
  df_plot <- data.frame(array(0,c(0,7)))
  colnames(df_plot) <- c("x","m","ub","lb","M","Outcome","Meaneff")
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      scale <- 81
      load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat"))
      
      m <- res$est_eval *scale
      V <- diag(res$V_eval)*scale^2
      x <- res$specification$eval_values
      ub <- m+qnorm(c(0.975))*sqrt(V)
      lb <- m+qnorm(c(0.025))*sqrt(V)
      meaneffect <- res$mean_effect*scale
      
      
      df_tmp <- cbind.data.frame(x = x,m,ub,lb,M = mm,Outcome = oo,Meaneff = meaneffect)
      df_plot <- rbind.data.frame(df_plot,df_tmp)
      
      
    }
    
  }
  
  
  df_plot$Aid <- as.factor(df_plot$x)
  levels(df_plot$Aid) <- c("Without aid", "With aid")
  df_plot$M <- as.factor(df_plot$M)
  levels(df_plot$M) <-  paste0("M=",1:10)
  
  p <- ggplot(df_plot, aes(x = M, y = m,color = Aid)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    geom_segment(aes(x = as.numeric(M) - 0.3, xend = as.numeric(M) + 0.3, y = Meaneff, yend = Meaneff),
                 inherit.aes = FALSE, color = "#D55E00", linetype = "solid",linewidth=0.6) +
    facet_grid(Outcome ~ ., scales = "free_y") +
    labs(x = "", y = "CATE (District Effect)") +
    theme_bw()+
    theme(
      axis.text.x = element_text(size = 12),  # Increase size of x-axis tick labels
      axis.text.y = element_text(size = 12)   # Increase size of y-axis tick labels
    )+
    theme(legend.position = "top")+
    scale_color_manual(values = c("With aid" = "#2C7A7B", "Without aid" = "#4D4D4D"))
  p
  
  ggsave(file = paste0(save_path,"plot_binary_cate_",aid_lag_name[al],".pdf"), plot = p, width = 6, height = 3.5, units = "in", device = "pdf")
}



# plots
for (al in 1:3) {
  save_path <- paste0("~/Iraq_causal/result/application_final_surge/",aid_lag_name[al],"/")
  
  # make df for plot
  df_plot <- data.frame(array(0,c(0,5)))
  colnames(df_plot) <- c("m","ub","lb","M","Outcome")
  points.basis <- cbind(points,imedian,imedian)
  points.basis <- cbind(0,points.basis)
  for (oo in c("SAF","IED")) {
    for (mm in 1:10) {
      scale <- 81
      load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_inten.dat"))
      imedian <- sqrt(0.36)

      
      est_diff <- points.basis%*%res$est_beta
      est_diff_var <- points.basis%*%res$V_beta%*%t(points.basis)
      
      m <- est_diff*scale
      V <- est_diff_var*scale^2
      
      ub <- m+qnorm(c(0.975))*sqrt(V)
      lb <- m+qnorm(c(0.025))*sqrt(V)
      df_tmp <- cbind.data.frame(m,ub,lb,M = mm,Outcome = oo)
      
      
      
      df_plot <- rbind.data.frame(df_plot,df_tmp)
      
      
    }
    
  }
  
  
  
  df_plot$M <- as.factor(df_plot$M)
  levels(df_plot$M) <-  paste0("M=",1:10)
  
  p <- ggplot(df_plot, aes(x = M, y = m)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = lb, ymax = ub), width = 0.2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
    facet_grid(Outcome ~ ., scales = "free_y") +
    labs(x = "", y = "Expected difference") +
    theme_bw()+
    theme(
      axis.text.x = element_text(size = 12),  # Increase size of x-axis tick labels
      axis.text.y = element_text(size = 12)   # Increase size of y-axis tick labels
    )+
    theme(legend.position = "bottom")
  p
  
  #ggsave(file = paste0(save_path,"plot_continuous_inten",aid_lag_name[al],".png"), plot = p, width = 6, height = 3.5, units = "in", dpi = 300)
}



library(MASS)
x <- inten
b <- boxcox(lm(x ~ 1))
# Exact lambda
lambda <- b$x[which.max(b$y)]
lambda
0.02020202
new_x_exact <- (x ^ lambda - 1) / lambda
hist(new_x_exact)

# Assuming 'x' is your variable to be binned
num_bins <- 5  # Number of bins
bins <- cut(x, num_bins)
custom_breaks <- c(0, 10, 20, 30, 40, 50)  # Custom bin boundaries
bins <- cut(x, breaks = custom_breaks)
bins <- cut(x, num_bins, include.lowest = TRUE)
table(bins)  # Display frequency counts for each bin


q1 <- quantile(inten, 0.25,na.rm=T)
q3 <- quantile(inten, 0.75,na.rm=T)
iqr <- q3 - q1
# Identify outliers
outliers <- inten[inten < q1 - 1.5 * iqr | inten > q3 + 1.5 * iqr]
