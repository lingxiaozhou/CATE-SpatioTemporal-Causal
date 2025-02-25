
library(ggplot2)
library(tidyverse)

outpath <- "~/GitHub/CATE-SpatioTemporal-Causal/Applications/Results/"
save_plot <- FALSE


#--------------------------Figure 6, A.14, A.15--------------------

truncation_vec <- c(0.9,0.95,0.98)
for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    # make df for plot
    df_plot <- data.frame(array(0,c(0,7)))
    colnames(df_plot) <- c("x","m","ub","lb","M","Outcome","Meaneff")
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        scale <- 81
        load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        
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
      scale_color_manual(values = c("With aid" = "#2C7A7B", "Without aid" = "#F4A261"))
    print(p)
    
    if(save_plot){
      ggsave(file = paste0(save_path,"plot_binary_cate_",aid_lag_name[al],"_",truncation*100,".pdf"), plot = p, width = 6, height = 3.5, units = "in", device = "pdf")
    }
    
  }
}



for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
    # make df for plot
    df_plot <- data.frame(array(0,c(0,5)))
    colnames(df_plot) <- c("m","ub","lb","M","Outcome")
    for (oo in c("SAF","IED")) {
      for (mm in 1:10) {
        scale <- 81
        load(paste0(save_path, "res_Type1_M", mm,"_",oo,"_trun_aid_bi_",truncation*100,".dat"))
        m <- res$est_beta[2]*scale
        V <- diag(res$V_beta)[2]*scale^2
        
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
    
    print(p)
    if(save_plot){
      ggsave(file = paste0(save_path,"plot_binary_beta_",aid_lag_name[al],"_",truncation*100,".pdf"), plot = p, width = 6, height = 3, units = "in", device = "pdf")
    }
    
  }
  
}



#-------------------------Figure 8, A.16, A.17-------------------------

for (truncation in truncation_vec) {
  for (al in 1:length(aid_lag)) {
    save_path <- paste0(outpath,aid_lag_name[al],"/")
    
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
    
    if(save_plot){
      ggsave(file = paste0(save_path,"plot_continuous_",aid_lag_name[al],"_",truncation*100,".pdf"), plot = p, width = 7, height = 3.8, units = "in",device = "pdf")
    }
    
    
  }
}






for (al in 1:length(aid_lag)) {
  save_path <- paste0(outpath,aid_lag_name[al],"/")
  
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
  print(p)
  if(save_plot){
    ggsave(file = paste0(save_path,"plot_binary_cate_",aid_lag_name[al],".pdf"), plot = p, width = 6, height = 3.5, units = "in", device = "pdf")
  }
  
}



for (al in 1:length(aid_lag)) {
  save_path <- paste0(outpath,aid_lag_name[al],"/")
  
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
  print(p)
  if(save_plot){
    ggsave(file = paste0(save_path,"plot_continuous_inten",aid_lag_name[al],".png"), plot = p, width = 6, height = 3.5, units = "in", dpi = 300)
  }
  
}


