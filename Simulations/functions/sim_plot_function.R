draw_estimate_all <- function(df, Effect = "13",name = "s",save = FALSE,save_path = "~/"){
  
  
  
  line_styles <-rep("solid",7)
  markers <- c("circle", "triangle", "square", "triangle", "square")
  colors <- c("black","#e31a1c", "#ff7f00",  "#6a3d9a", "#1f78b4") 
  
  df$type <- factor(df$type)
  newx <- 1:10/10
  df <- df[df$x%in% c(0.01,newx),]
  
  
  if(name=="st"){
    a <- 0.002271924
    b <- 0.796932719
  }else{
    a <- 0.007691263
    b <- 0.558239696
  }
  
  
  for (i in c(2,5)) {
    
    if(i==2){
      which_est <- "IPW"
    }else{
      which_est <- "Hajek"
    }
    
    df <- df[!(df$curve%in%c("IPW with estimated PS (truncated)","Hajek with estimated PS (truncated)")),]
    
    # Plot using ggplot2
    p <-  ggplot(df[df$type==which_est & df$effect%in%Effect & df$time==500,], aes(y = y,x=x,color = curve, linetype = curve, shape = curve)) +
      geom_rect(data = data.frame(xmin = a, xmax = b, ymin = -Inf, ymax = Inf),
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, 
                fill = "gray85", 
                alpha = 0.3)+
      geom_vline(xintercept = a, linetype = "solid", color = "grey80") +
      geom_vline(xintercept = b , linetype = "solid", color = "grey80")+
      geom_line(alpha = 0.7, aes(group = curve)) +
      geom_point(size = 2, alpha = 0.7, aes(shape = curve)) +  # Adjusted size and transparency for points
      scale_linetype_manual(values = line_styles) +
      scale_shape_manual(values = markers[i:(i+2)]) +
      scale_color_manual(values = colors[i:(i+2)]) +
      labs(title = "",
           x = "",
           y = "") +
      theme_bw()+
      theme(legend.position = "top",plot.margin = margin(t = -10, b = -15))+ 
      theme(legend.title=element_blank() # Adjust the text size here
      )+
      theme(axis.text = element_text(size = 11),
            legend.text = element_text(size = 12),
            panel.spacing.x = unit(1.1, "lines"),
            plot.margin = margin(t = 10, r = 20, b = 10, l = -5))+
      facet_grid(stat~M,scales="free")+scale_y_continuous(labels = function(x) paste0("\u2009", sprintf("%.3f", x)))#+
      # geom_line(data = truth[truth$type=="Hajek" & truth$effect%in%Effect & truth$time==500,],
      #           aes(x = x, y = y, group = 1),  # Ensure `group` avoids unintended groupings
      #           inherit.aes = FALSE,           # Prevent inheritance of aesthetics from the main ggplot call
      #           color = "black",               # Use a consistent color for theoretical lines
      #           linetype = "dashed",           # Distinguish theoretical lines
      #           size = 0.8,
      #           show.legend = FALSE)           # Exclude from legend
    
    print(p)
    
    if(save)
      ggsave(filename = paste0(save_path,name,"_",which_est,"est_",Effect,".pdf"), plot = p, width = 7.61, height = 4.61, units = "in", device = "pdf")
  }
  
  
}




draw_sd_ratio <- function(df, M=NULL, Effect = "23", onlyh = TRUE, name="s", save = FALSE, save_path = "~/"){
  
  line_styles <- c("solid", "solid","solid")
  markers <- c("circle", "triangle", "circle", "triangle")
  colors <- c("#e31a1c","#0076B9",  "#6a3d9a", "#ff7f00") 
  
  newx <- 1:10/10
  df <- df[df$x %in% c(0.01, newx),]
  
  df <- df[df$curve == "Estimated SD Bound(Hajek)",]
  df$time <- as.factor(df$time)
  
  low_limit <- min(df$y)
  
  p <- ggplot(df, aes(y = y, x = x, color = time, linetype = time, shape = curve)) +
    geom_line() +  
    geom_point() +  
    scale_linetype_manual(values = line_styles, 
                          labels = c("T=200", "T=500", "T=1000")) +  
    scale_color_manual(values = colors, 
                       labels = c("T=200", "T=500", "T=1000")) +  
    scale_shape_manual(values = markers, guide = "none") +  # Remove black dot legend
    labs(title = "", x = "Moderator", y = "") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank()) +
    facet_grid(M ~ effect, scales = "free", labeller = labeller(effect = as_labeller(cyl_names, label_parsed))) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black")
  
  return(p)
}