load("~/Github/CATE-SpatioTemporal-Causal/Simulations/Results/save_estimate.dat")
source("~/Github/CATE-SpatioTemporal-Causal/Simulations/functions/sim_plot_function.R")
save_path <- "~/GitHub/CATE-SpatioTemporal-Causal/Simulations/Results/"


#------------------------Figure 5, A.2-A.6-----------

draw_estimate_all(merged_df_st, Effect = "13",name = "st",save = FALSE,save_path = save_path)
draw_estimate_all(merged_df_st, Effect = "23",name = "st",save =  FALSE,save_path = save_path)
draw_estimate_all(merged_df_st, Effect = "12",name = "st",save = FALSE,save_path = save_path)


draw_estimate_all(merged_df_s, Effect = "13",name = "s",save = FALSE,save_path = save_path)
draw_estimate_all(merged_df_s, Effect = "23",name = "s",save =  FALSE,save_path = save_path)
draw_estimate_all(merged_df_s, Effect = "12",name = "s",save = FALSE,save_path = save_path)

#------------------------Figure A.1--------------------
draw_sd_ratio(s_df_sd_ratio,name = "s",save = FALSE,save_path = save_path)
draw_sd_ratio(st_df_sd_ratio,name = "st",save = FALSE,save_path = save_path)



#----------------------Figure A.7-----------------------

s_df_sd_relative_ps %>% ggplot(aes(x = x,y=y,color = Estimator))+
  geom_point()+
  scale_shape_manual(values = markers)+
  scale_color_manual(values = default_colors)+
  geom_line()+ 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  facet_grid(effect~M,scales="free",labeller = labeller(effect = as_labeller(cyl_names, label_parsed)))+
  xlab("Moderator")+
  ylab("SD Ratio")+
  theme(
    strip.text.y = element_text(size = 12) # Adjust the text size here
  )+
  theme_bw()+
  theme(legend.position = "top")


st_df_sd_relative_ps %>% ggplot(aes(x = x,y=y,color = Estimator))+
  geom_point()+
  scale_shape_manual(values = markers)+
  scale_color_manual(values = default_colors)+
  geom_line()+ 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  facet_grid(effect~M,scales="free",labeller = labeller(effect = as_labeller(cyl_names, label_parsed)))+
  xlab("Moderator")+
  ylab("SD Ratio")+
  theme(
    strip.text.y = element_text(size = 12) # Adjust the text size here
  )+
  theme_bw()+
  theme(legend.position = "top")



#--------------------------- Figure A.8----------------------

s_df_sd_relative_haj %>% ggplot(aes(x = x,y=y,color = M))+
  geom_point()+
  geom_line()+ 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  facet_wrap(~effect,scales="free",labeller = labeller(effect = as_labeller(cyl_names, label_parsed)))+
  xlab("Moderator")+
  ylab("SD Ratio")+
  theme(
    strip.text.y = element_text(size = 12) # Adjust the text size here
  )+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "", shape = "")


st_df_sd_relative_haj %>% ggplot(aes(x = x,y=y,color = M))+
  geom_point()+
  geom_line()+ 
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  facet_wrap(~effect,scales="free",labeller = labeller(effect = as_labeller(cyl_names, label_parsed)))+
  xlab("Moderator")+
  ylab("SD Ratio")+
  theme(
    strip.text.y = element_text(size = 12) # Adjust the text size here
  )+
  theme_bw()+
  theme(legend.position = "top")+
  labs(color = "", shape = "")
