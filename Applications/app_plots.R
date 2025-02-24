


















#-------------------------------Figure A.9---------------------------------

load("~/Github/CATE-SpatioTemporal-Causal/Applications/Results/fitted_ps.dat")

plot_names <- c(
  `1` = 'Fitted',
  `2` = 'Residual',
  `3` = 'Fitted (80%)',
  `4` = 'Residual (80%)'
)

plot_data %>% ggplot()+
  geom_line(aes(x = time, y = actual_counts), color = color_actual, linewidth = 0.6)+
  geom_line(aes(x = time, y = fit1), color = color_dens_1, linewidth = 0.6)+
  geom_line(aes(x = time, y = fit2), color = color_dens_2, linewidth = 0.6)+
  geom_line(aes(x = time, y = residuals1), color = color_dens_1, linewidth = 0.6) + 
  geom_line(aes(x = time, y = residuals2), color = color_dens_2, linewidth = 0.6) + 
  geom_hline(aes(yintercept = hline), linetype = "dashed")+
  facet_grid(govern~plot,scales = "free",labeller = labeller(plot = as_labeller(plot_names)))+
  ylab("Count")+
  xlab("Days")+
  theme_bw()


