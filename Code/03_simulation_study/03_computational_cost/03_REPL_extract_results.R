#-------------------# #-------------------# #-------------------#
#                   EXTRACT TIME AND SHOW RESULTS
#-------------------# #-------------------# #-------------------#

library(salso)
library(bSCDCsampler)
library(tidyverse)

#------# #------# #------# #------# #------# #------# #------# #------# #------# 
a_low = 1
par_tau2 = 1
par_sigma2 = 2.5
par1 = 10
par2 = 1


seq_n = c(50,100,200,400) 
seq_T = c(100,200,400)

parameter_comb = expand.grid(x = seq_n, y = seq_T)
parameter_comb
replications_sim = 50


if(!file.exists("03_simulation_study/03_computational_cost/output_RDS/D2.RDS")){
  D <- data.frame()
  for(i in 1:nrow(parameter_comb)){
    
    n <- parameter_comb[i,1]
    TT <- parameter_comb[i,2]
    
    for(n_sim in 1:replications_sim){
      
      name2 <- paste0("03_simulation_study/03_computational_cost/results/time_gibbs_n",n,"_TT", TT, "_sim", n_sim, ".RDS")
      time <- readRDS(name2)
      
      D <- rbind(D, data.frame("n" = n, "T" = TT, repl = n_sim,
                               time =  round(as.numeric(time, units = "secs"),5))) # convert all in seconds
    }
  }
  
  D2 <- D %>% #rename("Number of neurons" = n) %>%
    mutate(Timest = paste0("Number of timestaps: ", T))
  
  D2$Timest[D2$Timest == "Number of timestaps: 100"] = "Number of timestamps: 100"
  D2$Timest[D2$Timest == "Number of timestaps: 200"] = "Number of timestamps: 200"
  D2$Timest[D2$Timest == "Number of timestaps: 400"] = "Number of timestamps: 400"
  
} else {
  D2 = readRDS("03_simulation_study/03_computational_cost/output_RDS/D2.RDS")
}


G = ggplot(D2)+
  theme_bw()+
  geom_boxplot(aes(y=time,x=factor(n),group=factor(n)))+
  facet_wrap(~Timest)+
  xlab("Number of neurons")+
  ylab("Seconds for 500 iterations")+
  # scale_y_log10()+ # optional
  theme(text = element_text(size=13))
G
ggsave("03_simulation_study/03_computational_cost/output_images/comput_time.pdf", G, width = 8, height = 4)
ggsave("03_simulation_study/03_computational_cost/output_images/comput_time.eps",G, width = 8, height = 4, device = cairo_ps)
